#!/usr/bin/env python3
"""
sequence_comparison.py
Extracts sequences from both PDB files, runs pairwise global alignment,
and generates a comparison table including:
  - Sequence identity / similarity %
  - Whether sequences are >=90% similar (same/orthologous protein)
  - Binding pocket residue mapping between structures
  - Summary table suitable for poster display

Usage: python sequence_comparison.py <alphafold.pdb> <bound.pdb> <results_dir> [ligand_code]
"""

from Bio.PDB import PDBParser, PPBuilder
from Bio.Align import PairwiseAligner
import os
import sys
import csv
import math

THREE_TO_ONE = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C',
    'GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I',
    'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P',
    'SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V',
    'SEC':'U','PYL':'O','MSE':'M','UNK':'X',
}

CONSERVATIVE_GROUPS = [
    {'I','L','V','M'},
    {'F','Y','W'},
    {'K','R','H'},
    {'D','E'},
    {'S','T'},
    {'N','Q'},
    {'G','A'},
]


def conservative_similar(a, b):
    for g in CONSERVATIVE_GROUPS:
        if a in g and b in g:
            return True
    return False


def extract_sequences(pdb_file):
    """Return {chain_id: {'seq': str, 'residues': list of (resnum, resname, aa)}}"""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('prot', pdb_file)
    ppb = PPBuilder()
    chains = {}
    for pp in ppb.build_peptides(structure):
        chain_id = pp[0].get_parent().id
        if chain_id not in chains:
            chains[chain_id] = {'seq': '', 'residues': []}
        for res in pp:
            aa = THREE_TO_ONE.get(res.get_resname().strip(), 'X')
            chains[chain_id]['seq'] += aa
            chains[chain_id]['residues'].append(
                (res.id[1], res.get_resname().strip(), aa))
    return chains


def pairwise_align(seq1, seq2):
    """Global alignment using PairwiseAligner; returns (aligned1, aligned2, identity, similarity)."""
    if not seq1 or not seq2:
        return '', '', 0.0, 0.0

    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score    =  2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score    = -5.0
    aligner.extend_gap_score  = -0.5

    alignments = aligner.align(seq1, seq2)
    # Take the top alignment
    try:
        best = next(iter(alignments))
    except StopIteration:
        return seq1, seq2, 0.0, 0.0

    aligned1 = str(best).split('\n')[0]
    aligned2 = str(best).split('\n')[2] if len(str(best).split('\n')) > 2 else seq2

    # Recalculate identity from aligned strings
    identical = similar = aligned_pos = 0
    for a, b in zip(aligned1, aligned2):
        if a == '-' or b == '-':
            continue
        aligned_pos += 1
        if a == b:
            identical += 1
            similar += 1
        elif conservative_similar(a, b):
            similar += 1

    if aligned_pos == 0:
        return aligned1, aligned2, 0.0, 0.0

    identity   = round(100.0 * identical / aligned_pos, 1)
    similarity = round(100.0 * similar   / aligned_pos, 1)
    return aligned1, aligned2, identity, similarity


def run_sequence_comparison(af_file, bound_file, results_dir, ligand_code=''):
    af_chains    = extract_sequences(af_file)
    bound_chains = extract_sequences(bound_file)

    os.makedirs(results_dir, exist_ok=True)

    all_chain_ids = sorted(set(list(af_chains.keys()) + list(bound_chains.keys())))
    comparison_rows = []
    summary_lines   = []

    summary_lines += [
        "Sequence Comparison: AlphaFold vs Drug-Bound Structure",
        "=" * 60,
        f"AlphaFold file:    {os.path.basename(af_file)}",
        f"Drug-bound file:   {os.path.basename(bound_file)}",
        f"Ligand code:       {ligand_code}",
        "",
    ]

    for chain_id in all_chain_ids:
        af_entry    = af_chains.get(chain_id,    {'seq': '', 'residues': []})
        bound_entry = bound_chains.get(chain_id, {'seq': '', 'residues': []})
        seq_af    = af_entry['seq']
        seq_bound = bound_entry['seq']

        if not seq_af and not seq_bound:
            continue

        aligned_af, aligned_bound, identity, similarity = pairwise_align(seq_af, seq_bound)

        row = {
            'chain':                  chain_id,
            'af_length':              len(seq_af),
            'bound_length':           len(seq_bound),
            'sequence_identity_pct':  identity,
            'sequence_similarity_pct': similarity,
            'high_similarity_ge90':   identity >= 90.0,
            'aligned_af_preview':     aligned_af[:80] + ('...' if len(aligned_af) > 80 else ''),
            'aligned_bound_preview':  aligned_bound[:80] + ('...' if len(aligned_bound) > 80 else ''),
        }
        comparison_rows.append(row)

        summary_lines += [
            f"Chain {chain_id}",
            f"  AlphaFold length  : {len(seq_af)} residues",
            f"  Drug-bound length : {len(seq_bound)} residues",
            f"  Sequence identity : {identity}%",
            f"  Sequence similarity (conservative subs): {similarity}%",
        ]
        if identity >= 90:
            summary_lines.append(f"  STATUS: >=90% identity — same or orthologue protein; "
                                  f"superposition deviation (RMSD) shown in Advanced Metrics")
        elif identity >= 50:
            summary_lines.append(f"  STATUS: Moderate identity ({identity}%) — related proteins")
        else:
            summary_lines.append(f"  STATUS: Low identity — may be divergent proteins")
        summary_lines.append(f"  Alignment preview (first 80 chars):")
        summary_lines.append(f"    AF:    {aligned_af[:80]}")
        summary_lines.append(f"    Bound: {aligned_bound[:80]}")
        summary_lines.append("")

    # ── Binding pocket residue correspondence table ──
    # Show which residues in AF map to residues in bound structure
    pocket_rows = []
    for chain_id in all_chain_ids:
        af_res    = {r[0]: r for r in af_chains.get(chain_id, {'residues': []})['residues']}
        bound_res = {r[0]: r for r in bound_chains.get(chain_id, {'residues': []})['residues']}
        shared_resnums = sorted(set(af_res.keys()) & set(bound_res.keys()))
        for rn in shared_resnums:
            af_aa    = af_res[rn][2]
            bound_aa = bound_res[rn][2]
            pocket_rows.append({
                'chain': chain_id,
                'resnum': rn,
                'af_resname': af_res[rn][1],
                'af_aa': af_aa,
                'bound_resname': bound_res[rn][1],
                'bound_aa': bound_aa,
                'conserved': af_aa == bound_aa,
                'conservative_sub': conservative_similar(af_aa, bound_aa) if af_aa != bound_aa else False,
            })

    # ── Write CSVs ──
    seq_csv = os.path.join(results_dir, 'sequence_comparison.csv')
    if comparison_rows:
        with open(seq_csv, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=list(comparison_rows[0].keys()))
            writer.writeheader()
            writer.writerows(comparison_rows)

    pocket_csv = os.path.join(results_dir, 'residue_correspondence.csv')
    if pocket_rows:
        with open(pocket_csv, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=list(pocket_rows[0].keys()))
            writer.writeheader()
            writer.writerows(pocket_rows)

    # ── Write text summary ──
    txt_path = os.path.join(results_dir, 'sequence_comparison_summary.txt')
    with open(txt_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(summary_lines))
        if pocket_rows:
            f.write("\nResidue correspondence (first 30 shared residues):\n")
            f.write(f"  {'Chain':5s} {'Res#':6s} {'AF aa':8s} {'Bound aa':10s} {'Status':20s}\n")
            f.write("  " + "-" * 55 + "\n")
            for pr in pocket_rows[:30]:
                if pr['conserved']:
                    status = 'Conserved'
                elif pr['conservative_sub']:
                    status = 'Conservative sub'
                else:
                    status = 'Mutation'
                f.write(f"  {pr['chain']:5s} {pr['resnum']:6d} "
                        f"{pr['af_resname']:8s} {pr['bound_resname']:10s} {status:20s}\n")

    print(f"Sequence comparison complete: {len(comparison_rows)} chain(s) analysed")
    for row in comparison_rows:
        print(f"  Chain {row['chain']}: {row['sequence_identity_pct']}% identity, "
              f"{row['sequence_similarity_pct']}% similarity")

    return comparison_rows


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print("Usage: python sequence_comparison.py <alphafold.pdb> <bound.pdb> <results_dir> [ligand_code]")
        sys.exit(1)
    lig = sys.argv[4] if len(sys.argv) > 4 else ''
    run_sequence_comparison(sys.argv[1], sys.argv[2], sys.argv[3], lig)
