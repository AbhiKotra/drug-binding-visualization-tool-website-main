#!/usr/bin/env python3
"""
blast_search.py
Queries NCBI BLAST (blastp vs PDB database) to find structurally similar
proteins from other species (>= identity_threshold % identity).

Outputs:
  blast_results.csv    — tabular hit data
  blast_summary.txt    — human-readable report

Usage: python blast_search.py <bound.pdb> <results_dir> [identity_threshold]
"""

from Bio.PDB import PDBParser, PPBuilder
from Bio.Blast import NCBIWWW, NCBIXML
import os
import sys
import csv

THREE_TO_ONE = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C',
    'GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I',
    'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P',
    'SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V',
    'SEC':'U','PYL':'O','MSE':'M','UNK':'X',
}


def extract_sequence(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('prot', pdb_file)
    ppb = PPBuilder()
    seq = ''
    for pp in ppb.build_peptides(structure):
        for res in pp:
            seq += THREE_TO_ONE.get(res.get_resname().strip(), 'X')
        if seq:
            break  # use first polypeptide chain only
    return seq


def run_blast(pdb_file, results_dir, identity_threshold=80):
    os.makedirs(results_dir, exist_ok=True)

    seq = extract_sequence(pdb_file)
    if len(seq) < 20:
        _write_empty(results_dir, 'Sequence too short for BLAST.')
        return []

    # Cap at 500 residues for speed; NCBI handles longer but it's slower
    blast_seq = seq[:500]

    print(f"BLAST: querying NCBI blastp/pdb with {len(blast_seq)}-residue sequence...", flush=True)

    try:
        result_handle = NCBIWWW.qblast(
            'blastp', 'pdb', blast_seq,
            hitlist_size=20,
            expect=0.001,
        )
    except Exception as e:
        _write_empty(results_dir, f'NCBI BLAST request failed: {e}')
        return []

    try:
        blast_records = list(NCBIXML.parse(result_handle))
    except Exception as e:
        _write_empty(results_dir, f'BLAST XML parse error: {e}')
        return []

    hits = []
    for record in blast_records:
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                if hsp.align_length == 0:
                    continue
                identity_pct = round(100.0 * hsp.identities / hsp.align_length, 1)
                if identity_pct < identity_threshold:
                    continue

                title = alignment.title or ''
                # Title format: "1ABC_A Chain A, Protein Name [Organism]"
                pdb_chain = title.split()[0] if title else ''
                pdb_id  = pdb_chain[:4].upper()  if len(pdb_chain) >= 4 else ''
                chain_id = pdb_chain[5:6].upper() if len(pdb_chain) >= 6 else 'A'

                organism = ''
                if '[' in title and ']' in title:
                    organism = title[title.rfind('[')+1 : title.rfind(']')]

                hits.append({
                    'pdb_id':            pdb_id,
                    'chain':             chain_id,
                    'organism':          organism,
                    'identity_pct':      identity_pct,
                    'similarity_pct':    round(100.0 * hsp.positives / hsp.align_length, 1),
                    'e_value':           f'{hsp.expect:.2e}',
                    'alignment_length':  hsp.align_length,
                    'identities':        hsp.identities,
                    'query_start':       hsp.query_start,
                    'query_end':         hsp.query_end,
                    'query_seq_preview': str(hsp.query)[:60],
                    'hit_seq_preview':   str(hsp.sbjct)[:60],
                    'description':       title[:100],
                })
                break  # best HSP per alignment only

    # ── Write CSV ──
    csv_path = os.path.join(results_dir, 'blast_results.csv')
    if hits:
        with open(csv_path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=list(hits[0].keys()))
            writer.writeheader()
            writer.writerows(hits)
    else:
        _write_empty(results_dir, f'No hits >= {identity_threshold}% identity found.')
        print('BLAST: no qualifying hits.')
        return []

    # ── Write text summary ──
    txt_path = os.path.join(results_dir, 'blast_summary.txt')
    with open(txt_path, 'w', encoding='utf-8') as f:
        f.write("NCBI BLAST Results — blastp vs PDB Database\n")
        f.write("=" * 55 + "\n\n")
        f.write(f"Query sequence length: {len(seq)} residues (first {len(blast_seq)} used)\n")
        f.write(f"Identity threshold:    >= {identity_threshold}%\n")
        f.write(f"Qualifying hits:       {len(hits)}\n\n")
        for h in hits:
            f.write(f"PDB {h['pdb_id']} Chain {h['chain']}\n")
            f.write(f"  Identity:   {h['identity_pct']}%\n")
            f.write(f"  Similarity: {h['similarity_pct']}%\n")
            f.write(f"  E-value:    {h['e_value']}\n")
            if h['organism']:
                f.write(f"  Organism:   {h['organism']}\n")
            f.write(f"  Query:  {h['query_seq_preview']}\n")
            f.write(f"  Hit:    {h['hit_seq_preview']}\n\n")

    print(f"BLAST complete: {len(hits)} hit(s) >= {identity_threshold}% identity")
    return hits


def _write_empty(results_dir, reason):
    csv_path = os.path.join(results_dir, 'blast_results.csv')
    txt_path = os.path.join(results_dir, 'blast_summary.txt')
    with open(csv_path, 'w', encoding='utf-8') as f:
        f.write('pdb_id,chain,organism,identity_pct,similarity_pct,e_value\n')
    with open(txt_path, 'w', encoding='utf-8') as f:
        f.write(f"BLAST: {reason}\n")


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python blast_search.py <bound.pdb> <results_dir> [identity_threshold]")
        sys.exit(1)
    thresh = int(sys.argv[3]) if len(sys.argv) > 3 else 80
    run_blast(sys.argv[1], sys.argv[2], thresh)
