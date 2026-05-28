#!/usr/bin/env python3
"""
analyze_occupancy.py
Parses PDB occupancy and B-factor (thermal ellipsoid) parameters.

Occupancy = 1.0  means the atom is in a single well-defined position.
Occupancy < 1.0  (e.g. 0.5) means two alternate conformations exist;
                 the high B-factor that often accompanies these indicates
                 high thermal motion / weak binding at that site.

Usage: python analyze_occupancy.py <bound.pdb> <ligand_code> <results_dir>
"""

from Bio.PDB import PDBParser
import numpy as np
from collections import defaultdict
import os
import sys
import csv


def analyze_occupancy(pdb_file, ligand_code, results_dir):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)

    all_bfactors = []
    residue_bfactors = defaultdict(list)
    low_occ_residues = {}
    ligand_atoms_data = []

    for model in structure:
        for chain in model:
            for residue in chain:
                hetflag  = residue.id[0]
                resname  = residue.get_resname().strip()
                resnum   = residue.id[1]
                chain_id = chain.id
                is_ligand  = (hetflag != ' ' and resname.upper() == ligand_code.upper())
                is_protein = (hetflag == ' ')

                for atom in residue:
                    occ  = atom.get_occupancy()
                    bfac = atom.get_bfactor()
                    aname = atom.get_name().strip()

                    if bfac is not None and bfac > 0:
                        all_bfactors.append(bfac)
                        key = (chain_id, resname, resnum)
                        residue_bfactors[key].append(bfac)

                    if occ is not None and occ < 1.0:
                        key = (chain_id, resname, resnum)
                        if key not in low_occ_residues:
                            low_occ_residues[key] = {
                                'chain': chain_id, 'residue': resname,
                                'resnum': resnum, 'occupancy': occ,
                                'atoms': [], 'max_bfactor': 0.0,
                                'type': 'ligand' if is_ligand else ('protein' if is_protein else 'hetatm'),
                            }
                        low_occ_residues[key]['atoms'].append(aname)
                        if bfac:
                            low_occ_residues[key]['max_bfactor'] = max(
                                low_occ_residues[key]['max_bfactor'], bfac)

                    if is_ligand:
                        ligand_atoms_data.append({
                            'atom': aname, 'occupancy': occ,
                            'b_factor': round(bfac, 2) if bfac else 0.0,
                        })

    # B-factor statistics
    if all_bfactors:
        mean_b = float(np.mean(all_bfactors))
        std_b  = float(np.std(all_bfactors))
        high_thresh = mean_b + 2.0 * std_b
    else:
        mean_b = std_b = high_thresh = 0.0

    # High-B-factor residues
    high_bfac_list = []
    for key, bfacs in residue_bfactors.items():
        avg = float(np.mean(bfacs))
        if avg > high_thresh:
            high_bfac_list.append({
                'chain': key[0], 'residue': key[1], 'resnum': key[2],
                'avg_bfactor': round(avg, 2),
                'max_bfactor': round(max(bfacs), 2),
                'threshold': round(high_thresh, 2),
            })
    high_bfac_list.sort(key=lambda x: -x['avg_bfactor'])

    os.makedirs(results_dir, exist_ok=True)

    # ── Low-occupancy CSV ──
    low_occ_csv = os.path.join(results_dir, 'low_occupancy_residues.csv')
    with open(low_occ_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'chain', 'residue', 'resnum', 'occupancy', 'max_bfactor', 'atoms', 'type'])
        writer.writeheader()
        for data in low_occ_residues.values():
            writer.writerow({
                'chain': data['chain'], 'residue': data['residue'],
                'resnum': data['resnum'], 'occupancy': data['occupancy'],
                'max_bfactor': round(data['max_bfactor'], 2),
                'atoms': '; '.join(data['atoms']), 'type': data['type'],
            })

    # ── High-B-factor CSV ──
    high_b_csv = os.path.join(results_dir, 'high_bfactor_residues.csv')
    with open(high_b_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'chain', 'residue', 'resnum', 'avg_bfactor', 'max_bfactor', 'threshold'])
        writer.writeheader()
        writer.writerows(high_bfac_list)

    # ── Text summary ──
    txt_path = os.path.join(results_dir, 'occupancy_bfactor_summary.txt')
    with open(txt_path, 'w', encoding='utf-8') as f:
        f.write("Occupancy & B-Factor (Thermal Ellipsoid) Analysis\n")
        f.write("=" * 55 + "\n\n")
        f.write(f"Overall B-factor statistics:\n")
        f.write(f"  Mean B-factor:   {mean_b:.2f} A^2\n")
        f.write(f"  Std B-factor:    {std_b:.2f} A^2\n")
        f.write(f"  High threshold (mean + 2*std): {high_thresh:.2f} A^2\n\n")

        f.write(f"Residues with occupancy < 1.0: {len(low_occ_residues)}\n")
        f.write("(Indicates dual conformations — atom exists in two alternative positions)\n\n")
        for data in sorted(low_occ_residues.values(), key=lambda x: (x['chain'], x['resnum'])):
            f.write(f"  Chain {data['chain']}: {data['residue']} {data['resnum']}"
                    f"  occ={data['occupancy']}  max_B={data['max_bfactor']:.1f} A^2"
                    f"  [{data['type']}]\n")
            f.write(f"    Affected atoms: {', '.join(data['atoms'])}\n")

        f.write(f"\nResidues with high B-factor (> {high_thresh:.1f} A^2): {len(high_bfac_list)}\n")
        f.write("(High B-factor = high thermal motion / structural flexibility / weak binding)\n\n")
        for res in high_bfac_list[:20]:
            f.write(f"  Chain {res['chain']}: {res['residue']} {res['resnum']}"
                    f"  avg_B={res['avg_bfactor']} A^2  max_B={res['max_bfactor']} A^2\n")

        if ligand_atoms_data:
            lig_bfacs = [a['b_factor'] for a in ligand_atoms_data if a['b_factor'] > 0]
            lig_occs  = [a['occupancy'] for a in ligand_atoms_data if a['occupancy'] is not None]
            f.write(f"\nLigand ({ligand_code}) analysis:\n")
            if lig_bfacs:
                lig_avg_b = float(np.mean(lig_bfacs))
                f.write(f"  Mean B-factor: {lig_avg_b:.2f} A^2\n")
                f.write(f"  Max  B-factor: {max(lig_bfacs):.2f} A^2\n")
                if lig_avg_b > high_thresh:
                    f.write(f"  WARNING: Ligand B-factor > threshold — possible weak/flexible binding\n")
                else:
                    f.write(f"  Ligand B-factor within normal range — stable binding conformation\n")
            if lig_occs:
                min_occ = min(lig_occs)
                if min_occ < 1.0:
                    f.write(f"  NOTE: Ligand occupancy = {min_occ} — multiple binding conformations\n")
                else:
                    f.write(f"  Occupancy = 1.0 (fully occupied binding site)\n")

    print(f"Occupancy/B-factor analysis complete")
    print(f"  Low-occupancy residues: {len(low_occ_residues)}")
    print(f"  High-B-factor residues: {len(high_bfac_list)}")
    return {
        'low_occupancy': list(low_occ_residues.values()),
        'high_bfactor': high_bfac_list,
        'mean_bfactor': round(mean_b, 2),
        'high_threshold': round(high_thresh, 2),
    }


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print("Usage: python analyze_occupancy.py <bound.pdb> <ligand_code> <results_dir>")
        sys.exit(1)
    analyze_occupancy(sys.argv[1], sys.argv[2], sys.argv[3])
