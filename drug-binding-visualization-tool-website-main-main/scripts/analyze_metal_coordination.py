#!/usr/bin/env python3
"""
analyze_metal_coordination.py
Detects metal ions and analyzes coordination bonds.

Handles Fe (heme/porphyrin: 4 equatorial N + axial ligands),
Zn, Mg, Ca, Cu, Mn, Co, Ni, Na, K, Pt, Hg, Cd.

Usage: python analyze_metal_coordination.py <bound.pdb> <ligand_code> <results_dir>
"""

from Bio.PDB import PDBParser
import numpy as np
import os
import sys
import csv
import math

# Metal element symbols and their typical max coordination distances (A)
METALS = {
    'FE':  {'max_dist': 2.4, 'name': 'Iron'},
    'FE2': {'max_dist': 2.4, 'name': 'Iron(II)'},
    'FE3': {'max_dist': 2.2, 'name': 'Iron(III)'},
    'ZN':  {'max_dist': 2.5, 'name': 'Zinc'},
    'MG':  {'max_dist': 2.5, 'name': 'Magnesium'},
    'CA':  {'max_dist': 2.8, 'name': 'Calcium'},
    'CU':  {'max_dist': 2.4, 'name': 'Copper'},
    'CU1': {'max_dist': 2.4, 'name': 'Copper(I)'},
    'MN':  {'max_dist': 2.5, 'name': 'Manganese'},
    'CO':  {'max_dist': 2.3, 'name': 'Cobalt'},
    'NI':  {'max_dist': 2.3, 'name': 'Nickel'},
    'NA':  {'max_dist': 2.8, 'name': 'Sodium'},
    'K':   {'max_dist': 3.0, 'name': 'Potassium'},
    'PT':  {'max_dist': 2.5, 'name': 'Platinum'},
    'HG':  {'max_dist': 2.6, 'name': 'Mercury'},
    'CD':  {'max_dist': 2.5, 'name': 'Cadmium'},
}

HEME_RESIDUES = {'HEM', 'HEC', 'MHM', 'HEB', 'HEA', 'SRM', 'CLF', 'BCL', 'MG3', 'CLA'}

COORD_ELEMENTS = {'O', 'N', 'S', 'CL', 'F', 'SE'}


def get_element(atom):
    elem = atom.element
    if elem:
        return elem.strip().upper()
    name = atom.get_name().strip()
    return name[0].upper() if name else ''


def classify_geometry(n):
    return {
        1: 'Single ligand',
        2: 'Linear (2-coordinate)',
        3: 'Trigonal planar',
        4: 'Tetrahedral / Square planar',
        5: 'Trigonal bipyramidal / Square pyramidal',
        6: 'Octahedral',
        7: 'Pentagonal bipyramidal',
        8: 'Cubic / Square antiprismatic',
    }.get(n, f'{n}-coordinate (unusual)')


def analyze_metal_coordination(bound_file, ligand_code, results_dir):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('complex', bound_file)

    # Collect all atoms and residues once
    all_residues = []
    for model in structure:
        for chain in model:
            for residue in chain:
                all_residues.append((residue, chain))

    # Find metal atoms
    metal_instances = []
    for res, chain in all_residues:
        for atom in res:
            elem = get_element(atom)
            aname = atom.get_name().strip().upper()
            metal_key = None
            if elem in METALS:
                metal_key = elem
            elif aname in METALS:
                metal_key = aname
            if metal_key:
                metal_instances.append({
                    'atom': atom, 'residue': res, 'chain': chain,
                    'metal_key': metal_key,
                    'metal_name': METALS[metal_key]['name'],
                    'max_dist': METALS[metal_key]['max_dist'],
                })

    results = []

    for mi in metal_instances:
        metal_atom = mi['atom']
        max_dist    = mi['max_dist']

        coord_atoms = []
        for res, chain in all_residues:
            for atom in res:
                if atom is metal_atom:
                    continue
                elem = get_element(atom)
                if elem not in COORD_ELEMENTS:
                    continue
                dist = metal_atom - atom
                if dist > max_dist:
                    continue

                is_lig   = (res.id[0] != ' ' and
                            res.get_resname().strip().upper() == ligand_code.upper())
                is_prot  = res.id[0] == ' '
                is_heme_ring = res.get_resname().strip().upper() in HEME_RESIDUES
                atype = ('protein' if is_prot
                         else 'ligand' if is_lig
                         else 'heme/cofactor' if is_heme_ring
                         else 'water/hetatm')

                coord_atoms.append({
                    'atom_name': atom.get_name().strip(),
                    'element': elem,
                    'residue': res.get_resname().strip(),
                    'resnum': res.id[1],
                    'chain': chain.id,
                    'distance_A': round(float(dist), 3),
                    'atom_type': atype,
                })

        coord_atoms.sort(key=lambda x: x['distance_A'])

        # Detect porphyrin: >= 4 N atoms from heme residue
        heme_N = [a for a in coord_atoms
                  if a['element'] == 'N' and a['atom_type'] == 'heme/cofactor']
        is_heme = len(heme_N) >= 4

        results.append({
            'metal_name': mi['metal_name'],
            'metal_atom': metal_atom.get_name().strip(),
            'metal_residue': mi['residue'].get_resname().strip(),
            'metal_resnum': mi['residue'].id[1],
            'metal_chain': mi['chain'].id,
            'coordination_number': len(coord_atoms),
            'geometry': classify_geometry(len(coord_atoms)),
            'is_heme_iron': is_heme,
            'coordinating_atoms': coord_atoms,
        })

    os.makedirs(results_dir, exist_ok=True)

    # ── CSV: one row per coordinating bond ──
    csv_rows = []
    for r in results:
        for ca in r['coordinating_atoms']:
            csv_rows.append({
                'metal_name': r['metal_name'],
                'metal_atom': r['metal_atom'],
                'metal_residue': r['metal_residue'],
                'metal_resnum': r['metal_resnum'],
                'metal_chain': r['metal_chain'],
                'coordination_number': r['coordination_number'],
                'geometry': r['geometry'],
                'is_heme_iron': r['is_heme_iron'],
                'coord_atom': ca['atom_name'],
                'coord_element': ca['element'],
                'coord_residue': ca['residue'],
                'coord_resnum': ca['resnum'],
                'coord_chain': ca['chain'],
                'bond_distance_A': ca['distance_A'],
                'ligand_type': ca['atom_type'],
            })

    csv_path = os.path.join(results_dir, 'metal_coordination.csv')
    if csv_rows:
        with open(csv_path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=list(csv_rows[0].keys()))
            writer.writeheader()
            writer.writerows(csv_rows)

    # ── Text summary ──
    txt_path = os.path.join(results_dir, 'metal_coordination_summary.txt')
    with open(txt_path, 'w', encoding='utf-8') as f:
        f.write("Metal Ion Coordination Analysis\n")
        f.write("=" * 50 + "\n\n")
        if not results:
            f.write("No metal ions detected in this structure.\n")
        else:
            f.write(f"Metal ions found: {len(results)}\n\n")
            for r in results:
                f.write("-" * 45 + "\n")
                f.write(f"Metal:   {r['metal_name']} atom [{r['metal_atom']}]\n")
                f.write(f"Site:    {r['metal_residue']} {r['metal_resnum']}  Chain {r['metal_chain']}\n")
                f.write(f"Coord #: {r['coordination_number']}  Geometry: {r['geometry']}\n")
                if r['is_heme_iron']:
                    f.write(f"NOTE:    Fe in porphyrin/heme — 4 equatorial N + axial ligands\n")
                f.write(f"\nCoordinating atoms:\n")
                for ca in r['coordinating_atoms']:
                    f.write(f"  [{ca['atom_type']:12s}]  {ca['residue']:4s} {ca['resnum']:4d}"
                            f"  {ca['atom_name']:4s} ({ca['element']})  "
                            f"{ca['distance_A']:.3f} A  Chain {ca['chain']}\n")
                f.write("\n")

    print(f"Metal coordination analysis: {len(results)} metal ion(s) found")
    return results


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print("Usage: python analyze_metal_coordination.py <bound.pdb> <ligand_code> <results_dir>")
        sys.exit(1)
    analyze_metal_coordination(sys.argv[1], sys.argv[2], sys.argv[3])
