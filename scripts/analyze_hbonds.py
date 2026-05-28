#!/usr/bin/env python3
"""
analyze_hbonds.py
Detects protein-ligand hydrogen bond interactions.

For structures WITH explicit H atoms:
  - Reports D-H covalent bond distance (~0.96 A for N-H, ~0.82 A for O-H in X-ray)
  - Reports H...A distance
  - Computes D-H...A angle: >90 deg = weak, >120 deg = moderate, >160 deg = strong

For structures WITHOUT H atoms (typical X-ray):
  - Uses D...A heavy-atom distance criterion (<= 3.5 A)
  - Classifies by distance: <2.5 strong, 2.5-3.2 moderate, 3.2-3.5 weak

Usage: python analyze_hbonds.py <bound.pdb> <ligand_code> <results_dir>
"""

from Bio.PDB import PDBParser
import numpy as np
import os
import sys
import csv
import math

ACCEPTOR_ELEMENTS = {'O', 'N', 'S', 'F'}
DONOR_ELEMENTS = {'O', 'N'}

# Maximum D...A distance to consider
DA_MAX_DIST = 3.5
# Maximum H...A distance when H is present
HA_MAX_DIST = 2.5
# Maximum D-H covalent bond distance
DH_MAX_DIST = 1.2


def get_element(atom):
    elem = atom.element
    if elem:
        return elem.strip().upper()
    name = atom.get_name().strip()
    return name[0].upper() if name else ''


def find_bonded_h(donor_atom, all_atoms):
    """Find H atoms covalently bonded to a donor atom (within ~1.2 A)."""
    bonded = []
    for a in all_atoms:
        if get_element(a) == 'H':
            dist = donor_atom - a
            if dist < DH_MAX_DIST:
                bonded.append(a)
    return bonded


def calc_angle_deg(coord_d, coord_h, coord_a):
    """Calculate D-H...A angle in degrees (angle at H)."""
    vec_hd = np.array(coord_d) - np.array(coord_h)
    vec_ha = np.array(coord_a) - np.array(coord_h)
    norm_hd = np.linalg.norm(vec_hd)
    norm_ha = np.linalg.norm(vec_ha)
    if norm_hd == 0 or norm_ha == 0:
        return None
    cos_a = np.dot(vec_hd, vec_ha) / (norm_hd * norm_ha)
    cos_a = float(np.clip(cos_a, -1.0, 1.0))
    return math.degrees(math.acos(cos_a))


def classify_by_distance(d_a_dist):
    if d_a_dist < 2.5:
        return 'Strong'
    elif d_a_dist < 3.2:
        return 'Moderate'
    elif d_a_dist <= 3.5:
        return 'Weak'
    return None


def classify_by_angle(angle_deg):
    if angle_deg > 160:
        return 'Strong'
    elif angle_deg > 120:
        return 'Moderate'
    elif angle_deg > 90:
        return 'Weak'
    return None


def analyze_hbonds(bound_file, ligand_code, results_dir):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('complex', bound_file)

    all_atoms = list(structure.get_atoms())

    ligand_atoms = []
    protein_atoms = []

    for model in structure:
        for chain in model:
            for residue in chain:
                hetflag = residue.id[0]
                if hetflag != ' ':
                    if residue.get_resname().strip().upper() == ligand_code.upper():
                        for atom in residue:
                            ligand_atoms.append((atom, residue, chain))
                elif hetflag == ' ':
                    for atom in residue:
                        protein_atoms.append((atom, residue, chain))

    has_h = any(get_element(a) == 'H' for a in all_atoms)

    hbonds = []
    seen = set()

    for p_atom, p_res, p_chain in protein_atoms:
        p_elem = get_element(p_atom)
        if p_elem not in DONOR_ELEMENTS and p_elem not in ACCEPTOR_ELEMENTS:
            continue

        for l_atom, l_res, l_chain in ligand_atoms:
            l_elem = get_element(l_atom)
            if l_elem not in ACCEPTOR_ELEMENTS and l_elem not in DONOR_ELEMENTS:
                continue

            da_dist = p_atom - l_atom
            if da_dist > DA_MAX_DIST:
                continue

            # Determine donor / acceptor direction
            if p_elem in DONOR_ELEMENTS and l_elem in ACCEPTOR_ELEMENTS:
                direction = 'Protein -> Ligand'
                donor_atom, donor_res_name, donor_resnum, donor_chain_id = (
                    p_atom, p_res.get_resname().strip(), p_res.id[1], p_chain.id)
                acc_atom, acc_res_name, acc_resnum, acc_chain_id = (
                    l_atom, l_res.get_resname().strip(), l_res.id[1], l_chain.id)
                donor_label = f"{donor_res_name} {donor_resnum} (Chain {donor_chain_id})"
                acc_label = f"{acc_res_name} (Ligand)"
            elif l_elem in DONOR_ELEMENTS and p_elem in ACCEPTOR_ELEMENTS:
                direction = 'Ligand -> Protein'
                donor_atom, donor_res_name, donor_resnum, donor_chain_id = (
                    l_atom, l_res.get_resname().strip(), l_res.id[1], l_chain.id)
                acc_atom, acc_res_name, acc_resnum, acc_chain_id = (
                    p_atom, p_res.get_resname().strip(), p_res.id[1], p_chain.id)
                donor_label = f"{donor_res_name} (Ligand)"
                acc_label = f"{acc_res_name} {acc_resnum} (Chain {acc_chain_id})"
            else:
                continue

            dedup_key = (donor_label, donor_atom.get_name().strip(),
                         acc_label, acc_atom.get_name().strip())
            if dedup_key in seen:
                continue
            seen.add(dedup_key)

            entry = {
                'direction': direction,
                'donor_residue': donor_label,
                'donor_atom': donor_atom.get_name().strip(),
                'donor_element': get_element(donor_atom),
                'acceptor_residue': acc_label,
                'acceptor_atom': acc_atom.get_name().strip(),
                'acceptor_element': get_element(acc_atom),
                'da_distance_A': round(da_dist, 3),
                'has_explicit_H': False,
                'dh_distance_A': 'N/A',
                'ha_distance_A': 'N/A',
                'dha_angle_deg': 'N/A',
                'strength': classify_by_distance(da_dist),
                'geometry_note': 'D...A distance only (no H in structure)',
            }

            # Try to refine with explicit H
            if has_h:
                bonded_h = find_bonded_h(donor_atom, all_atoms)
                for h_atom in bonded_h:
                    ha_dist = h_atom - acc_atom
                    if ha_dist > HA_MAX_DIST:
                        continue
                    dh_dist = donor_atom - h_atom
                    angle = calc_angle_deg(donor_atom.coord, h_atom.coord, acc_atom.coord)
                    if angle is None:
                        continue
                    strength_angle = classify_by_angle(angle)
                    if strength_angle is None:
                        continue
                    entry.update({
                        'has_explicit_H': True,
                        'dh_distance_A': round(dh_dist, 3),
                        'ha_distance_A': round(ha_dist, 3),
                        'dha_angle_deg': round(angle, 1),
                        'strength': strength_angle,
                        'geometry_note': 'X-H...Y angle geometry (explicit H)',
                    })
                    break

            hbonds.append(entry)

    os.makedirs(results_dir, exist_ok=True)

    fields = [
        'direction', 'donor_residue', 'donor_atom', 'donor_element',
        'acceptor_residue', 'acceptor_atom', 'acceptor_element',
        'da_distance_A', 'has_explicit_H', 'dh_distance_A', 'ha_distance_A',
        'dha_angle_deg', 'strength', 'geometry_note',
    ]

    csv_path = os.path.join(results_dir, 'hbond_interactions.csv')
    with open(csv_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(hbonds)

    strong   = [h for h in hbonds if h['strength'] == 'Strong']
    moderate = [h for h in hbonds if h['strength'] == 'Moderate']
    weak     = [h for h in hbonds if h['strength'] == 'Weak']

    txt_path = os.path.join(results_dir, 'hbond_summary.txt')
    with open(txt_path, 'w', encoding='utf-8') as f:
        f.write(f"Hydrogen Bond Analysis: {ligand_code} <-> Protein\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Explicit H atoms present in structure: {has_h}\n")
        f.write(f"Total H-bond interactions: {len(hbonds)}\n")
        f.write(f"  Strong   (>160 deg / D...A < 2.5 A):   {len(strong)}\n")
        f.write(f"  Moderate (>120 deg / D...A 2.5-3.2 A): {len(moderate)}\n")
        f.write(f"  Weak     (>90  deg / D...A 3.2-3.5 A): {len(weak)}\n\n")

        for hb in hbonds:
            f.write(f"[{hb['strength'].upper()}] {hb['direction']}\n")
            f.write(f"  Donor:    {hb['donor_residue']}  Atom: {hb['donor_atom']} ({hb['donor_element']})\n")
            f.write(f"  Acceptor: {hb['acceptor_residue']}  Atom: {hb['acceptor_atom']} ({hb['acceptor_element']})\n")
            f.write(f"  D...A distance: {hb['da_distance_A']} A\n")
            if hb['has_explicit_H']:
                f.write(f"  D-H bond: {hb['dh_distance_A']} A  |  H...A: {hb['ha_distance_A']} A  |  D-H...A angle: {hb['dha_angle_deg']} deg\n")
            f.write(f"  {hb['geometry_note']}\n\n")

    print(f"H-bond analysis complete: {len(hbonds)} interactions (strong={len(strong)}, moderate={len(moderate)}, weak={len(weak)})")
    return hbonds


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print("Usage: python analyze_hbonds.py <bound.pdb> <ligand_code> <results_dir>")
        sys.exit(1)
    analyze_hbonds(sys.argv[1], sys.argv[2], sys.argv[3])
