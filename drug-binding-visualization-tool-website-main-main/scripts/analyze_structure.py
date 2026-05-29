from Bio.PDB import PDBParser, Superimposer
import os
import pandas as pd
import sys
import math
import numpy as np


def count_atoms_residues(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    atom_count = 0
    residue_count = 0
    chain_ids = set()

    for model in structure:
        for chain in model:
            chain_ids.add(chain.id)
            for residue in chain:
                residue_count += 1
                for atom in residue:
                    atom_count += 1

    return {
        "file": os.path.basename(pdb_file),
        "chains": len(chain_ids),
        "residues": residue_count,
        "atoms": atom_count
    }


def get_ca_atoms(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    ca_atoms = []

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == " " and "CA" in residue:
                    ca_atoms.append(residue["CA"])

    return ca_atoms


def get_ca_atom_map(pdb_file):
    """
    Creates a dictionary of CA atoms by:
    chain ID + residue number + insertion code

    Example key:
    ('A', 858, ' ')
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    ca_map = {}

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == " " and "CA" in residue:
                    chain_id = chain.id
                    residue_number = residue.id[1]
                    insertion_code = residue.id[2]
                    key = (chain_id, residue_number, insertion_code)
                    ca_map[key] = residue["CA"]

    return ca_map


def get_ligand_atoms(bound_file, ligand_code):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("bound", bound_file)

    ligand_atoms = []

    for model in structure:
        for chain in model:
            for residue in chain:
                hetflag = residue.id[0]

                if hetflag != " ":
                    if residue.get_resname().strip().upper() == ligand_code.upper():
                        for atom in residue:
                            ligand_atoms.append(atom)

    return ligand_atoms


def get_residues_near_ligand(bound_file, ligand_code, radius):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("bound", bound_file)

    ligand_atoms = get_ligand_atoms(bound_file, ligand_code)

    nearby_keys = set()

    if len(ligand_atoms) == 0:
        return nearby_keys

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != " ":
                    continue

                found = False

                for atom in residue:
                    for ligand_atom in ligand_atoms:

                        distance = atom - ligand_atom

                        if distance <= radius:

                            chain_id = chain.id
                            residue_number = residue.id[1]
                            insertion_code = residue.id[2]

                            key = (chain_id, residue_number, insertion_code)

                            nearby_keys.add(key)

                            found = True
                            break

                    if found:
                        break

    return nearby_keys


def calculate_global_rmsd(file1, file2):
    atoms1 = get_ca_atoms(file1)
    atoms2 = get_ca_atoms(file2)

    min_len = min(len(atoms1), len(atoms2))

    if min_len < 3:
        return None

    atoms1 = atoms1[:min_len]
    atoms2 = atoms2[:min_len]

    super_imposer = Superimposer()
    super_imposer.set_atoms(atoms1, atoms2)

    return round(super_imposer.rms, 3)


def calculate_true_local_rmsd(alphafold_file, bound_file, ligand_code, radius):
    """
    True local RMSD:
    1. Find residues in the bound PDB within radius of ligand.
    2. Match those residues to AlphaFold by chain + residue number.
    3. Superimpose only those local CA atoms.
    4. Return RMSD for that local binding-site region.
    """
    alphafold_map = get_ca_atom_map(alphafold_file)
    bound_map = get_ca_atom_map(bound_file)

    nearby_keys = get_residues_near_ligand(bound_file, ligand_code, radius)

    alphafold_atoms = []
    bound_atoms = []

    for key in nearby_keys:
        if key in alphafold_map and key in bound_map:
            alphafold_atoms.append(alphafold_map[key])
            bound_atoms.append(bound_map[key])

    matched_residue_count = len(alphafold_atoms)

    if matched_residue_count < 3:
        return None, matched_residue_count

    super_imposer = Superimposer()
    super_imposer.set_atoms(alphafold_atoms, bound_atoms)

    return round(super_imposer.rms, 3), matched_residue_count


def calculate_radius_of_gyration(pdb_file):
    atoms = get_ca_atoms(pdb_file)

    if len(atoms) == 0:
        return None

    coords = np.array([atom.coord for atom in atoms])
    center = coords.mean(axis=0)

    squared_distances = np.sum((coords - center) ** 2, axis=1)
    rg = math.sqrt(squared_distances.mean())

    return round(rg, 3)


if len(sys.argv) < 5:
    print("Usage: python analyze_structure.py <alphafold.pdb> <bound.pdb> <results_dir> <ligand_code>")
    sys.exit(1)

alphafold_path = sys.argv[1]
bound_path = sys.argv[2]
results_dir = sys.argv[3]
ligand_code = sys.argv[4].strip().upper()

os.makedirs(results_dir, exist_ok=True)

files = [alphafold_path, bound_path]
structure_rows = []

for f in files:
    if os.path.exists(f):
        structure_rows.append(count_atoms_residues(f))
    else:
        print(f"Missing file: {f}")

structure_df = pd.DataFrame(structure_rows)

structure_csv = os.path.join(results_dir, "structure_summary.csv")
structure_df.to_csv(structure_csv, index=False)

global_rmsd = calculate_global_rmsd(alphafold_path, bound_path)

local_5A, matched_5A = calculate_true_local_rmsd(
    alphafold_path, bound_path, ligand_code, 5.0
)

local_8A, matched_8A = calculate_true_local_rmsd(
    alphafold_path, bound_path, ligand_code, 8.0
)

local_10A, matched_10A = calculate_true_local_rmsd(
    alphafold_path, bound_path, ligand_code, 10.0
)

alphafold_rg = calculate_radius_of_gyration(alphafold_path)
bound_rg = calculate_radius_of_gyration(bound_path)

metrics = [{
    "alphafold_file": os.path.basename(alphafold_path),
    "bound_file": os.path.basename(bound_path),
    "ligand_code": ligand_code,
    "global_rmsd": global_rmsd,
    "local_rmsd_5A": local_5A,
    "local_rmsd_8A": local_8A,
    "local_rmsd_10A": local_10A,
    "matched_residues_5A": matched_5A,
    "matched_residues_8A": matched_8A,
    "matched_residues_10A": matched_10A,
    "radius_of_gyration_alphafold": alphafold_rg,
    "radius_of_gyration_bound": bound_rg
}]

metrics_df = pd.DataFrame(metrics)

metrics_csv = os.path.join(results_dir, "advanced_metrics.csv")
metrics_df.to_csv(metrics_csv, index=False)

print("Basic structure summary:")
print(structure_df)

print("\nAdvanced metrics:")
print(metrics_df)

print("\nSaved results to:")
print(structure_csv)
print(metrics_csv)