from Bio.PDB import PDBParser
import os
import sys

# ─────────────────────────────────────────────────────────────────
# This script now accepts file paths as arguments from app.py.
# Usage: python find_binding_pocket.py <bound.pdb> <ligand_code> <results_dir>
# ─────────────────────────────────────────────────────────────────

if len(sys.argv) < 4:
    print("Usage: python find_binding_pocket.py <bound.pdb> <ligand_code> <results_dir>")
    sys.exit(1)

pdb_file    = sys.argv[1]   # e.g. /tmp/uploads/bound.pdb
ligand_code = sys.argv[2]   # e.g. "IRE" or "gefitinib residue name"
results_dir = sys.argv[3]   # e.g. /tmp/results/

output_file = os.path.join(results_dir, "binding_pocket_residues.txt")

parser    = PDBParser(QUIET=True)
structure = parser.get_structure("protein_ligand", pdb_file)

ligand_atoms    = []
protein_residues = []

for model in structure:
    for chain in model:
        for residue in chain:
            hetflag = residue.id[0]

            if hetflag != " ":
                if residue.get_resname() == ligand_code:
                    for atom in residue:
                        ligand_atoms.append(atom)
            else:
                protein_residues.append(residue)

nearby_residues = set()

for residue in protein_residues:
    for atom in residue:
        for lig_atom in ligand_atoms:
            distance = atom - lig_atom
            if distance <= 5.0:
                resname  = residue.get_resname()
                resid    = residue.id[1]
                chain_id = residue.get_parent().id
                nearby_residues.add((chain_id, resname, resid))
                break
        else:
            continue
        break

nearby_residues = sorted(nearby_residues, key=lambda x: (x[0], x[2]))

os.makedirs(results_dir, exist_ok=True)
with open(output_file, "w") as f:
    f.write(f"Binding pocket residues within 5 Å of ligand {ligand_code}:\n\n")
    for chain_id, resname, resid in nearby_residues:
        line = f"Chain {chain_id}: {resname} {resid}"
        print(line)
        f.write(line + "\n")

print("\nSaved binding pocket residues to:", output_file)
