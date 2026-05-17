# Drug–Protein Binding Visualization Tool

**Live Demo:** [drug-binding-visualizer.onrender.com](https://drug-binding-visualizer.onrender.com)

A web-based computational tool for analyzing drug–protein interactions. Upload any two PDB files — an AlphaFold predicted structure and an experimental drug-bound structure — and the tool compares them, identifies the binding pocket, and renders an interactive 3D model in the browser.

---

## Overview

This project analyzes how a drug molecule binds to a protein using structural biology data. The tool was built to make protein-ligand analysis accessible through a browser interface, without requiring any local software installation.

The workflow combines AlphaFold predicted protein structures, experimental Protein Data Bank (PDB) structures, Python structural analysis using Biopython, and interactive 3D visualization using NGL Viewer.

The primary test case uses EGFR (Epidermal Growth Factor Receptor) and gefitinib, a targeted cancer therapy. EGFR plays a key role in cell signaling and growth, and mutations in the receptor are associated with several cancers. Understanding how drugs bind to EGFR helps researchers design more effective targeted therapies.

---

## Features

- Upload any AlphaFold PDB and drug-bound PDB for comparison
- Auto-fetch bound structures directly from RCSB PDB by entering a PDB ID
- Detects all protein residues within 5 Angstroms of the ligand
- Interactive 3D viewer with adjustable slider to isolate protein, drug, or both simultaneously
- Generates structure comparison graphs (atom count, residue count, chain count)
- Displays binding pocket residues in the browser

---

## Technologies

- Python, Flask
- Biopython, Pandas, Matplotlib
- NGL Viewer (browser-based 3D molecular visualization)
- RCSB PDB REST API
- Deployed on Render

---

## How to Use

**Option 1 — Live demo**

Visit the live URL above. Upload an AlphaFold PDB file and a drug-bound PDB file, enter the ligand residue code, and click Run Analysis.

**Option 2 — Run locally**

Clone the repository and install dependencies:

```
git clone https://github.com/AbhiKotra/drug-binding-visualization-tool-website.git
cd drug-binding-visualization-tool-website
pip install -r requirements.txt
python app.py
```

Then open `http://localhost:5000` in your browser.

---

## Project Workflow

1. Upload or fetch two PDB files (unbound AlphaFold structure + drug-bound experimental structure)
2. The tool parses both files and counts atoms, residues, and chains
3. Residues within 5 Angstroms of the ligand are identified as the binding pocket
4. A structure comparison chart is generated
5. Results are displayed alongside an interactive 3D viewer

---

## Project Structure

```
drug-binding-visualization-tool-website
│
├── app.py                    Flask web server
├── requirements.txt
│
├── scripts
│   ├── analyze_structure.py      Counts atoms, residues, chains
│   ├── find_binding_pocket.py    Identifies residues within 5A of ligand
│   └── plot_summary.py           Generates comparison bar charts
│
├── templates
│   └── index.html                Frontend web interface
│
└── data
    ├── alphafold                 AlphaFold PDB input files
    ├── pdb                       Drug-bound PDB input files
    ├── uploads                   Temporary uploaded files
    ├── results                   CSV and TXT output files
    └── images                    Generated graph PNGs
```

---

## Output Files

Each analysis generates:

- `structure_summary.csv` — atom, residue, and chain counts for both structures
- `binding_pocket_residues.txt` — list of residues within 5 Angstroms of the ligand
- `structure_comparison.png` — bar chart comparing the two structures

---

## Future Improvements

- Automatic ligand detection without requiring a residue code
- Support for comparing multiple inhibitors simultaneously
- Per-residue RMSD calculation between bound and unbound states
- Machine learning analysis of binding site characteristics

---

## Author

Abhi Kotra
Biomedical Engineering, UNC Charlotte
github.com/AbhiKotra
