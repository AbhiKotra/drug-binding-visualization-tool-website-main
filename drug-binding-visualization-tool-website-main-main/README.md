# Drug–Protein Binding Visualization Tool

**Live Demo:** [drug-binding-visualizer.onrender.com](https://drug-binding-visualizer.onrender.com)

A web-based computational tool for analyzing drug–protein interactions. Upload any two PDB files — an AlphaFold predicted structure and an experimental drug-bound structure — and the tool compares them, identifies the binding pocket, calculates interaction geometry, and renders an interactive 3D model in the browser.

---

## Overview

This project analyzes how a drug molecule binds to a protein using structural biology data. The tool was built to make protein-ligand analysis accessible through a browser interface, without requiring any local software installation.

The workflow combines AlphaFold predicted protein structures, experimental Protein Data Bank (PDB) structures, Python structural analysis using Biopython, machine learning pipelines, and interactive 3D visualization using NGL Viewer.

The primary test case uses EGFR (Epidermal Growth Factor Receptor) and gefitinib, a targeted cancer therapy. EGFR plays a key role in cell signaling and growth, and mutations in the receptor are associated with several cancers. Understanding how drugs bind to EGFR helps researchers design more effective targeted therapies.

---

## Features

### Structure Analysis
- Upload any AlphaFold PDB and drug-bound PDB for comparison
- Auto-fetch bound structures directly from RCSB PDB by entering a PDB ID
- Global RMSD and local binding-site RMSD (5 Å, 8 Å, 10 Å cutoffs)
- Radius of gyration calculation for both structures
- Detects all protein residues within 5 Å of the ligand (binding pocket)
- Structure comparison bar charts (atoms, residues, chains, RMSD, radius of gyration)

### Hydrogen Bond Analysis
- Detects all protein ↔ ligand hydrogen bond interactions
- For structures with explicit H atoms: reports covalent D–H bond distance (~0.9 Å), H···A distance, and D–H···A angle
- For structures without H atoms: uses D···A heavy-atom distance criterion (≤ 3.5 Å)
- Classifies bonds by geometry: Weak (> 90°), Moderate (> 120°), Strong (> 160°)
- H-bond donor/acceptor residues highlighted magenta in the 3D viewer

### Occupancy & Thermal Ellipsoid (B-Factor) Analysis
- Parses occupancy values directly from PDB ATOM/HETATM records
- Flags residues with occupancy < 1.0 — indicates dual conformations of protein side chains or ligands
- Calculates mean and per-residue B-factors; flags high B-factor residues (mean + 2σ threshold)
- High B-factor correlates with thermal motion, structural flexibility, and weak binding
- Reports ligand occupancy and B-factor separately to assess binding stability

### Metal Ion Coordination
- Detects common metal ions: Fe, Zn, Mg, Ca, Cu, Mn, Co, Ni, Na, K, Pt, Hg, Cd
- Reports all coordinating N/O/S atoms with bond distances (Å)
- Classifies coordination geometry: linear, tetrahedral, square planar, octahedral, etc.
- Special handling for Fe in heme/porphyrin rings: identifies 4 equatorial N atoms from the porphyrin ring plus axial ligands (e.g., histidine, O₂)

### Ramachandran Plot
- Calculates φ (phi) and ψ (psi) backbone dihedral angles for every residue
- Generates separate plots for AlphaFold and drug-bound structures
- Classifies residues into α-helix, β-sheet, left-handed helix, and outlier regions
- Reports percentage of residues in favoured regions — used to validate ML-predicted structures

### Sequence Comparison
- Extracts amino acid sequences from both PDB files
- Performs global pairwise alignment (PairwiseAligner, match +2 / mismatch −1 / gap −5)
- Reports sequence identity % and conservative-substitution similarity %
- Generates residue-by-residue correspondence table: conserved / conservative sub / mutation
- Flags structures with ≥ 90% identity as same protein or orthologue

### NCBI BLAST Orthologue Search
- Queries NCBI BLAST (blastp vs PDB database) for similar structures from other species
- Configurable identity threshold (default ≥ 80%)
- Returns PDB IDs, chain IDs, organism names, identity %, E-values, alignment previews
- Runs as a separate request so the main analysis stays fast (BLAST takes 30–90 s)

### Machine Learning Pipeline
- Aggregates structural metrics across all sessions into a dataset
- Random Forest classification of structural impact (global RMSD ≥ 20 Å threshold)
- Feature importance analysis — identifies which structural metrics are most predictive
- PCA clustering (k=3) of all analysed structures
- Mutation similarity heatmap (session-to-session structural correlation)
- Algorithm comparison: Random Forest vs Logistic Regression vs K-Means
- Auto-generated plain-English ML summary report

### 3D Viewer
- Interactive NGL Viewer with left-click rotate, scroll zoom, right-click pan
- AlphaFold structure: blue cartoon
- Drug-bound protein: teal cartoon
- Ligand / drug: orange ball+stick
- Binding pocket residues: yellow ball+stick
- H-bond residues: magenta ball+stick
- 3-state opacity slider: protein only ↔ both ↔ drug only
- Toggle between cartoon, surface, and ball+stick representations

---

## Technologies

- Python, Flask, Gunicorn
- Biopython (PDB parsing, PPBuilder, PairwiseAligner, NCBI BLAST)
- Pandas, NumPy, Matplotlib, scikit-learn
- NGL Viewer v2.0 (browser-based 3D molecular visualization)
- RCSB PDB REST API
- Deployed on Render

---

## How to Use

**Option 1 — Live demo**

Visit the live URL above. Upload an AlphaFold PDB file and a drug-bound PDB file, enter the ligand residue code, and click Run Analysis.

**Option 2 — Run locally**

```
git clone https://github.com/AbhiKotra/drug-binding-visualization-tool-website-main.git
cd drug-binding-visualization-tool-website-main
pip install -r requirements.txt
python app.py
```

Then open `http://localhost:5000` in your browser.

---

## Project Workflow

1. Upload or auto-fetch two PDB files (AlphaFold unbound + experimental drug-bound)
2. Structural metrics computed: atoms, residues, chains, global/local RMSD, radius of gyration
3. Binding pocket identified: all residues within 5 Å of the ligand
4. Hydrogen bond interactions detected with distance and angle geometry
5. Occupancy and B-factor values parsed; dual conformations and flexible regions flagged
6. Metal ion coordination bonds detected and classified
7. Ramachandran φ/ψ plots generated for both structures
8. Sequence alignment run between the two structures
9. Results displayed in the browser with interactive 3D viewer
10. Optional: BLAST search finds orthologous structures from other species

---

## Project Structure

```
drug-binding-visualization-tool-website-main
│
├── app.py                         Flask web server + all API routes
├── requirements.txt
├── .gitignore
│
├── scripts
│   ├── analyze_structure.py           RMSD, radius of gyration, atom/residue counts
│   ├── find_binding_pocket.py         Residues within 5 A of ligand
│   ├── plot_summary.py                Structure comparison bar charts
│   ├── analyze_hbonds.py              H-bond detection with distance and angle geometry
│   ├── analyze_occupancy.py           Occupancy and B-factor (thermal ellipsoid) analysis
│   ├── analyze_metal_coordination.py  Metal ion coordination bonds (Fe, Zn, Mg, Ca, Cu...)
│   ├── ramachandran_plot.py           Phi/psi backbone dihedral angle plots
│   ├── sequence_comparison.py         Pairwise sequence alignment and residue correspondence
│   ├── blast_search.py                NCBI BLAST orthologue search vs PDB database
│   ├── build_ml_dataset.py            Aggregates metrics from all sessions
│   ├── clean_ml_dataset.py            Data cleaning and imputation
│   ├── train_mutation_impact_model.py Random Forest structural impact classifier
│   ├── feature_importance_analysis.py Feature importance ranking
│   ├── pca_clustering.py              PCA + K-Means clustering
│   ├── mutation_similarity_heatmap.py Session-to-session correlation heatmap
│   ├── feature_correlation_matrix.py  Pearson correlation matrix
│   ├── compare_algorithms.py          RF vs LR vs K-Means benchmark
│   └── generate_ml_summary.py         Plain-English ML narrative report
│
├── templates
│   └── index.html                     Full frontend (HTML + CSS + JavaScript)
│
└── data
    ├── alphafold                      Reference AlphaFold PDB files
    ├── pdb                            Reference drug-bound PDB files
    ├── uploads                        Temporary uploaded files (gitignored)
    ├── results                        Per-session output files (gitignored)
    └── images                         Generated graph PNGs (gitignored)
```

---

## Output Files (per session)

| File | Description |
|---|---|
| `structure_summary.csv` | Atom, residue, chain counts for both structures |
| `advanced_metrics.csv` | Global/local RMSD, matched residues, radius of gyration |
| `binding_pocket_residues.txt` | Residues within 5 Å of ligand |
| `hbond_interactions.csv` | All H-bond interactions with distance and angle |
| `hbond_summary.txt` | H-bond analysis report |
| `low_occupancy_residues.csv` | Residues with occupancy < 1.0 |
| `high_bfactor_residues.csv` | Residues above B-factor threshold |
| `occupancy_bfactor_summary.txt` | Occupancy and thermal ellipsoid report |
| `metal_coordination.csv` | Metal ion coordination bond table |
| `metal_coordination_summary.txt` | Metal coordination report |
| `phi_psi_alphafold.csv` | Ramachandran angles for AlphaFold structure |
| `phi_psi_drug-bound.csv` | Ramachandran angles for drug-bound structure |
| `sequence_comparison.csv` | Pairwise alignment identity and similarity |
| `residue_correspondence.csv` | Residue-by-residue conservation map |
| `blast_results.csv` | NCBI BLAST hits from other species |
| `structure_comparison.png` | Bar chart: atoms/residues/chains |
| `rmsd_metrics.png` | RMSD comparison chart |
| `radius_of_gyration.png` | Radius of gyration chart |
| `ramachandran_alphafold.png` | Ramachandran plot for AlphaFold |
| `ramachandran_drug-bound.png` | Ramachandran plot for drug-bound |

---

## Author

Abhi Kotra
Biomedical Engineering, UNC Charlotte
github.com/AbhiKotra
