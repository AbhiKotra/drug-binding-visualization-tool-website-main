#!/usr/bin/env python3
"""
ramachandran_plot.py
Calculates phi/psi backbone dihedral angles and generates Ramachandran plots
for both AlphaFold and drug-bound structures (useful for validating ML-built models).

Usage: python ramachandran_plot.py <alphafold.pdb> <bound.pdb> <images_dir> <results_dir>
"""

from Bio.PDB import PDBParser, PPBuilder
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os
import sys
import csv
import math


def classify_region(phi, psi):
    """Classify phi/psi into Ramachandran regions."""
    # Core alpha-helix
    if -120 <= phi <= -30 and -80 <= psi <= 50:
        return 'alpha_helix'
    # Core beta-sheet (including right-twisted)
    if -170 <= phi <= -50 and 60 <= psi <= 180:
        return 'beta_sheet'
    # Left-handed helix (common for Gly)
    if 30 <= phi <= 100 and -20 <= psi <= 90:
        return 'left_helix'
    return 'outlier'


def get_phi_psi(structure):
    """Return list of dicts with chain, resname, resnum, phi, psi."""
    ppb = PPBuilder()
    data = []
    for pp in ppb.build_peptides(structure):
        phi_psi_list = pp.get_phi_psi_list()
        for i, (phi, psi) in enumerate(phi_psi_list):
            if phi is None or psi is None:
                continue
            res = pp[i]
            data.append({
                'chain':   res.get_parent().id,
                'residue': res.get_resname().strip(),
                'resnum':  res.id[1],
                'phi':     round(math.degrees(phi), 2),
                'psi':     round(math.degrees(psi), 2),
            })
    return data


def draw_ramachandran(phi_psi_data, label, images_dir, results_dir, dot_color):
    if not phi_psi_data:
        print(f"  No phi/psi angles for {label} — skipping plot")
        return None

    phi_vals   = [d['phi'] for d in phi_psi_data]
    psi_vals   = [d['psi'] for d in phi_psi_data]
    regions    = [classify_region(p, s) for p, s in zip(phi_vals, psi_vals)]

    region_color = {
        'alpha_helix': '#FF6B35',
        'beta_sheet':  '#00e5ff',
        'left_helix':  '#22d3a4',
        'outlier':     '#71717a',
    }
    colors = [region_color[r] for r in regions]

    counts = {k: regions.count(k) for k in region_color}
    total  = len(regions)
    pct_favored = round(100 * (counts['alpha_helix'] + counts['beta_sheet']) / total, 1) if total else 0.0

    # ── Plot ──
    fig, ax = plt.subplots(figsize=(8, 7))
    fig.patch.set_facecolor('#0d0d10')
    ax.set_facecolor('#111114')

    # Shaded allowed-region backgrounds
    ax.add_patch(plt.Rectangle((-120, -80),  90, 130,
                                color='#FF6B35', alpha=0.10, zorder=0, label='_nolegend_'))
    ax.add_patch(plt.Rectangle((-170,  60), 130, 120,
                                color='#00e5ff', alpha=0.10, zorder=0, label='_nolegend_'))
    ax.add_patch(plt.Rectangle(( 30,  -20),  70, 110,
                                color='#22d3a4', alpha=0.08, zorder=0, label='_nolegend_'))

    # Reference axes
    ax.axhline(0, color='#444450', linewidth=0.6, linestyle='--', zorder=1)
    ax.axvline(0, color='#444450', linewidth=0.6, linestyle='--', zorder=1)

    ax.scatter(phi_vals, psi_vals, c=colors, s=16, alpha=0.80, zorder=5, linewidths=0)

    ax.set_xlabel('φ (Phi) / degrees', color='#e4e4e7', fontsize=12)
    ax.set_ylabel('ψ (Psi) / degrees', color='#e4e4e7', fontsize=12)
    ax.set_title(f'Ramachandran Plot — {label}', color='#e4e4e7', fontsize=14, fontweight='bold')
    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.set_xticks(range(-180, 181, 60))
    ax.set_yticks(range(-180, 181, 60))
    ax.tick_params(colors='#71717a', labelsize=9)
    for spine in ax.spines.values():
        spine.set_edgecolor('#333340')

    legend_handles = [
        mpatches.Patch(color='#FF6B35', label=f'α-helix ({counts["alpha_helix"]})'),
        mpatches.Patch(color='#00e5ff', label=f'β-sheet ({counts["beta_sheet"]})'),
        mpatches.Patch(color='#22d3a4', label=f'Left helix ({counts["left_helix"]})'),
        mpatches.Patch(color='#71717a', label=f'Outlier ({counts["outlier"]})'),
    ]
    ax.legend(handles=legend_handles, loc='upper right',
              facecolor='#1a1a1e', edgecolor='#333340',
              labelcolor='#e4e4e7', fontsize=9)

    ax.text(0.02, 0.02,
            f'Total residues: {total}\nFavoured: {pct_favored}%',
            transform=ax.transAxes, color='#71717a', fontsize=8,
            verticalalignment='bottom', fontfamily='monospace')

    plt.tight_layout()

    safe_label = label.lower().replace(' ', '_').replace('/', '_')
    plot_path  = os.path.join(images_dir, f'ramachandran_{safe_label}.png')
    plt.savefig(plot_path, dpi=150, bbox_inches='tight', facecolor='#0d0d10')
    plt.close()

    # ── CSV ──
    csv_path = os.path.join(results_dir, f'phi_psi_{safe_label}.csv')
    with open(csv_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=['chain', 'residue', 'resnum', 'phi', 'psi', 'region'])
        writer.writeheader()
        for d, reg in zip(phi_psi_data, regions):
            d['region'] = reg
            writer.writerow(d)

    print(f"  Ramachandran ({label}): {total} residues, {pct_favored}% favoured — {plot_path}")
    return {
        'plot_path': plot_path, 'total': total,
        'pct_favored': pct_favored, 'counts': counts,
    }


if __name__ == '__main__':
    if len(sys.argv) < 5:
        print("Usage: python ramachandran_plot.py <alphafold.pdb> <bound.pdb> <images_dir> <results_dir>")
        sys.exit(1)

    af_path    = sys.argv[1]
    bound_path = sys.argv[2]
    images_dir = sys.argv[3]
    results_dir = sys.argv[4]

    os.makedirs(images_dir,  exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)

    parser = PDBParser(QUIET=True)

    print("Generating Ramachandran plots...")
    af_struct    = parser.get_structure('af',    af_path)
    bound_struct = parser.get_structure('bound', bound_path)

    af_data    = get_phi_psi(af_struct)
    bound_data = get_phi_psi(bound_struct)

    draw_ramachandran(af_data,    'AlphaFold',  images_dir, results_dir, '#2196F3')
    draw_ramachandran(bound_data, 'Drug-Bound', images_dir, results_dir, '#00e5cc')
    print("Ramachandran plots complete.")
