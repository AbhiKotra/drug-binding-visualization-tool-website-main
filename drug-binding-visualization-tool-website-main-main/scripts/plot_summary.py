import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys


if len(sys.argv) < 3:
    print("Usage: python plot_summary.py <results_dir> <images_dir>")
    sys.exit(1)

results_dir = sys.argv[1]
images_dir = sys.argv[2]

structure_csv = os.path.join(results_dir, "structure_summary.csv")
metrics_csv = os.path.join(results_dir, "advanced_metrics.csv")

os.makedirs(images_dir, exist_ok=True)


# Plot 1: original atom/residue/chain comparison
structure_df = pd.read_csv(structure_csv)

output_image = os.path.join(images_dir, "structure_comparison.png")

fig, axes = plt.subplots(1, 3, figsize=(14, 5))
fig.suptitle("Structure Comparison Summary", fontsize=14, fontweight='bold')

colors = ['#2196F3', '#FF5722']

axes[0].bar(structure_df["file"], structure_df["atoms"], color=colors)
axes[0].set_title("Atom Count")
axes[0].set_ylabel("Number of Atoms")
axes[0].tick_params(axis='x', rotation=15)

axes[1].bar(structure_df["file"], structure_df["residues"], color=colors)
axes[1].set_title("Residue Count")
axes[1].set_ylabel("Number of Residues")
axes[1].tick_params(axis='x', rotation=15)

axes[2].bar(structure_df["file"], structure_df["chains"], color=colors)
axes[2].set_title("Chain Count")
axes[2].set_ylabel("Number of Chains")
axes[2].tick_params(axis='x', rotation=15)

plt.tight_layout()
plt.savefig(output_image, dpi=150)
plt.close()

print("Saved plot to:", output_image)


# Plot 2: RMSD metrics
if os.path.exists(metrics_csv):
    metrics_df = pd.read_csv(metrics_csv)

    rmsd_values = [
        metrics_df.loc[0, "global_rmsd"],
        metrics_df.loc[0, "local_rmsd_5A"],
        metrics_df.loc[0, "local_rmsd_8A"],
        metrics_df.loc[0, "local_rmsd_10A"]
    ]

    rmsd_labels = [
        "Global RMSD",
        "Local 5A",
        "Local 8A",
        "Local 10A"
    ]

    plt.figure(figsize=(8, 5))
    plt.bar(rmsd_labels, rmsd_values)
    plt.title("RMSD Structural Deviation Metrics")
    plt.ylabel("RMSD")
    plt.xticks(rotation=20)
    plt.tight_layout()

    rmsd_image = os.path.join(images_dir, "rmsd_metrics.png")
    plt.savefig(rmsd_image, dpi=150)
    plt.close()

    print("Saved plot to:", rmsd_image)


    # Plot 3: Radius of gyration
    rg_values = [
        metrics_df.loc[0, "radius_of_gyration_alphafold"],
        metrics_df.loc[0, "radius_of_gyration_bound"]
    ]

    rg_labels = [
        "AlphaFold",
        "Bound/PDB"
    ]

    plt.figure(figsize=(7, 5))
    plt.bar(rg_labels, rg_values)
    plt.title("Radius of Gyration Comparison")
    plt.ylabel("Radius of Gyration")
    plt.tight_layout()

    rg_image = os.path.join(images_dir, "radius_of_gyration.png")
    plt.savefig(rg_image, dpi=150)
    plt.close()

    print("Saved plot to:", rg_image)