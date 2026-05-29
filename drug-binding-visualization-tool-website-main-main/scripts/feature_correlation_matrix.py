import os
import pandas as pd
import matplotlib.pyplot as plt


BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

RESULTS_DIR = os.path.join(BASE_DIR, "data", "results")
IMAGES_DIR = os.path.join(BASE_DIR, "data", "images")

INPUT_FILE = os.path.join(RESULTS_DIR, "ml_dataset_clean.csv")

OUTPUT_IMAGE = os.path.join(
    IMAGES_DIR,
    "feature_correlation_matrix.png"
)


def main():

    if not os.path.exists(INPUT_FILE):
        print("Clean dataset not found.")
        return

    df = pd.read_csv(INPUT_FILE)

    feature_columns = [
        "global_rmsd",
        "local_rmsd_5A",
        "local_rmsd_8A",
        "local_rmsd_10A",
        "matched_residues_5A",
        "matched_residues_8A",
        "matched_residues_10A",
        "radius_of_gyration_alphafold",
        "radius_of_gyration_bound"
    ]

    correlation_matrix = df[feature_columns].corr()

    plt.figure(figsize=(10, 8))

    im = plt.imshow(
        correlation_matrix,
        interpolation='nearest',
        aspect='auto'
    )

    plt.colorbar(im, label="Correlation")

    plt.xticks(
        range(len(feature_columns)),
        feature_columns,
        rotation=90
    )

    plt.yticks(
        range(len(feature_columns)),
        feature_columns
    )

    # Add correlation values inside cells
    for i in range(len(feature_columns)):
        for j in range(len(feature_columns)):

            value = correlation_matrix.iloc[i, j]

            plt.text(
                j,
                i,
                f"{value:.2f}",
                ha="center",
                va="center",
                fontsize=8
            )

    plt.title("Structural Feature Correlation Matrix")

    plt.tight_layout()

    os.makedirs(IMAGES_DIR, exist_ok=True)

    plt.savefig(
        OUTPUT_IMAGE,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    print("Feature correlation matrix created.")
    print(f"Saved to: {OUTPUT_IMAGE}")


if __name__ == "__main__":
    main()