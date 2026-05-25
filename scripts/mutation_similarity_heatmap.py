import os
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler


BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

RESULTS_DIR = os.path.join(BASE_DIR, "data", "results")
IMAGES_DIR = os.path.join(BASE_DIR, "data", "images")

INPUT_FILE = os.path.join(RESULTS_DIR, "ml_dataset_clean.csv")

OUTPUT_IMAGE = os.path.join(
    IMAGES_DIR,
    "mutation_similarity_heatmap.png"
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

    X = df[feature_columns]

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    similarity_matrix = pd.DataFrame(
        X_scaled,
        index=df["session_id"]
    ).T.corr()

    plt.figure(figsize=(10, 8))

    plt.imshow(similarity_matrix, aspect='auto')

    plt.colorbar(label="Structural Similarity")

    short_labels = [sid[:6] for sid in df["session_id"]]

    plt.xticks(
        range(len(short_labels)),
        short_labels,
        rotation=90,
        fontsize=8
    )

    plt.yticks(
        range(len(short_labels)),
        short_labels,
        fontsize=8
    )

    plt.title("Mutation / Structure Similarity Heatmap")

    plt.tight_layout()

    os.makedirs(IMAGES_DIR, exist_ok=True)

    plt.savefig(
        OUTPUT_IMAGE,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    print("Mutation similarity heatmap created.")
    print(f"Saved to: {OUTPUT_IMAGE}")


if __name__ == "__main__":
    main()