import os
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans


BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(BASE_DIR, "data", "results")
IMAGES_DIR = os.path.join(BASE_DIR, "data", "images")

INPUT_FILE = os.path.join(RESULTS_DIR, "ml_dataset_clean.csv")
OUTPUT_IMAGE = os.path.join(IMAGES_DIR, "pca_clusters_improved.png")
OUTPUT_CSV = os.path.join(RESULTS_DIR, "pca_cluster_results.csv")


def main():
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

    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_scaled)

    kmeans = KMeans(n_clusters=3, random_state=42, n_init=10)
    clusters = kmeans.fit_predict(X_scaled)

    df["PCA1"] = X_pca[:, 0]
    df["PCA2"] = X_pca[:, 1]
    df["cluster"] = clusters

    df.to_csv(OUTPUT_CSV, index=False)

    plt.figure(figsize=(11, 8))

    for cluster_id in sorted(df["cluster"].unique()):
        cluster_df = df[df["cluster"] == cluster_id]

        plt.scatter(
            cluster_df["PCA1"],
            cluster_df["PCA2"],
            s=140,
            label=f"Cluster {cluster_id}"
        )

    for i, row in df.iterrows():
        label = row["session_id"][:6]

        plt.annotate(
            label,
            (row["PCA1"], row["PCA2"]),
            textcoords="offset points",
            xytext=(6, 6),
            fontsize=8
        )

    plt.title("PCA Clustering of Protein Structural Metrics", fontsize=14)
    plt.xlabel(f"Principal Component 1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)")
    plt.ylabel(f"Principal Component 2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)")

    plt.legend(title="K-Means Cluster")
    plt.grid(True, alpha=0.3)

    os.makedirs(IMAGES_DIR, exist_ok=True)

    plt.savefig(OUTPUT_IMAGE, dpi=300, bbox_inches="tight")
    plt.close()

    print("Improved PCA clustering plot created.")
    print(f"Saved image to: {OUTPUT_IMAGE}")
    print(f"Saved results to: {OUTPUT_CSV}")
    print("Explained variance ratio:", pca.explained_variance_ratio_)


if __name__ == "__main__":
    main()