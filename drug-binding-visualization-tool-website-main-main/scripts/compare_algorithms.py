import os
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, silhouette_score


BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(BASE_DIR, "data", "results")
IMAGES_DIR = os.path.join(BASE_DIR, "data", "images")

INPUT_FILE = os.path.join(RESULTS_DIR, "ml_dataset_clean.csv")
OUTPUT_CSV = os.path.join(RESULTS_DIR, "algorithm_comparison.csv")
OUTPUT_IMAGE = os.path.join(IMAGES_DIR, "algorithm_comparison.png")


def main():
    if not os.path.exists(INPUT_FILE):
        print("Clean dataset not found.")
        return

    df = pd.read_csv(INPUT_FILE)

    df["significant_structural_change"] = (
        df["global_rmsd"] >= 20
    ).astype(int)

    feature_columns = [
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
    y = df["significant_structural_change"]

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    results = []

    # Random Forest
    rf = RandomForestClassifier(n_estimators=100, random_state=42)
    rf.fit(X, y)
    rf_preds = rf.predict(X)
    rf_acc = accuracy_score(y, rf_preds)

    results.append({
        "algorithm": "Random Forest",
        "task_type": "classification",
        "score_type": "training_accuracy",
        "score": rf_acc,
        "notes": "Supervised model predicting high vs low structural impact"
    })

    # Logistic Regression
    if y.nunique() >= 2:
        lr = LogisticRegression(max_iter=1000)
        lr.fit(X_scaled, y)
        lr_preds = lr.predict(X_scaled)
        lr_acc = accuracy_score(y, lr_preds)

        results.append({
            "algorithm": "Logistic Regression",
            "task_type": "classification",
            "score_type": "training_accuracy",
            "score": lr_acc,
            "notes": "Linear supervised classifier for structural impact"
        })

    # K-Means
    kmeans = KMeans(n_clusters=3, random_state=42, n_init=10)
    clusters = kmeans.fit_predict(X_scaled)

    if len(set(clusters)) > 1 and len(df) > 3:
        sil_score = silhouette_score(X_scaled, clusters)
    else:
        sil_score = 0

    results.append({
        "algorithm": "K-Means",
        "task_type": "clustering",
        "score_type": "silhouette_score",
        "score": sil_score,
        "notes": "Unsupervised clustering of structural similarity"
    })

    results_df = pd.DataFrame(results)
    results_df.to_csv(OUTPUT_CSV, index=False)

    plt.figure(figsize=(9, 6))
    plt.bar(results_df["algorithm"], results_df["score"])

    plt.title("ML Algorithm Comparison")
    plt.ylabel("Score")
    plt.xlabel("Algorithm")
    plt.ylim(0, 1.05)

    for i, value in enumerate(results_df["score"]):
        plt.text(
            i,
            value + 0.02,
            f"{value:.2f}",
            ha="center"
        )

    os.makedirs(IMAGES_DIR, exist_ok=True)

    plt.tight_layout()
    plt.savefig(OUTPUT_IMAGE, dpi=300, bbox_inches="tight")
    plt.close()

    print("Algorithm comparison complete.")
    print(f"Saved CSV to: {OUTPUT_CSV}")
    print(f"Saved image to: {OUTPUT_IMAGE}")
    print(results_df)


if __name__ == "__main__":
    main()