import os
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.ensemble import RandomForestClassifier


BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(BASE_DIR, "data", "results")
IMAGES_DIR = os.path.join(BASE_DIR, "data", "images")

INPUT_FILE = os.path.join(RESULTS_DIR, "ml_dataset_clean.csv")
OUTPUT_CSV = os.path.join(RESULTS_DIR, "feature_importance.csv")
OUTPUT_IMAGE = os.path.join(IMAGES_DIR, "feature_importance.png")


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

    if y.nunique() < 2:
        print("Only one class found. Cannot compute feature importance.")
        return

    model = RandomForestClassifier(
        n_estimators=200,
        random_state=42
    )

    model.fit(X, y)

    importance_df = pd.DataFrame({
        "feature": feature_columns,
        "importance": model.feature_importances_
    }).sort_values(by="importance", ascending=False)

    importance_df.to_csv(OUTPUT_CSV, index=False)

    plt.figure(figsize=(10, 6))
    plt.barh(
        importance_df["feature"],
        importance_df["importance"]
    )
    plt.gca().invert_yaxis()

    plt.title("Feature Importance for Structural Impact Classification")
    plt.xlabel("Importance Score")
    plt.ylabel("Structural Feature")

    os.makedirs(IMAGES_DIR, exist_ok=True)

    plt.tight_layout()
    plt.savefig(OUTPUT_IMAGE, dpi=300, bbox_inches="tight")
    plt.close()

    print("Feature importance analysis complete.")
    print(f"Saved CSV to: {OUTPUT_CSV}")
    print(f"Saved image to: {OUTPUT_IMAGE}")
    print(importance_df)


if __name__ == "__main__":
    main()