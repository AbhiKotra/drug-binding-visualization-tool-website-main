import os
import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, accuracy_score


BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(BASE_DIR, "data", "results")

INPUT_FILE = os.path.join(RESULTS_DIR, "ml_dataset_clean.csv")
OUTPUT_FILE = os.path.join(RESULTS_DIR, "mutation_impact_predictions.csv")


def main():
    if not os.path.exists(INPUT_FILE):
        print("Clean dataset not found. Run clean_ml_dataset.py first.")
        return

    df = pd.read_csv(INPUT_FILE)

    # Create simple target label
    df["significant_structural_change"] = (
        df["global_rmsd"] >= 20
    ).astype(int)

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
    y = df["significant_structural_change"]

    # If dataset is too small or only one class exists, skip split
    class_counts = y.value_counts()

    if len(df) < 8 or y.nunique() < 2 or class_counts.min() < 2:
        print("Dataset too small or only one class found.")
        print("Creating rule-based predictions instead.")

        df["predicted_impact"] = df["significant_structural_change"]
        df["impact_label"] = df["predicted_impact"].map({
            0: "Low structural impact",
            1: "High structural impact"
        })

        df.to_csv(OUTPUT_FILE, index=False)

        print(f"Saved rule-based predictions to: {OUTPUT_FILE}")
        print(df[["session_id", "global_rmsd", "impact_label"]])
        return

    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=0.3,
        random_state=42,
        stratify=y
    )

    model = RandomForestClassifier(
        n_estimators=100,
        random_state=42
    )

    model.fit(X_train, y_train)

    y_pred = model.predict(X_test)

    print("Model accuracy:")
    print(accuracy_score(y_test, y_pred))

    print("\nClassification report:")
    print(classification_report(y_test, y_pred))

    df["predicted_impact"] = model.predict(X)
    df["impact_label"] = df["predicted_impact"].map({
        0: "Low structural impact",
        1: "High structural impact"
    })

    df.to_csv(OUTPUT_FILE, index=False)

    print(f"\nSaved predictions to: {OUTPUT_FILE}")
    print(df[["session_id", "global_rmsd", "impact_label"]])


if __name__ == "__main__":
    main()