import os
import pandas as pd

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

RESULTS_DIR = os.path.join(BASE_DIR, "data", "results")

INPUT_FILE = os.path.join(RESULTS_DIR, "ml_dataset.csv")
OUTPUT_FILE = os.path.join(RESULTS_DIR, "ml_dataset_clean.csv")


def main():

    if not os.path.exists(INPUT_FILE):
        print("ml_dataset.csv not found.")
        return

    df = pd.read_csv(INPUT_FILE)

    print("\nOriginal dataset shape:")
    print(df.shape)

    # Remove duplicate sessions
    df = df.drop_duplicates(subset=["session_id"])

    # Convert columns to numeric
    numeric_columns = [
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

    for col in numeric_columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    # Fill missing local RMSD values
    df["local_rmsd_5A"] = df["local_rmsd_5A"].fillna(0)
    df["local_rmsd_8A"] = df["local_rmsd_8A"].fillna(0)
    df["local_rmsd_10A"] = df["local_rmsd_10A"].fillna(0)

    # Fill missing matched residues
    df["matched_residues_5A"] = df["matched_residues_5A"].fillna(0)
    df["matched_residues_8A"] = df["matched_residues_8A"].fillna(0)
    df["matched_residues_10A"] = df["matched_residues_10A"].fillna(0)

    # Remove rows missing global RMSD
    df = df.dropna(subset=["global_rmsd"])

    # Reset index
    df = df.reset_index(drop=True)

    # Save cleaned dataset
    df.to_csv(OUTPUT_FILE, index=False)

    print("\nCleaned dataset shape:")
    print(df.shape)

    print("\nSaved cleaned dataset:")
    print(OUTPUT_FILE)

    print("\nPreview:")
    print(df.head())


if __name__ == "__main__":
    main()