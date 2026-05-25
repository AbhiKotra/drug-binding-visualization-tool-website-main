import os
import pandas as pd

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(BASE_DIR, "data", "results")

OUTPUT_FILE = os.path.join(RESULTS_DIR, "ml_dataset.csv")


def main():

    dataset_rows = []

    for session_id in os.listdir(RESULTS_DIR):

        session_path = os.path.join(RESULTS_DIR, session_id)

        if not os.path.isdir(session_path):
            continue

        metrics_path = os.path.join(session_path, "advanced_metrics.csv")

        if not os.path.exists(metrics_path):
            continue

        try:
            metrics_df = pd.read_csv(metrics_path)

            if metrics_df.empty:
                continue

            metrics = metrics_df.iloc[0].to_dict()

            dataset_rows.append({
                "session_id": session_id,

                "global_rmsd": metrics.get("global_rmsd"),
                "local_rmsd_5A": metrics.get("local_rmsd_5A"),
                "local_rmsd_8A": metrics.get("local_rmsd_8A"),
                "local_rmsd_10A": metrics.get("local_rmsd_10A"),

                "matched_residues_5A": metrics.get("matched_residues_5A"),
                "matched_residues_8A": metrics.get("matched_residues_8A"),
                "matched_residues_10A": metrics.get("matched_residues_10A"),

                "radius_of_gyration_alphafold":
                    metrics.get("radius_of_gyration_alphafold"),

                "radius_of_gyration_bound":
                    metrics.get("radius_of_gyration_bound")
            })

        except Exception as e:
            print(f"Skipping {session_id}: {e}")

    ml_df = pd.DataFrame(dataset_rows)

    ml_df.to_csv(OUTPUT_FILE, index=False)

    print("\nML dataset created successfully.")
    print(f"Saved to: {OUTPUT_FILE}")
    print(ml_df.head())


if __name__ == "__main__":
    main()