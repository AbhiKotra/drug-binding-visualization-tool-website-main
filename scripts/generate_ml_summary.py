import os
import pandas as pd


BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(BASE_DIR, "data", "results")

DATASET_FILE = os.path.join(RESULTS_DIR, "ml_dataset_clean.csv")
PREDICTIONS_FILE = os.path.join(RESULTS_DIR, "mutation_impact_predictions.csv")
FEATURE_FILE = os.path.join(RESULTS_DIR, "feature_importance.csv")
PCA_FILE = os.path.join(RESULTS_DIR, "pca_cluster_results.csv")

OUTPUT_FILE = os.path.join(RESULTS_DIR, "ml_summary_report.txt")


def main():
    report_lines = []

    report_lines.append("BioStruct-AI Machine Learning Summary Report")
    report_lines.append("=" * 55)
    report_lines.append("")

    if os.path.exists(DATASET_FILE):
        df = pd.read_csv(DATASET_FILE)

        report_lines.append("Dataset Overview")
        report_lines.append("-" * 20)
        report_lines.append(f"Total analysis sessions: {len(df)}")
        report_lines.append(f"Average global RMSD: {df['global_rmsd'].mean():.3f}")
        report_lines.append(f"Maximum global RMSD: {df['global_rmsd'].max():.3f}")
        report_lines.append(f"Minimum global RMSD: {df['global_rmsd'].min():.3f}")
        report_lines.append("")

    if os.path.exists(PREDICTIONS_FILE):
        pred_df = pd.read_csv(PREDICTIONS_FILE)

        report_lines.append("Structural Impact Classification")
        report_lines.append("-" * 35)

        if "impact_label" in pred_df.columns:
            counts = pred_df["impact_label"].value_counts()

            for label, count in counts.items():
                report_lines.append(f"{label}: {count} sessions")

        report_lines.append("")

    if os.path.exists(FEATURE_FILE):
        feature_df = pd.read_csv(FEATURE_FILE)

        report_lines.append("Feature Importance Interpretation")
        report_lines.append("-" * 35)

        top_features = feature_df.head(5)

        for _, row in top_features.iterrows():
            report_lines.append(
                f"{row['feature']} contributed an importance score of {row['importance']:.3f}."
            )

        top_feature = feature_df.iloc[0]["feature"]

        report_lines.append("")
        report_lines.append(
            f"The strongest predictor in the current dataset was {top_feature}."
        )
        report_lines.append("")

    if os.path.exists(PCA_FILE):
        pca_df = pd.read_csv(PCA_FILE)

        report_lines.append("PCA and Clustering Summary")
        report_lines.append("-" * 30)

        if "cluster" in pca_df.columns:
            cluster_counts = pca_df["cluster"].value_counts().sort_index()

            for cluster_id, count in cluster_counts.items():
                report_lines.append(f"Cluster {cluster_id}: {count} sessions")

        report_lines.append("")

    report_lines.append("Scientific Interpretation")
    report_lines.append("-" * 30)
    report_lines.append(
        "The current machine learning pipeline suggests that local structural deviation "
        "metrics and radius of gyration are useful indicators of protein structural impact."
    )
    report_lines.append(
        "PCA and clustering analyses help identify structurally similar protein-ligand "
        "cases and highlight outlier conformations for further biological investigation."
    )
    report_lines.append(
        "Because the dataset is still small, the current prediction results should be "
        "interpreted as exploratory rather than definitive."
    )

    with open(OUTPUT_FILE, "w", encoding="utf-8") as f:
        f.write("\n".join(report_lines))

    print("ML summary report created.")
    print(f"Saved to: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()