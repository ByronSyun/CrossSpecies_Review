import argparse
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score

def evaluate_performance(predictions_path, ground_truth_path, output_path, threshold):
    """
    Evaluates SpliceAI performance by comparing its predictions against a ground truth file.
    Calculates both threshold-based metrics and threshold-independent metrics (AUROC/AUPRC).

    Args:
        predictions_path (str): Path to the TSV file with SpliceAI scores (e.g., spliceai_scores.tsv).
        ground_truth_path (str): Path to the ground truth TSV file (e.g., splice_variants_for_evo2.tsv).
        output_path (str): Path to save the merged and evaluated results.
        threshold (float): The delta score threshold to classify a variant as 'splice-altering'.
    """
    print(f"Loading predictions from: {predictions_path}")
    try:
        predictions_df = pd.read_csv(predictions_path, sep='\t')
        # Standardize chromosome column if it exists
        if 'CHROM' in predictions_df.columns:
            predictions_df['CHROM'] = predictions_df['CHROM'].astype(str).str.replace('^chr', '', regex=True)
    except FileNotFoundError:
        print(f"Error: Predictions file not found at {predictions_path}")
        return

    print(f"Loading ground truth from: {ground_truth_path}")
    try:
        ground_truth_df = pd.read_csv(ground_truth_path, sep='\t', usecols=['chrom', 'pos_0based', 'ref', 'alt', 'label'])
        ground_truth_df.rename(columns={'label': 'True_Label'}, inplace=True)
        # Standardize chromosome column
        ground_truth_df['chrom'] = ground_truth_df['chrom'].astype(str).str.replace('^chr', '', regex=True)
        # Create a 1-based position for merging with VCF-based data
        ground_truth_df['POS'] = ground_truth_df['pos_0based'] + 1
    except FileNotFoundError:
        print(f"Error: Ground truth file not found at {ground_truth_path}")
        return

    print("Merging prediction scores with true labels using a robust key (CHROM, POS, REF, ALT)...")
    # A robust inner merge ensures we only evaluate variants present in both files.
    # The key now includes CHROM for accuracy.
    merge_cols = ['CHROM', 'POS', 'REF', 'ALT']
    # To prepare for the merge, we rename ground_truth_df columns to match predictions_df
    ground_truth_df.rename(columns={'chrom': 'CHROM', 'ref': 'REF', 'alt': 'ALT'}, inplace=True)
    
    # Ensure dtypes are compatible for merging
    predictions_df['POS'] = predictions_df['POS'].astype(int)
    ground_truth_df['POS'] = ground_truth_df['POS'].astype(int)
    predictions_df['CHROM'] = predictions_df['CHROM'].astype(str)
    ground_truth_df['CHROM'] = ground_truth_df['CHROM'].astype(str)

    eval_df = pd.merge(
        predictions_df, 
        ground_truth_df,
        how='inner',
        on=merge_cols
    )

    if eval_df.empty:
        print("Error: No matching variants found between prediction and ground truth files.")
        print("Please check if chromosome and position formats are consistent (e.g., 'chr1' vs '1').")
        return

    # --- Threshold-based Evaluation (for a specific operating point) ---
    eval_df['Predicted_Label'] = eval_df['MAX_DS'].apply(
        lambda score: 'splice-altering' if score >= threshold else 'normal'
    )
    is_true_altering = eval_df['True_Label'] == 'splice-altering'
    is_pred_altering = eval_df['Predicted_Label'] == 'splice-altering'
    true_positives = len(eval_df[is_true_altering & is_pred_altering])
    false_negatives = len(eval_df[is_true_altering & ~is_pred_altering])
    true_negatives = len(eval_df[~is_true_altering & ~is_pred_altering])
    false_positives = len(eval_df[~is_true_altering & is_pred_altering])
    recall = true_positives / (true_positives + false_negatives) if (true_positives + false_negatives) > 0 else 0
    precision = true_positives / (true_positives + false_positives) if (true_positives + false_positives) > 0 else 0
    accuracy = (true_positives + true_negatives) / len(eval_df) if len(eval_df) > 0 else 0

    # --- Threshold-independent Evaluation (Overall Performance) ---
    y_true = eval_df['True_Label'].apply(lambda x: 1 if x == 'splice-altering' else 0)
    y_score = eval_df['MAX_DS']
    auroc = roc_auc_score(y_true, y_score)
    auprc = average_precision_score(y_true, y_score)


    print("\n--- SpliceAI Performance Evaluation ---")
    print(f"Total variants matched and evaluated: {len(eval_df)}")
    print(f"  - Ground Truth Positives (splice-altering): {y_true.sum()}")
    print(f"  - Ground Truth Negatives (normal): {len(y_true) - y_true.sum()}")
    
    print("\n--- Overall Performance (Threshold-Independent) ---")
    print(f"  - AUROC (Area Under ROC Curve): {auroc:.4f}")
    print(f"  - AUPRC (Area Under PR Curve):  {auprc:.4f}")

    print(f"\n--- Performance at Threshold = {threshold} ---")
    print("Confusion Matrix:")
    print(f"  - True Positives (TP): {true_positives}")
    print(f"  - False Negatives (FN): {false_negatives}")
    print(f"  - True Negatives (TN): {true_negatives}")
    print(f"  - False Positives (FP): {false_positives}")
    print("Metrics at Threshold:")
    print(f"  - Recall (Sensitivity): {recall:.4f}")
    print(f"  - Precision: {precision:.4f}")
    print(f"  - Accuracy: {accuracy:.4f}")
    print("---------------------------------------\n")

    eval_df.to_csv(output_path, sep='\t', index=False)
    print(f"Detailed results saved to: {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Evaluate SpliceAI predictions against a ground truth dataset.")
    parser.add_argument("-P", "--predictions", required=True, help="Path to the parsed SpliceAI scores TSV file.")
    parser.add_argument("-G", "--ground_truth", required=True, help="Path to the ground truth TSV file.")
    parser.add_argument("-O", "--output", required=True, help="Path for the output file with merged results.")
    parser.add_argument("-T", "--threshold", type=float, default=0.2, help="Delta score threshold for specific metrics (default: 0.2).")

    args = parser.parse_args()

    evaluate_performance(args.predictions, args.ground_truth, args.output, args.threshold) 