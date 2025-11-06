import pandas as pd
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, average_precision_score, confusion_matrix
import numpy as np
import os

def evaluate_predictions(predictions_file, ground_truth_file):
    """
    Evaluates model predictions against ground truth labels using a robust coordinate-based merge.

    Args:
        predictions_file (str): Path to the TSV file with model predictions (e.g., evo2_predictions.tsv).
                                Expected columns: 'chrom', 'pos_0based', 'ref', 'alt', 'label', 'delta_logp'.
        ground_truth_file (str): Path to the TSV file used as input for predictions (e.g., splice_variants_for_evo2.tsv).
                                 Expected columns: 'chrom', 'pos_0based', 'ref', 'alt', 'label'.
    """
    try:
        # Load the prediction and ground truth files
        df_preds = pd.read_csv(predictions_file, sep='\t')
        df_gt = pd.read_csv(ground_truth_file, sep='\t')
        print("Files loaded successfully.")
        print(f"Prediction file has {len(df_preds)} rows.")
        print(f"Ground truth input file has {len(df_gt)} rows.")

        # --- Data Preprocessing ---
        # For clarity, rename the label columns before merging to avoid conflicts.
        df_preds.rename(columns={'label': 'pred_label'}, inplace=True)
        df_gt.rename(columns={'label': 'gt_label'}, inplace=True)

        # The key for merging will be the genomic coordinates, which is robust.
        merge_cols = ['chrom', 'pos_0based', 'ref', 'alt']
        
        # --- Merging Data ---
        # Merge based on the robust genomic coordinates.
        df_merged = pd.merge(df_preds, df_gt, on=merge_cols, how='inner')
        
        print(f"Successfully merged {len(df_merged)} variants for comparison.")

        if len(df_merged) == 0:
            print("\nError: No variants could be matched between the prediction and ground truth input files.")
            print("Please check the coordinate columns ('chrom', 'pos_0based', 'ref', 'alt') in both files.")
            return

        # --- Performance Evaluation ---
        # y_true is the ground truth from the input file.
        y_true = df_merged['gt_label'].apply(lambda x: 1 if x.lower() == 'splice-altering' else 0)
        
        # y_score is the raw score from the model. For AUPRC, higher score should mean higher probability of being the positive class (1).
        # In the data, a large negative delta_logp indicates "splice-altering".
        # We therefore use the negative of delta_logp as our score.
        y_score = -df_merged['delta_logp']
        
        # y_pred is the binary prediction derived from the score, based on a logical threshold.
        # Based on the data, any delta_logp < 0 indicates a splice-altering prediction by the model.
        y_pred = df_merged['delta_logp'].apply(lambda x: 1 if x < 0 else 0)
        
        accuracy = accuracy_score(y_true, y_pred)
        # Use zero_division=0 to avoid errors if a class is never predicted.
        precision = precision_score(y_true, y_pred, zero_division=0)
        recall = recall_score(y_true, y_pred, zero_division=0)
        f1 = f1_score(y_true, y_pred, zero_division=0)
        auprc = average_precision_score(y_true, y_score)
        cm = confusion_matrix(y_true, y_pred)

        print("\n--- Model Performance Evaluation ---")
        print(f"Accuracy:  {accuracy:.4f}")
        print(f"Precision: {precision:.4f}")
        print(f"Recall:    {recall:.4f}")
        print(f"F1 Score:  {f1:.4f}")
        print(f"AUPRC:     {auprc:.4f} (calculated using 'delta_logp' as score)")
        
        print("\nConfusion Matrix:")
        print("             Predicted Normal  Predicted Splice-Altering")
        # Ensure there are two classes before trying to access cm[0] and cm[1]
        if cm.shape == (2, 2):
            print(f"Actual Normal      {cm[0][0]:<15} {cm[0][1]:<15}")
            print(f"Actual Splice-Alt  {cm[1][0]:<15} {cm[1][1]:<15}")
        else:
            print("Confusion matrix is not 2x2, printing raw matrix:")
            print(cm)

    except FileNotFoundError as e:
        print(f"Error: {e}. Please ensure the file paths are correct.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    # --- Configuration ---
    # Using a robust method to define paths relative to the script's location,
    # which makes the script independent of the current working directory.
    try:
        SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
    except NameError:
        # Fallback for interactive environments like Jupyter notebooks
        SCRIPT_DIR = os.getcwd()

    PREDICTIONS_FILE_PATH = os.path.join(SCRIPT_DIR, 'splicevardb_data/evo2_predictions_evo2_splicevardb_dataset_dedup.tsv')
    GROUND_TRUTH_FILE_PATH = os.path.join(SCRIPT_DIR, 'splicevardb_data/evo2_splicevardb_dataset_dedup.tsv')
    
    print("Starting evaluation with final benchmark files...")
    evaluate_predictions(PREDICTIONS_FILE_PATH, GROUND_TRUTH_FILE_PATH)
    print("\nEvaluation finished.") 