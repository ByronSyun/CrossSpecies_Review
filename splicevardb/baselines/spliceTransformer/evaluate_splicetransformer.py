import pandas as pd
from sklearn.metrics import (
    roc_auc_score,
    precision_recall_curve,
    auc,
    confusion_matrix,
    accuracy_score,
    precision_score,
    recall_score,
    f1_score,
)


def evaluate_splicetransformer_predictions(predictions_file, ground_truth_file, threshold=0.2):
    """Evaluate SpliceTransformer predictions against SpliceVarDB labels."""
    print(f"Loading predictions: {predictions_file}")
    try:
        pred_df = pd.read_csv(predictions_file, index_col=0)
    except Exception as e:
        print(f"Failed to read predictions: {e}")
        return

    print(f"Loading ground truth: {ground_truth_file}")
    try:
        gt_df = pd.read_csv(ground_truth_file, sep='\t', low_memory=False)
    except Exception as e:
        print(f"Failed to read ground truth: {e}")
        return

    # Normalise IDs and build merge keys
    pred_df['#CHROM'] = pred_df['#CHROM'].apply(lambda x: x if 'chr' in str(x) else 'chr' + str(x))
    pred_df['variant_id'] = (
        pred_df['#CHROM'].astype(str) + '-' +
        pred_df['POS'].astype(str) + '-' +
        pred_df['REF'].astype(str) + '-' +
        pred_df['ALT'].astype(str)
    )

    gt_df['label'] = gt_df['label'].apply(lambda x: 1 if x == 'splice-altering' else 0)
    gt_df['POS'] = gt_df['pos_0based'] + 1
    gt_df['variant_id'] = (
        gt_df['chrom'].astype(str) + '-' +
        gt_df['POS'].astype(str) + '-' +
        gt_df['ref'].astype(str) + '-' +
        gt_df['alt'].astype(str)
    )
    gt_df = gt_df[['variant_id', 'label']]

    merged_df = pd.merge(pred_df, gt_df, on='variant_id', how='inner')
    if merged_df.empty:
        print("Empty merge; check variant_id formats.")
        return
    print(f"Merged variants: {len(merged_df)}")

    y_true = merged_df['label']
    y_scores = merged_df['score']
    y_pred = (y_scores >= threshold).astype(int)

    auroc = roc_auc_score(y_true, y_scores)
    precision, recall, _ = precision_recall_curve(y_true, y_scores)
    auprc = auc(recall, precision)

    acc = accuracy_score(y_true, y_pred)
    prec = precision_score(y_true, y_pred)
    rec = recall_score(y_true, y_pred)
    f1 = f1_score(y_true, y_pred)

    print(f"AUROC: {auroc:.4f} | AUPRC: {auprc:.4f}")
    print(f"Accuracy: {acc:.4f} | Precision: {prec:.4f} | Recall: {rec:.4f} | F1: {f1:.4f}")
    print("Confusion matrix:\n", confusion_matrix(y_true, y_pred))


if __name__ == "__main__":
    PREDICTIONS_FILE = '/mnt/userdata4/splicing/SpliceTransformer/data/SpliceT_predictions.tsv'
    GROUND_TRUTH_FILE = '/mnt/userdata4/splicing/Evo2Splicevardb/splicevardb_data/evo2_splicevardb_dataset.tsv'
    evaluate_splicetransformer_predictions(PREDICTIONS_FILE, GROUND_TRUTH_FILE)