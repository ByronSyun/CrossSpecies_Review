#!/usr/bin/env python3
"""
Evaluate Evo2 predictions on ratGTEx binary classification task.
Computes AUROC, AUPRC, and other classification metrics.
"""

import pandas as pd
import numpy as np
import argparse
import json
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import (
    roc_auc_score, average_precision_score, roc_curve, 
    precision_recall_curve, confusion_matrix
)
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')


def load_predictions(predictions_file):
    """Load Evo2 predictions"""
    logging.info(f"Loading predictions from: {predictions_file}")
    df = pd.read_csv(predictions_file, sep='\t')
    logging.info(f"Loaded {len(df)} predictions")
    logging.info(f"Columns: {list(df.columns)}")
    
    # Check for required columns
    required = ['variant_id', 'delta_logp', 'label']
    for col in required:
        if col not in df.columns:
            raise ValueError(f"Missing required column: {col}")
    
    return df


def evaluate_performance(df, output_dir):
    """
    Evaluate Evo2 performance on rat data.
    
    Args:
        df: DataFrame with predictions and labels
        output_dir: Output directory for results
    """
    logging.info("\n" + "=" * 60)
    logging.info("Evaluating Evo2 performance on ratGTEx")
    logging.info("=" * 60)
    
    # Extract scores and labels
    scores = df['delta_logp'].values
    labels = df['label'].values
    
    # Handle NaN values
    valid_mask = ~np.isnan(scores)
    if not valid_mask.all():
        n_nan = (~valid_mask).sum()
        logging.warning(f"Found {n_nan} NaN scores. Removing them.")
        scores = scores[valid_mask]
        labels = labels[valid_mask]
    
    if len(scores) == 0:
        logging.error("No valid scores left after NaN removal.")
        return None
    
    # Calculate metrics
    try:
        auroc = roc_auc_score(labels, scores)
        auprc = average_precision_score(labels, scores)
    except Exception as e:
        logging.error(f"Error calculating metrics: {e}")
        return None
    
    # Correlation
    from scipy.stats import pearsonr
    corr, p_val = pearsonr(scores, labels)
    
    # Score distribution
    pos_scores = scores[labels == 1]
    neg_scores = scores[labels == 0]
    
    logging.info(f"\n{'=' * 60}")
    logging.info("PERFORMANCE METRICS")
    logging.info(f"{'=' * 60}")
    logging.info(f"Total variants:     {len(scores):,}")
    logging.info(f"Positive samples:   {np.sum(labels):,} ({np.mean(labels)*100:.1f}%)")
    logging.info(f"Negative samples:   {(len(labels) - np.sum(labels)):,} ({(1-np.mean(labels))*100:.1f}%)")
    logging.info(f"")
    logging.info(f"Score statistics:")
    logging.info(f"  Range: [{np.min(scores):.6f}, {np.max(scores):.6f}]")
    logging.info(f"  Mean:  {np.mean(scores):.6f}")
    logging.info(f"  Std:   {np.std(scores):.6f}")
    logging.info(f"")
    logging.info(f"{'=' * 60}")
    logging.info(f"AUROC: {auroc:.4f}")
    logging.info(f"AUPRC: {auprc:.4f}")
    logging.info(f"{'=' * 60}")
    logging.info(f"")
    logging.info(f"Pearson correlation: r={corr:.4f}, p={p_val:.4e}")
    logging.info(f"")
    logging.info(f"Score distribution by label:")
    logging.info(f"  Normal (0):         mean={np.mean(neg_scores):.6f}, std={np.std(neg_scores):.6f}")
    logging.info(f"  Splice-altering (1): mean={np.mean(pos_scores):.6f}, std={np.std(pos_scores):.6f}")
    
    # Plot ROC and PR curves
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # ROC Curve
    fpr, tpr, roc_thresholds = roc_curve(labels, scores)
    ax1.plot(fpr, tpr, label=f'AUROC = {auroc:.4f}', linewidth=2)
    ax1.plot([0, 1], [0, 1], 'k--', label='Random', alpha=0.5)
    ax1.set_xlabel('False Positive Rate', fontsize=12)
    ax1.set_ylabel('True Positive Rate', fontsize=12)
    ax1.set_title('ROC Curve (Evo2 on Rat)', fontsize=14)
    ax1.legend(loc='lower right', fontsize=11)
    ax1.grid(alpha=0.3)
    
    # PR Curve
    precision, recall, pr_thresholds = precision_recall_curve(labels, scores)
    ax2.plot(recall, precision, label=f'AUPRC = {auprc:.4f}', linewidth=2)
    ax2.axhline(y=np.mean(labels), color='k', linestyle='--', 
                label=f'Baseline = {np.mean(labels):.2f}', alpha=0.5)
    ax2.set_xlabel('Recall', fontsize=12)
    ax2.set_ylabel('Precision', fontsize=12)
    ax2.set_title('Precision-Recall Curve', fontsize=14)
    ax2.legend(loc='best', fontsize=11)
    ax2.grid(alpha=0.3)
    
    plt.tight_layout()
    
    # Save figure
    fig_path = output_dir / 'evo2_rat_curves.png'
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    logging.info(f"\nSaved performance curves to: {fig_path}")
    plt.close()
    
    # Calculate confusion matrix at optimal threshold (Youden's J)
    j_scores = tpr - fpr
    optimal_idx = np.argmax(j_scores)
    optimal_threshold = roc_thresholds[optimal_idx]
    
    y_pred_optimal = (scores >= optimal_threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(labels, y_pred_optimal).ravel()
    
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
    ppv = tp / (tp + fp) if (tp + fp) > 0 else 0
    npv = tn / (tn + fn) if (tn + fn) > 0 else 0
    
    logging.info(f"\nOptimal threshold (Youden's J): {optimal_threshold:.6f}")
    logging.info(f"Confusion Matrix:")
    logging.info(f"  TN: {tn:,}  FP: {fp:,}")
    logging.info(f"  FN: {fn:,}  TP: {tp:,}")
    logging.info(f"")
    logging.info(f"Sensitivity: {sensitivity:.4f}")
    logging.info(f"Specificity: {specificity:.4f}")
    logging.info(f"PPV:         {ppv:.4f}")
    logging.info(f"NPV:         {npv:.4f}")
    
    # Plot score distributions
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Histogram
    ax1.hist(neg_scores, bins=50, alpha=0.7, label='Normal', density=True, color='blue')
    ax1.hist(pos_scores, bins=50, alpha=0.7, label='Splice-altering', density=True, color='red')
    ax1.set_xlabel('Delta LogP Score', fontsize=12)
    ax1.set_ylabel('Density', fontsize=12)
    ax1.set_title('Score Distribution by Class', fontsize=14)
    ax1.legend(fontsize=11)
    ax1.grid(alpha=0.3)
    
    # Box plot
    data_to_plot = [neg_scores, pos_scores]
    ax2.boxplot(data_to_plot, labels=['Normal', 'Splice-altering'])
    ax2.set_ylabel('Delta LogP Score', fontsize=12)
    ax2.set_title('Score Distribution Box Plot', fontsize=14)
    ax2.grid(alpha=0.3)
    
    plt.tight_layout()
    dist_path = output_dir / 'evo2_rat_score_distribution.png'
    plt.savefig(dist_path, dpi=300, bbox_inches='tight')
    logging.info(f"Saved score distribution plot to: {dist_path}")
    plt.close()
    
    # Save results to JSON
    results = {
        'n_variants': int(len(scores)),
        'n_positive': int(np.sum(labels)),
        'n_negative': int(len(labels) - np.sum(labels)),
        'positive_rate': float(np.mean(labels)),
        'auroc': float(auroc),
        'auprc': float(auprc),
        'pearson_r': float(corr),
        'pearson_p': float(p_val),
        'optimal_threshold': float(optimal_threshold),
        'true_negatives': int(tn),
        'false_positives': int(fp),
        'false_negatives': int(fn),
        'true_positives': int(tp),
        'sensitivity': float(sensitivity),
        'specificity': float(specificity),
        'ppv': float(ppv),
        'npv': float(npv),
        'score_mean': float(np.mean(scores)),
        'score_std': float(np.std(scores)),
        'score_min': float(np.min(scores)),
        'score_max': float(np.max(scores))
    }
    
    results_path = output_dir / 'evo2_rat_results.json'
    with open(results_path, 'w') as f:
        json.dump(results, f, indent=2)
    logging.info(f"Saved detailed results to: {results_path}")
    
    return results


def main():
    parser = argparse.ArgumentParser(
        description="Evaluate Evo2 predictions on ratGTEx data"
    )
    parser.add_argument(
        '--predictions', 
        type=str, 
        default='results/evo2_rat_predictions.tsv',
        help='Path to Evo2 prediction TSV file'
    )
    parser.add_argument(
        '--output_dir', 
        type=str, 
        default='results',
        help='Directory to save evaluation results'
    )
    
    args = parser.parse_args()
    
    # Load predictions (labels are already in the predictions file)
    df = load_predictions(args.predictions)
    
    # Evaluate
    results = evaluate_performance(df, args.output_dir)
    
    if results:
        logging.info(f"\n{'=' * 60}")
        logging.info("Evaluation complete!")
        logging.info(f"Results saved to: {args.output_dir}")
        logging.info(f"{'=' * 60}")


if __name__ == '__main__':
    main()

