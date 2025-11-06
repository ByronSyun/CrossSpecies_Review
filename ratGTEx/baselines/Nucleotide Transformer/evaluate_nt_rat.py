#!/usr/bin/env python3
"""
Evaluate Nucleotide Transformer predictions on rat (ratGTEx) data.
"""

import pandas as pd
import numpy as np
import argparse
import json
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import (
    roc_auc_score, average_precision_score, 
    roc_curve, precision_recall_curve,
    confusion_matrix, accuracy_score, f1_score
)

def load_nt_scores(scores_csv):
    """Load NT prediction scores"""
    print(f"\nLoading NT scores from: {scores_csv}")
    df = pd.read_csv(scores_csv)
    print(f"Loaded {len(df):,} predictions")
    print(f"Columns: {list(df.columns)}")
    return df

def load_rat_labels(rat_sequences_tsv):
    """
    Load rat variant labels from the sequences TSV.
    This file contains the ground truth labels.
    """
    print(f"\nLoading rat labels from: {rat_sequences_tsv}")
    # Rat TSV has no header, format: variant_id ref_seq alt_seq label chrom
    df = pd.read_csv(rat_sequences_tsv, sep=r'\s+', header=None,
                     names=['variant_id', 'ref_seq', 'alt_seq', 'label', 'chrom'])
    print(f"Loaded {len(df):,} labeled variants")
    
    # Keep only variant_id and label
    df_labels = df[['variant_id', 'label']].copy()
    
    # Ensure label is binary (0/1)
    unique_labels = df_labels['label'].unique()
    print(f"Label values: {sorted(unique_labels)}")
    if not set(unique_labels).issubset({0, 1}):
        raise ValueError(f"Expected binary labels (0/1), got: {unique_labels}")
    
    print(f"Label distribution: 0={np.sum(df_labels['label']==0):,}, 1={np.sum(df_labels['label']==1):,}")
    
    return df_labels

def merge_predictions_with_labels(df_predictions, df_labels):
    """Merge NT predictions with ground truth labels"""
    print(f"\nMerging predictions with labels...")
    
    # Remove duplicates from labels (rat data may have multi-tissue replicates)
    print(f"  Labels before dedup: {len(df_labels):,}")
    df_labels_dedup = df_labels.drop_duplicates(subset=['variant_id'], keep='first')
    print(f"  Labels after dedup:  {len(df_labels_dedup):,}")
    
    df_merged = df_predictions.merge(df_labels_dedup, on='variant_id', how='inner')
    
    print(f"Successfully merged {len(df_merged):,} variants")
    n_pos = int(df_merged['label'].sum())
    n_neg = int((~df_merged['label'].astype(bool)).sum())
    pos_rate = df_merged['label'].mean() * 100
    
    print(f"  Positive (splice-altering): {n_pos:,} ({pos_rate:.1f}%)")
    print(f"  Negative (normal):          {n_neg:,} ({100-pos_rate:.1f}%)")
    
    # Check for missing variants
    missing_count = len(df_labels) - len(df_merged)
    if missing_count > 0:
        print(f"  Warning: {missing_count} variants in labels not found in predictions")
    
    return df_merged

def evaluate_performance(df_merged, score_column, output_dir):
    """
    Evaluate NT performance using specified score column.
    
    Args:
        df_merged: Merged predictions and labels
        score_column: Name of score column to evaluate ('score_cosine' or 'score_l2')
        output_dir: Output directory for results
    """
    print(f"\n{'='*60}")
    print(f"Evaluating performance using: {score_column}")
    print(f"{'='*60}")
    
    # Extract scores and labels
    scores = df_merged[score_column].values
    labels = df_merged['label'].values
    
    # Handle NaN values
    valid_mask = ~np.isnan(scores)
    if not valid_mask.all():
        n_nan = (~valid_mask).sum()
        print(f"Warning: Found {n_nan} NaN scores. Removing them.")
        scores = scores[valid_mask]
        labels = labels[valid_mask]
    
    if len(scores) == 0:
        print("ERROR: No valid scores for evaluation.")
        return None
    
    # Calculate metrics
    auroc = roc_auc_score(labels, scores)
    auprc = average_precision_score(labels, scores)
    
    print(f"\n{'='*60}")
    print(f"PERFORMANCE METRICS ({score_column})")
    print(f"{'='*60}")
    print(f"Total variants:     {len(scores):,}")
    print(f"Positive samples:   {np.sum(labels):,} ({np.mean(labels)*100:.1f}%)")
    print(f"Negative samples:   {(len(labels) - np.sum(labels)):,} ({(1-np.mean(labels))*100:.1f}%)")
    print(f"\nScore statistics:")
    print(f"  Range: [{np.min(scores):.6e}, {np.max(scores):.6e}]")
    print(f"  Mean:  {np.mean(scores):.6e}")
    print(f"  Std:   {np.std(scores):.6e}")
    print(f"\n{'='*60}")
    print(f"AUROC: {auroc:.4f}")
    print(f"AUPRC: {auprc:.4f}")
    print(f"{'='*60}")
    
    # Calculate correlation with labels
    from scipy.stats import pearsonr
    corr, pval = pearsonr(scores, labels)
    print(f"\nPearson correlation with labels: r={corr:.4f}, p={pval:.4e}")
    
    # Score distribution by label
    print(f"\nScore distribution by label:")
    print(f"  Normal (0):         mean={np.mean(scores[labels==0]):.6e}, std={np.std(scores[labels==0]):.6e}")
    print(f"  Splice-altering (1): mean={np.mean(scores[labels==1]):.6e}, std={np.std(scores[labels==1]):.6e}")
    
    # Plot ROC and PR curves
    plt.figure(figsize=(12, 5))
    
    # ROC Curve
    fpr, tpr, roc_thresholds = roc_curve(labels, scores)
    ax1 = plt.subplot(1, 2, 1)
    ax1.plot(fpr, tpr, label=f'AUROC = {auroc:.4f}', linewidth=2)
    ax1.plot([0, 1], [0, 1], 'k--', label='Random', alpha=0.3)
    ax1.set_xlabel('False Positive Rate', fontsize=12)
    ax1.set_ylabel('True Positive Rate', fontsize=12)
    ax1.set_title(f'ROC Curve - NT Rat ({score_column})', fontsize=13, fontweight='bold')
    ax1.legend(loc='lower right', fontsize=11)
    ax1.grid(alpha=0.3)
    
    # PR Curve
    precision, recall, pr_thresholds = precision_recall_curve(labels, scores)
    ax2 = plt.subplot(1, 2, 2)
    ax2.plot(recall, precision, label=f'AUPRC = {auprc:.4f}', linewidth=2)
    ax2.set_xlabel('Recall', fontsize=12)
    ax2.set_ylabel('Precision', fontsize=12)
    ax2.set_title(f'Precision-Recall Curve - NT Rat ({score_column})', fontsize=13, fontweight='bold')
    ax2.legend(loc='lower left', fontsize=11)
    ax2.grid(alpha=0.3)
    
    plt.tight_layout()
    
    # Save figure
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    fig_path = output_dir / f'nt_rat_{score_column}_curves.png'
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    print(f"\nSaved performance curves to: {fig_path}")
    plt.close()
    
    # Calculate optimal threshold (Youden's J)
    j_scores = tpr - fpr
    optimal_idx = np.argmax(j_scores)
    optimal_threshold = roc_thresholds[optimal_idx]
    
    y_pred_optimal = (scores >= optimal_threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(labels, y_pred_optimal).ravel()
    
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
    ppv = tp / (tp + fp) if (tp + fp) > 0 else 0
    npv = tn / (tn + fn) if (tn + fn) > 0 else 0
    accuracy = accuracy_score(labels, y_pred_optimal)
    f1 = f1_score(labels, y_pred_optimal)
    
    print(f"\nOptimal threshold (Youden's J): {optimal_threshold:.6e}")
    print(f"Confusion Matrix:")
    print(f"  TN: {tn:,}  FP: {fp:,}")
    print(f"  FN: {fn:,}  TP: {tp:,}")
    print(f"\nClassification metrics at optimal threshold:")
    print(f"  Accuracy:    {accuracy:.4f}")
    print(f"  Sensitivity: {sensitivity:.4f}")
    print(f"  Specificity: {specificity:.4f}")
    print(f"  PPV:         {ppv:.4f}")
    print(f"  NPV:         {npv:.4f}")
    print(f"  F1 Score:    {f1:.4f}")
    
    # Save detailed results
    results = {
        'score_column': score_column,
        'n_variants': int(len(scores)),
        'n_positive': int(np.sum(labels)),
        'n_negative': int(len(labels) - np.sum(labels)),
        'positive_rate': float(np.mean(labels)),
        'score_stats': {
            'mean': float(np.mean(scores)),
            'std': float(np.std(scores)),
            'min': float(np.min(scores)),
            'max': float(np.max(scores))
        },
        'correlation': {
            'pearson_r': float(corr),
            'p_value': float(pval)
        },
        'auroc': float(auroc),
        'auprc': float(auprc),
        'optimal_threshold': float(optimal_threshold),
        'confusion_matrix': {
            'true_negatives': int(tn),
            'false_positives': int(fp),
            'false_negatives': int(fn),
            'true_positives': int(tp)
        },
        'classification_metrics': {
            'accuracy': float(accuracy),
            'sensitivity': float(sensitivity),
            'specificity': float(specificity),
            'ppv': float(ppv),
            'npv': float(npv),
            'f1_score': float(f1)
        }
    }
    
    results_path = output_dir / f'nt_rat_{score_column}_results.json'
    with open(results_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved detailed results to: {results_path}")
    
    return results

def main():
    parser = argparse.ArgumentParser(description="Evaluate Nucleotide Transformer on rat data")
    parser.add_argument('--scores_csv', type=str,
                        default='results/nt_rat_scores.csv',
                        help='Path to NT prediction scores CSV')
    parser.add_argument('--sequences_tsv', type=str,
                        default='rat_sequences_8192bp.tsv',
                        help='Path to rat sequences TSV with labels')
    parser.add_argument('--output_dir', type=str,
                        default='results',
                        help='Directory to save evaluation results')
    
    args = parser.parse_args()
    
    print("="*60)
    print("Nucleotide Transformer Rat Evaluation")
    print("="*60)
    
    # Load data
    df_scores = load_nt_scores(args.scores_csv)
    df_labels = load_rat_labels(args.sequences_tsv)
    
    # Merge
    df_merged = merge_predictions_with_labels(df_scores, df_labels)
    if len(df_merged) == 0:
        print("ERROR: No overlapping variants after merging. Exiting.")
        return
    
    # Save merged data
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    merged_path = output_dir / 'nt_rat_merged.csv'
    df_merged.to_csv(merged_path, index=False)
    print(f"\nSaved merged data to: {merged_path}")
    
    # Evaluate both score types
    all_results = {}
    
    # Cosine distance (primary metric used in human evaluation)
    print("\n" + "="*60)
    print("PRIMARY EVALUATION: Cosine Distance")
    print("="*60)
    results_cosine = evaluate_performance(df_merged, 'score_cosine', args.output_dir)
    if results_cosine:
        all_results['cosine'] = results_cosine
    
    # L2 distance (secondary metric)
    print("\n" + "="*60)
    print("SECONDARY EVALUATION: L2 Distance")
    print("="*60)
    results_l2 = evaluate_performance(df_merged, 'score_l2', args.output_dir)
    if results_l2:
        all_results['l2'] = results_l2
    
    # Save comparison
    if all_results:
        comparison_path = output_dir / 'nt_rat_comparison.json'
        with open(comparison_path, 'w') as f:
            json.dump(all_results, f, indent=2)
        print(f"\nSaved score comparison to: {comparison_path}")
        
        print(f"\n{'='*60}")
        print(f"SCORE COMPARISON SUMMARY")
        print(f"{'='*60}")
        print(f"{'Metric':<15} {'Cosine':<12} {'L2':<12}")
        print(f"{'-'*40}")
        print(f"{'AUROC':<15} {all_results['cosine']['auroc']:.4f}       {all_results['l2']['auroc']:.4f}")
        print(f"{'AUPRC':<15} {all_results['cosine']['auprc']:.4f}       {all_results['l2']['auprc']:.4f}")
        print(f"{'='*60}")
        
        # Comparison with human results
        print(f"\n{'='*60}")
        print(f"CROSS-SPECIES COMPARISON")
        print(f"{'='*60}")
        print(f"Nucleotide Transformer Performance:")
        print(f"  Human (SpliceVarDB): AUROC 0.5073, AUPRC 0.5319")
        print(f"  Rat (ratGTEx):       AUROC {all_results['cosine']['auroc']:.4f}, AUPRC {all_results['cosine']['auprc']:.4f}")
        print(f"\nConclusion:")
        if abs(all_results['cosine']['auroc'] - 0.5073) < 0.02:
            print("  ✓ Consistent near-random performance across species")
            print("  → Confirms NT's failure is task-specific, not species-specific")
        else:
            print("  ⚠ Performance differs from human results")
        print(f"{'='*60}")

if __name__ == '__main__':
    main()
