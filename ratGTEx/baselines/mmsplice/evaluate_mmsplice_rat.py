#!/usr/bin/env python3
"""
Evaluate MMSplice performance on ratGTEx dataset.

Following the same methodology as SpliceVarDB evaluation:
- Aggregate multi-exon/transcript predictions per variant (max absolute effect)
- Evaluate pathogenicity and absolute delta_logit_psi
- Compare with SpliceAI (AUROC 0.57), SpliceTransformer (AUROC 0.50) on rat
"""

import pandas as pd
import numpy as np
from sklearn.metrics import (
    roc_auc_score,
    average_precision_score,
    roc_curve,
    precision_recall_curve,
    confusion_matrix
)
import matplotlib.pyplot as plt
import seaborn as sns
import json
import argparse
from pathlib import Path


def load_mmsplice_predictions(predictions_csv):
    """
    Load and aggregate MMSplice predictions.
    
    Args:
        predictions_csv: Path to MMSplice output CSV
        
    Returns:
        DataFrame with aggregated predictions per variant
    """
    print(f"Loading MMSplice predictions from: {predictions_csv}")
    df = pd.read_csv(predictions_csv)
    print(f"Loaded {len(df):,} predictions for {df['ID'].nunique():,} unique variants")
    
    # Aggregate multiple predictions per variant
    # Following MMSplice best practices: use maximum absolute effect across all exons
    print("\nAggregating predictions (taking max absolute value per variant)...")
    
    agg_dict = {
        'delta_logit_psi': lambda x: x.loc[x.abs().idxmax()],  # Max absolute value
        'pathogenicity': 'max',  # Max pathogenicity score
        'efficiency': lambda x: x.loc[x.abs().idxmax()],  # Max absolute value
        'region': 'first',  # Keep first region type for reference
        'gene_name': 'first'  # Keep first gene name
    }
    
    df_agg = df.groupby('ID').agg(agg_dict).reset_index()
    print(f"Aggregated to {len(df_agg):,} unique variants")
    
    return df_agg


def load_ratgtex_labels(ratgtex_tsv):
    """
    Load ratGTEx ground truth labels.
    
    Args:
        ratgtex_tsv: Path to ratGTEx TSV file
        
    Returns:
        DataFrame with variant_id and label columns
    """
    print(f"\nLoading ratGTEx labels from: {ratgtex_tsv}")
    # RatGTEx format: space-separated, no header
    # Column 1: variant_id (e.g., '15:55290237_A/G')
    # Column 2: ref_seq
    # Column 3: alt_seq
    # Column 4: label (0 or 1)
    # Column 5: chrom
    df = pd.read_csv(ratgtex_tsv, sep=r'\s+', header=None,
                     names=['variant_id_orig', 'ref_seq', 'alt_seq', 'label', 'chrom'])
    print(f"Loaded {len(df):,} labeled variants")
    
    # Convert variant_id format from '15:55290237_A/G' to '15:55290237:A>G'
    # to match MMSplice output format (CHR:POS:REF>ALT)
    df['variant_id'] = df['variant_id_orig'].str.replace('_', ':', n=1).str.replace('/', '>')
    
    # Keep only variant_id and label
    df_labels = df[['variant_id', 'label']].copy()
    
    # Label should already be binary (0/1), verify
    unique_labels = df_labels['label'].unique()
    print(f"Label values: {unique_labels}")
    if not set(unique_labels).issubset({0, 1}):
        raise ValueError(f"Expected binary labels (0/1), got: {unique_labels}")
    
    print(f"Created {len(df_labels):,} variant IDs")
    print(f"Label distribution: 0={np.sum(df_labels['label']==0):,}, 1={np.sum(df_labels['label']==1):,}")
    
    return df_labels


def merge_predictions_with_labels(df_predictions, df_labels):
    """
    Merge MMSplice predictions with ground truth labels.
    
    Args:
        df_predictions: Aggregated MMSplice predictions
        df_labels: ratGTEx labels
        
    Returns:
        Merged DataFrame
    """
    print("\nMerging predictions with labels...")
    
    # Rename ID column to variant_id for merging
    df_predictions = df_predictions.rename(columns={'ID': 'variant_id'})
    
    # Merge
    df_merged = df_predictions.merge(df_labels, on='variant_id', how='inner')
    
    print(f"Successfully merged {len(df_merged):,} variants")
    n_pos = int(df_merged['label'].sum())
    n_neg = int(len(df_merged) - n_pos)
    pos_rate = df_merged['label'].mean() * 100
    neg_rate = (1 - df_merged['label'].mean()) * 100
    print(f"Positive samples: {n_pos:,} ({pos_rate:.1f}%)")
    print(f"Negative samples: {n_neg:,} ({neg_rate:.1f}%)")
    
    # Check for missing variants
    missing_count = len(df_labels) - len(df_merged)
    if missing_count > 0:
        print(f"Warning: {missing_count} variants in ratGTEx not found in MMSplice predictions")
    
    return df_merged


def evaluate_performance(df_merged, score_column, output_dir, use_absolute=False):
    """
    Evaluate MMSplice performance.
    
    Args:
        df_merged: Merged predictions and labels
        score_column: Name of score column to evaluate
        output_dir: Output directory for results
        use_absolute: If True, use absolute value of scores
    """
    suffix = "_abs" if use_absolute else ""
    print(f"\n{'='*60}")
    print(f"Evaluating performance using: {score_column}{suffix}")
    print(f"{'='*60}")
    
    # Extract scores and labels
    scores = df_merged[score_column].values
    labels = df_merged['label'].values
    
    # Use absolute value if specified
    if use_absolute:
        print(f"Using absolute value (magnitude of change)")
        scores = np.abs(scores)
    
    # Handle NaN values
    valid_mask = ~np.isnan(scores)
    if not valid_mask.all():
        n_nan = (~valid_mask).sum()
        print(f"Warning: {n_nan} NaN scores found, excluding from evaluation")
        scores = scores[valid_mask]
        labels = labels[valid_mask]
    
    # Calculate metrics
    auroc = roc_auc_score(labels, scores)
    auprc = average_precision_score(labels, scores)
    
    print(f"\n{'='*60}")
    print(f"PERFORMANCE METRICS ({score_column}{suffix})")
    print(f"{'='*60}")
    print(f"Total variants: {len(scores):,}")
    print(f"Positive samples: {np.sum(labels):,} ({np.mean(labels)*100:.1f}%)")
    print(f"Negative samples: {len(labels) - np.sum(labels):,} ({(1-np.mean(labels))*100:.1f}%)")
    print(f"\nScore range: [{np.min(scores):.4f}, {np.max(scores):.4f}]")
    print(f"Score mean: {np.mean(scores):.4f}")
    print(f"Score std: {np.std(scores):.4f}")
    print(f"\n{'='*60}")
    print(f"AUROC: {auroc:.4f}")
    print(f"AUPRC: {auprc:.4f}")
    print(f"{'='*60}")
    
    # ROC and PR curves
    fpr, tpr, roc_thresholds = roc_curve(labels, scores)
    precision, recall, pr_thresholds = precision_recall_curve(labels, scores)
    
    # Create visualization
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # ROC curve
    axes[0].plot(fpr, tpr, color='darkblue', lw=2, label=f'ROC (AUC = {auroc:.4f})')
    axes[0].plot([0, 1], [0, 1], 'k--', lw=1, label='Random')
    axes[0].set_xlabel('False Positive Rate', fontsize=12)
    axes[0].set_ylabel('True Positive Rate', fontsize=12)
    axes[0].set_title(f'ROC Curve - MMSplice Rat ({score_column}{suffix})', fontsize=14, fontweight='bold')
    axes[0].legend(loc='lower right', fontsize=11)
    axes[0].grid(alpha=0.3)
    
    # PR curve
    axes[1].plot(recall, precision, color='darkgreen', lw=2, label=f'PR (AUC = {auprc:.4f})')
    baseline = np.mean(labels)
    axes[1].axhline(baseline, color='k', linestyle='--', lw=1, label=f'Baseline ({baseline:.3f})')
    axes[1].set_xlabel('Recall', fontsize=12)
    axes[1].set_ylabel('Precision', fontsize=12)
    axes[1].set_title(f'Precision-Recall Curve - MMSplice Rat ({score_column}{suffix})', fontsize=14, fontweight='bold')
    axes[1].legend(loc='best', fontsize=11)
    axes[1].grid(alpha=0.3)
    
    plt.tight_layout()
    
    # Save figure
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    fig_path = output_dir / f'mmsplice_rat_{score_column}{suffix}_curves.png'
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    print(f"\nSaved performance curves to: {fig_path}")
    plt.close()
    
    # Calculate confusion matrix at optimal threshold (Youden's J)
    j_scores = tpr - fpr
    optimal_idx = np.argmax(j_scores)
    optimal_threshold = roc_thresholds[optimal_idx]
    
    y_pred = (scores >= optimal_threshold).astype(int)
    cm = confusion_matrix(labels, y_pred)
    
    tn, fp, fn, tp = cm.ravel()
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
    ppv = tp / (tp + fp) if (tp + fp) > 0 else 0
    npv = tn / (tn + fn) if (tn + fn) > 0 else 0
    
    print(f"\nOptimal threshold (Youden's J): {optimal_threshold:.4f}")
    print(f"Confusion Matrix:")
    print(f"  TN: {tn:,}  FP: {fp:,}")
    print(f"  FN: {fn:,}  TP: {tp:,}")
    print(f"\nSensitivity: {sensitivity:.4f}")
    print(f"Specificity: {specificity:.4f}")
    print(f"PPV: {ppv:.4f}")
    print(f"NPV: {npv:.4f}")
    
    # Save results
    results = {
        'score_column': score_column + suffix,
        'use_absolute': use_absolute,
        'n_variants': int(len(scores)),
        'n_positive': int(np.sum(labels)),
        'n_negative': int(len(labels) - np.sum(labels)),
        'positive_rate': float(np.mean(labels)),
        'auroc': float(auroc),
        'auprc': float(auprc),
        'score_range': [float(np.min(scores)), float(np.max(scores))],
        'score_mean': float(np.mean(scores)),
        'score_std': float(np.std(scores)),
        'optimal_threshold': float(optimal_threshold),
        'confusion_matrix': {
            'tn': int(tn),
            'fp': int(fp),
            'fn': int(fn),
            'tp': int(tp)
        },
        'sensitivity': float(sensitivity),
        'specificity': float(specificity),
        'ppv': float(ppv),
        'npv': float(npv)
    }
    
    results_path = output_dir / f'mmsplice_rat_{score_column}{suffix}_results.json'
    with open(results_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved detailed results to: {results_path}")
    
    return results


def main():
    parser = argparse.ArgumentParser(
        description="Evaluate MMSplice performance on ratGTEx"
    )
    parser.add_argument(
        '--predictions',
        type=str,
        default='results/mmsplice_ratgtex_scores.csv',
        help='Path to MMSplice predictions CSV'
    )
    parser.add_argument(
        '--labels',
        type=str,
        default='/Users/byronsun/Desktop/AS_复现模型/BIB_review/data/processed_data/ratGTEx/ratgtex_silver_benchmark_balanced.tsv',
        help='Path to ratGTEx labels TSV'
    )
    parser.add_argument(
        '--output_dir',
        type=str,
        default='results',
        help='Output directory for results'
    )
    
    args = parser.parse_args()
    
    print("="*60)
    print("MMSplice ratGTEx Evaluation")
    print("="*60)
    
    # Load data
    df_predictions = load_mmsplice_predictions(args.predictions)
    df_labels = load_ratgtex_labels(args.labels)
    
    # Merge
    df_merged = merge_predictions_with_labels(df_predictions, df_labels)
    
    # Save merged data
    merged_path = Path(args.output_dir) / 'mmsplice_ratgtex_merged.csv'
    df_merged.to_csv(merged_path, index=False)
    print(f"\nSaved merged data to: {merged_path}")
    
    # Evaluate each score type
    # For delta_logit_psi and efficiency: use absolute value (magnitude matters)
    # For pathogenicity: use original value (higher = more pathogenic)
    all_results = {}
    
    # Pathogenicity (original)
    results = evaluate_performance(df_merged, 'pathogenicity', args.output_dir, use_absolute=False)
    all_results['pathogenicity'] = results
    
    # Delta_logit_psi (absolute value - magnitude matters)
    results = evaluate_performance(df_merged, 'delta_logit_psi', args.output_dir, use_absolute=True)
    all_results['delta_logit_psi_abs'] = results
    
    # Efficiency (absolute value - magnitude matters)
    results = evaluate_performance(df_merged, 'efficiency', args.output_dir, use_absolute=True)
    all_results['efficiency_abs'] = results
    
    # Save summary comparison
    summary_path = Path(args.output_dir) / 'mmsplice_rat_scores_comparison.json'
    with open(summary_path, 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f"\n{'='*60}")
    print(f"Saved score comparison to: {summary_path}")
    print(f"{'='*60}")
    
    # Print summary table
    print(f"\n{'='*60}")
    print("SCORE COMPARISON SUMMARY (Rat)")
    print(f"{'='*60}")
    print(f"{'Score':<25} {'AUROC':<10} {'AUPRC':<10}")
    print(f"{'-'*45}")
    for score_col, res in all_results.items():
        print(f"{score_col:<25} {res['auroc']:<10.4f} {res['auprc']:<10.4f}")
    print(f"{'='*60}")
    
    # Comparison with other models
    print(f"\n{'='*60}")
    print("CROSS-MODEL COMPARISON (Rat)")
    print(f"{'='*60}")
    print(f"{'Model':<30} {'AUROC':<10}")
    print(f"{'-'*40}")
    print(f"{'SpliceAI':<30} {'0.5696':<10}")
    print(f"{'SpliceTransformer':<30} {'0.5011':<10}")
    print(f"{'Pangolin':<30} {'0.5100':<10}")
    print(f"{'MMSplice (pathogenicity)':<30} {all_results['pathogenicity']['auroc']:<10.4f}")
    print(f"{'MMSplice (delta_logit_psi_abs)':<30} {all_results['delta_logit_psi_abs']['auroc']:<10.4f}")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
