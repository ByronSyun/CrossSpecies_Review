#!/usr/bin/env python3
"""
Evaluate SpliceBERT variant effect scores against ground truth labels
"""

import argparse
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import (
    roc_auc_score, average_precision_score, accuracy_score,
    precision_score, recall_score, f1_score, confusion_matrix,
    precision_recall_curve, roc_curve
)
import os

def load_and_merge_data(labels_file, scores_file):
    """Load labels and scores, merge on variant_id"""
    
    # Load ground truth labels
    df_labels = pd.read_csv(labels_file, sep='\t')
    print(f"Loaded {len(df_labels)} labeled variants")
    
    # Load SpliceBERT scores
    df_scores = pd.read_csv(scores_file)
    print(f"Loaded {len(df_scores)} scored variants")
    
    # Prepare labels
    if 'label' in df_labels.columns:
        df_labels['y'] = (df_labels['label'] == 'splice-altering').astype(int)
    elif 'is_splice_altering' in df_labels.columns:
        df_labels['y'] = df_labels['is_splice_altering'].astype(int)
    else:
        raise ValueError("No suitable label column found. Expected 'label' or 'is_splice_altering'")
    
    # Create variant ID for matching - use original variant_id if available
    if 'variant_id' in df_labels.columns:
        df_labels['vid'] = df_labels['variant_id']
    elif {'chrom', 'pos_0based', 'ref', 'alt'}.issubset(df_labels.columns):
        chrom = df_labels['chrom'].astype(str).str.replace('^chr', '', regex=True)
        df_labels['vid'] = chrom + ':' + (df_labels['pos_0based'] + 1).astype(str) + '_' + df_labels['ref'] + '/' + df_labels['alt']
    else:
        raise ValueError("No suitable columns found for creating variant IDs")
    
    df_scores['vid'] = df_scores['variant_id']
    
    # Merge datasets
    merged = df_labels[['vid', 'y']].merge(df_scores[['vid', 'kl_context_score']], on='vid', how='inner')
    print(f"Successfully merged {len(merged)} variants")
    print(f"Positive samples: {merged['y'].sum()}, Negative samples: {(merged['y'] == 0).sum()}")
    
    return merged

def compute_metrics(y_true, y_scores, threshold=0.0):
    """Compute classification metrics at given threshold"""
    y_pred = (y_scores >= threshold).astype(int)
    
    metrics = {
        'threshold': threshold,
        'auroc': roc_auc_score(y_true, y_scores),
        'auprc': average_precision_score(y_true, y_scores),
        'accuracy': accuracy_score(y_true, y_pred),
        'precision': precision_score(y_true, y_pred, zero_division=0),
        'recall': recall_score(y_true, y_pred, zero_division=0),
        'f1': f1_score(y_true, y_pred, zero_division=0),
        'specificity': recall_score(1 - y_true, 1 - y_pred, zero_division=0)
    }
    
    # Confusion matrix
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    metrics.update({
        'true_negatives': int(tn),
        'false_positives': int(fp),
        'false_negatives': int(fn),
        'true_positives': int(tp)
    })
    
    return metrics

def find_optimal_threshold(y_true, y_scores):
    """Find F1-optimal threshold"""
    precision_vals, recall_vals, thresholds = precision_recall_curve(y_true, y_scores)
    f1_scores = 2 * (precision_vals * recall_vals) / (precision_vals + recall_vals + 1e-8)
    
    # Find threshold that maximizes F1
    optimal_idx = np.argmax(f1_scores)
    optimal_threshold = thresholds[optimal_idx] if optimal_idx < len(thresholds) else thresholds[-1]
    
    return optimal_threshold, f1_scores[optimal_idx]

def create_visualizations(y_true, y_scores, out_dir):
    """Create performance visualizations"""
    os.makedirs(out_dir, exist_ok=True)
    
    # Set style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # 1. Score distribution
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Score distribution by class
    pos_scores = y_scores[y_true == 1]
    neg_scores = y_scores[y_true == 0]
    
    ax1.hist(neg_scores, bins=50, alpha=0.7, label=f'Normal (n={len(neg_scores)})', density=True)
    ax1.hist(pos_scores, bins=50, alpha=0.7, label=f'Splice-altering (n={len(pos_scores)})', density=True)
    ax1.set_xlabel('SpliceBERT KL-context Score')
    ax1.set_ylabel('Density')
    ax1.set_title('Score Distribution by Class')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Box plot
    data_for_box = [neg_scores, pos_scores]
    ax2.boxplot(data_for_box, labels=['Normal', 'Splice-altering'])
    ax2.set_ylabel('SpliceBERT KL-context Score')
    ax2.set_title('Score Distribution by Class')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{out_dir}/score_distribution.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. ROC Curve
    fpr, tpr, _ = roc_curve(y_true, y_scores)
    auroc = roc_auc_score(y_true, y_scores)
    
    plt.figure(figsize=(8, 6))
    plt.plot(fpr, tpr, linewidth=2, label=f'SpliceBERT (AUROC = {auroc:.3f})')
    plt.plot([0, 1], [0, 1], 'k--', alpha=0.5, label='Random')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curve')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(f"{out_dir}/roc_curve.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. Precision-Recall Curve
    precision_vals, recall_vals, _ = precision_recall_curve(y_true, y_scores)
    auprc = average_precision_score(y_true, y_scores)
    
    plt.figure(figsize=(8, 6))
    plt.plot(recall_vals, precision_vals, linewidth=2, label=f'SpliceBERT (AUPRC = {auprc:.3f})')
    plt.axhline(y=np.mean(y_true), color='k', linestyle='--', alpha=0.5, 
                label=f'Random (baseline = {np.mean(y_true):.3f})')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Precision-Recall Curve')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(f"{out_dir}/precision_recall_curve.png", dpi=300, bbox_inches='tight')
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Evaluate SpliceBERT variant effect scores")
    parser.add_argument('--labels', required=True, help='Ground truth labels file (TSV)')
    parser.add_argument('--scores', required=True, help='SpliceBERT scores file (CSV)')
    parser.add_argument('--out_dir', required=True, help='Output directory for results')
    
    args = parser.parse_args()
    
    print("=== SpliceBERT Score Evaluation ===")
    
    # Load and merge data
    merged = load_and_merge_data(args.labels, args.scores)
    y_true = merged['y'].values
    y_scores = merged['kl_context_score'].values
    
    # Handle NaN values
    nan_mask = np.isnan(y_scores)
    n_nan = np.sum(nan_mask)
    print(f"\nNaN handling:")
    print(f"  Total variants: {len(y_scores)}")
    print(f"  NaN scores: {n_nan} ({n_nan/len(y_scores)*100:.1f}%)")
    
    if n_nan > 0:
        # Remove NaN values
        valid_mask = ~nan_mask
        y_scores = y_scores[valid_mask]
        y_true = y_true[valid_mask]
        print(f"  Valid scores for evaluation: {len(y_scores)}")
    
    # SpliceBERT uses KL-divergence: higher values indicate MORE difference (more splicing-altering)
    # But our scores are negative, so we need to invert them
    y_scores = -y_scores
    print(f"  Inverted scores (higher = more splice-altering)")
    
    # Basic score statistics
    print(f"\nScore Statistics:")
    if len(y_scores) > 0:
        print(f"  Range: [{np.min(y_scores):.4f}, {np.max(y_scores):.4f}]")
        print(f"  Mean: {np.mean(y_scores):.4f}")
        print(f"  Std: {np.std(y_scores):.4f}")
        print(f"  Median: {np.median(y_scores):.4f}")
    else:
        print("  No valid scores available!")
        return
    
    # Find optimal threshold
    optimal_thresh, optimal_f1 = find_optimal_threshold(y_true, y_scores)
    print(f"\nOptimal F1 threshold: {optimal_thresh:.4f} (F1 = {optimal_f1:.4f})")
    
    # Evaluate at multiple thresholds
    thresholds = [0.0, optimal_thresh]
    
    # Add some additional thresholds based on score distribution
    score_percentiles = np.percentile(y_scores, [25, 50, 75, 90, 95])
    thresholds.extend(score_percentiles)
    thresholds = sorted(list(set(thresholds)))  # Remove duplicates and sort
    
    print(f"\nEvaluating at {len(thresholds)} thresholds...")
    
    results = []
    for thresh in thresholds:
        metrics = compute_metrics(y_true, y_scores, thresh)
        results.append(metrics)
        
        if thresh == 0.0 or thresh == optimal_thresh:
            print(f"\nThreshold {thresh:.4f}:")
            print(f"  AUROC: {metrics['auroc']:.4f}")
            print(f"  AUPRC: {metrics['auprc']:.4f}")
            print(f"  Accuracy: {metrics['accuracy']:.4f}")
            print(f"  Precision: {metrics['precision']:.4f}")
            print(f"  Recall: {metrics['recall']:.4f}")
            print(f"  F1: {metrics['f1']:.4f}")
            print(f"  Specificity: {metrics['specificity']:.4f}")
    
    # Save results
    os.makedirs(args.out_dir, exist_ok=True)
    
    # Save detailed results
    results_df = pd.DataFrame(results)
    results_df.to_csv(f"{args.out_dir}/splicebert_evaluation_results.csv", index=False)
    
    # Save summary JSON
    summary = {
        'model': 'SpliceBERT',
        'dataset_size': len(merged),
        'positive_samples': int(np.sum(y_true)),
        'negative_samples': int(len(y_true) - np.sum(y_true)),
        'auroc': float(results[0]['auroc']),  # At threshold 0.0
        'auprc': float(results[0]['auprc']),
        'optimal_threshold': float(optimal_thresh),
        'optimal_f1': float(optimal_f1),
        'score_range': [float(np.min(y_scores)), float(np.max(y_scores))],
        'score_mean': float(np.mean(y_scores)),
        'score_std': float(np.std(y_scores))
    }
    
    with open(f"{args.out_dir}/splicebert_evaluation_summary.json", 'w') as f:
        json.dump(summary, f, indent=2)
    
    # Create visualizations
    print(f"\nCreating visualizations...")
    create_visualizations(y_true, y_scores, args.out_dir)
    
    print(f"\n✅ Evaluation complete!")
    print(f"Results saved to: {args.out_dir}")
    print(f"  - Detailed results: splicebert_evaluation_results.csv")
    print(f"  - Summary: splicebert_evaluation_summary.json")
    print(f"  - Visualizations: score_distribution.png, roc_curve.png, precision_recall_curve.png")

if __name__ == '__main__':
    main()
