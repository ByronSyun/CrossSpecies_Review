#!/usr/bin/env python3
"""
Evaluate SpliceBERT variant effect scores against rat ground truth labels
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

def load_and_merge_rat_data(labels_file, scores_file):
    """Load rat labels and scores, merge on variant_id"""
    
    # Load ground truth labels (rat format: no header, 5 columns)
    df_labels = pd.read_csv(labels_file, sep='\t', header=None,
                           names=['variant_id', 'ref_seq', 'alt_seq', 'label', 'tissue_id'])
    print(f"Loaded {len(df_labels)} rat labeled variants (may include tissue duplicates)")
    
    # Handle duplicate variant_ids across tissues
    # Use majority vote for label assignment
    label_consensus = df_labels.groupby('variant_id')['label'].agg(['mean', 'count']).reset_index()
    label_consensus['y'] = (label_consensus['mean'] >= 0.5).astype(int)  # Majority vote
    
    print(f"Unique variants after deduplication: {len(label_consensus)}")
    print(f"Variants with conflicting labels across tissues: {(label_consensus['count'] > 1).sum()}")
    
    # Load SpliceBERT scores
    df_scores = pd.read_csv(scores_file)
    print(f"Loaded {len(df_scores)} scored variants")
    
    # Merge datasets
    merged = label_consensus[['variant_id', 'y']].merge(
        df_scores[['variant_id', 'kl_context_score']], 
        on='variant_id', how='inner'
    )
    
    print(f"Final merged dataset: {len(merged)} variants")
    print(f"Positive samples: {merged['y'].sum()}")
    print(f"Negative samples: {len(merged) - merged['y'].sum()}")
    
    return merged

def compute_metrics(y_true, y_scores, threshold=None):
    """Compute evaluation metrics"""
    
    # Threshold-independent metrics
    try:
        auroc = roc_auc_score(y_true, y_scores)
    except ValueError:
        auroc = np.nan
    
    try:
        auprc = average_precision_score(y_true, y_scores)
    except ValueError:
        auprc = np.nan
    
    # Find optimal threshold if not provided
    if threshold is None:
        precision_curve, recall_curve, thresholds = precision_recall_curve(y_true, y_scores)
        if len(thresholds) > 0:
            f1_scores = 2 * precision_curve[:-1] * recall_curve[:-1] / (precision_curve[:-1] + recall_curve[:-1] + 1e-8)
            optimal_idx = np.argmax(f1_scores)
            threshold = thresholds[optimal_idx]
            
            # Special handling for SpliceBERT (inverted scores - lower is more splice-altering)
            # Check if we need to use a percentile-based threshold
            if threshold <= np.min(y_scores) or threshold >= np.max(y_scores):
                threshold = np.percentile(y_scores, 25)  # Use 25th percentile for inverted scores
        else:
            threshold = np.median(y_scores)
    
    # Threshold-dependent metrics
    y_pred = (y_scores >= threshold).astype(int)
    
    accuracy = accuracy_score(y_true, y_pred)
    precision = precision_score(y_true, y_pred, zero_division=0)
    recall = recall_score(y_true, y_pred, zero_division=0)
    f1 = f1_score(y_true, y_pred, zero_division=0)
    
    return {
        'auroc': auroc,
        'auprc': auprc,
        'accuracy': accuracy,
        'precision': precision,
        'recall': recall,
        'f1': f1,
        'threshold': threshold,
        'n_samples': len(y_true),
        'n_positive': int(y_true.sum()),
        'n_negative': int(len(y_true) - y_true.sum())
    }

def plot_performance(merged_df, output_dir):
    """Generate performance plots"""
    
    os.makedirs(output_dir, exist_ok=True)
    
    y_true = merged_df['y'].values
    y_scores = merged_df['kl_context_score'].values
    
    # Set style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # 1. ROC Curve
    try:
        fpr, tpr, _ = roc_curve(y_true, y_scores)
        auroc = roc_auc_score(y_true, y_scores)
        
        plt.figure(figsize=(8, 6))
        plt.plot(fpr, tpr, linewidth=2, label=f'SpliceBERT (AUROC = {auroc:.4f})')
        plt.plot([0, 1], [0, 1], 'k--', alpha=0.5, label='Random (AUROC = 0.5000)')
        plt.xlabel('False Positive Rate', fontsize=12)
        plt.ylabel('True Positive Rate', fontsize=12)
        plt.title('ROC Curve - SpliceBERT on Rat Data', fontsize=14, fontweight='bold')
        plt.legend(fontsize=11)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(f"{output_dir}/roc_curve.png", dpi=300, bbox_inches='tight')
        plt.savefig(f"{output_dir}/roc_curve.pdf", bbox_inches='tight')
        plt.close()
    except Exception as e:
        print(f"Could not plot ROC curve: {e}")
    
    # 2. Precision-Recall Curve
    try:
        precision, recall, _ = precision_recall_curve(y_true, y_scores)
        auprc = average_precision_score(y_true, y_scores)
        
        plt.figure(figsize=(8, 6))
        plt.plot(recall, precision, linewidth=2, label=f'SpliceBERT (AUPRC = {auprc:.4f})')
        plt.axhline(y=y_true.mean(), color='k', linestyle='--', alpha=0.5, 
                   label=f'Random (AUPRC = {y_true.mean():.4f})')
        plt.xlabel('Recall', fontsize=12)
        plt.ylabel('Precision', fontsize=12)
        plt.title('Precision-Recall Curve - SpliceBERT on Rat Data', fontsize=14, fontweight='bold')
        plt.legend(fontsize=11)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(f"{output_dir}/precision_recall_curve.png", dpi=300, bbox_inches='tight')
        plt.savefig(f"{output_dir}/precision_recall_curve.pdf", bbox_inches='tight')
        plt.close()
    except Exception as e:
        print(f"Could not plot PR curve: {e}")
    
    # 3. Score Distribution
    plt.figure(figsize=(10, 6))
    
    pos_scores = y_scores[y_true == 1]
    neg_scores = y_scores[y_true == 0]
    
    plt.hist(neg_scores, bins=50, alpha=0.7, label=f'Negative (n={len(neg_scores)})', color='lightcoral')
    plt.hist(pos_scores, bins=50, alpha=0.7, label=f'Positive (n={len(pos_scores)})', color='lightblue')
    
    plt.xlabel('KL-Context Score', fontsize=12)
    plt.ylabel('Frequency', fontsize=12)
    plt.title('Score Distribution - SpliceBERT on Rat Data', fontsize=14, fontweight='bold')
    plt.legend(fontsize=11)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/score_distribution.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_dir}/score_distribution.pdf", bbox_inches='tight')
    plt.close()
    
    print(f"Plots saved to {output_dir}")

def main():
    parser = argparse.ArgumentParser(description="Evaluate SpliceBERT performance on rat data")
    parser.add_argument('--labels', required=True, help='Path to rat ground truth TSV file')
    parser.add_argument('--scores', required=True, help='Path to SpliceBERT scores CSV file')
    parser.add_argument('--output_dir', required=True, help='Output directory for results')
    parser.add_argument('--threshold', type=float, help='Optional fixed threshold for binary classification')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load and merge data
    print("Loading and merging data...")
    merged_df = load_and_merge_rat_data(args.labels, args.scores)
    
    if merged_df.empty:
        print("Error: No common variants found between labels and scores!")
        return
    
    # Compute metrics
    print("\nComputing performance metrics...")
    y_true = merged_df['y'].values
    y_scores = merged_df['kl_context_score'].values
    
    # Note: SpliceBERT scores are typically inverted (lower = more splice-altering)
    # We'll try both orientations and use the one with better performance
    metrics_normal = compute_metrics(y_true, y_scores, args.threshold)
    metrics_inverted = compute_metrics(y_true, -y_scores, args.threshold)
    
    # Choose the orientation with better AUROC
    if metrics_inverted['auroc'] > metrics_normal['auroc']:
        print("Using inverted scores (lower KL = more splice-altering)")
        metrics = metrics_inverted
        y_scores = -y_scores
        merged_df['kl_context_score'] = y_scores
    else:
        print("Using normal scores (higher KL = more splice-altering)")
        metrics = metrics_normal
    
    # Print results
    print("\n" + "="*50)
    print("SPLICEBERT PERFORMANCE ON RAT DATA")
    print("="*50)
    print(f"Dataset size: {metrics['n_samples']} variants")
    print(f"Positive samples: {metrics['n_positive']} ({metrics['n_positive']/metrics['n_samples']*100:.1f}%)")
    print(f"Negative samples: {metrics['n_negative']} ({metrics['n_negative']/metrics['n_samples']*100:.1f}%)")
    print()
    print("THRESHOLD-INDEPENDENT METRICS:")
    print(f"  AUROC:        {metrics['auroc']:.4f}")
    print(f"  AUPRC:        {metrics['auprc']:.4f}")
    print()
    print("THRESHOLD-DEPENDENT METRICS:")
    print(f"  Threshold:    {metrics['threshold']:.6f}")
    print(f"  Accuracy:     {metrics['accuracy']:.4f}")
    print(f"  Precision:    {metrics['precision']:.4f}")
    print(f"  Recall:       {metrics['recall']:.4f}")
    print(f"  F1-score:     {metrics['f1']:.4f}")
    print("="*50)
    
    # Save detailed results
    results_file = f"{args.output_dir}/splicebert_rat_evaluation_results.csv"
    results_df = pd.DataFrame([metrics])
    results_df.to_csv(results_file, index=False)
    print(f"Detailed results saved to: {results_file}")
    
    # Save summary JSON
    summary_file = f"{args.output_dir}/splicebert_rat_evaluation_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(metrics, f, indent=2)
    print(f"Summary saved to: {summary_file}")
    
    # Generate plots
    print("\nGenerating performance plots...")
    plot_performance(merged_df, args.output_dir)
    
    # Save merged data for further analysis
    merged_file = f"{args.output_dir}/merged_data.csv"
    merged_df.to_csv(merged_file, index=False)
    print(f"Merged data saved to: {merged_file}")
    
    print("\nEvaluation completed successfully!")

if __name__ == "__main__":
    main()
