#!/usr/bin/env python3
"""
Enhanced Pangolin Evaluation Script with Detailed Visualizations
Generates confusion matrix, precision-recall curves, ROC curves, and comprehensive metrics
for Pangolin performance on SpliceVarDB dataset.

Author: AI Assistant
Date: September 2025
"""

import pandas as pd
import numpy as np
import json
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import (
    roc_auc_score, average_precision_score, accuracy_score,
    precision_score, recall_score, f1_score, confusion_matrix,
    roc_curve, precision_recall_curve, classification_report
)
from datetime import datetime
import os
import re
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - [%(levelname)s] - %(message)s')

def parse_pangolin_score(pangolin_field):
    """Extract maximum absolute score from Pangolin field"""
    if pd.isna(pangolin_field) or pangolin_field == '.' or pangolin_field == '':
        return 0.0
    
    # Format: ENSG00000156876.10|-7:0.85|-1:-0.87|Warnings:
    pattern = r'([+-]?\d+):([+-]?\d+\.?\d*)'
    matches = re.findall(pattern, str(pangolin_field))
    
    if not matches:
        return 0.0
    
    max_abs_score = 0.0
    for pos, score_str in matches:
        try:
            score = float(score_str)
            abs_score = abs(score)
            if abs_score > max_abs_score:
                max_abs_score = abs_score
        except ValueError:
            continue
    
    return max_abs_score

def load_splicevardb_labels(labels_file):
    """Load labels from SpliceVarDB TSV"""
    logging.info(f"Loading labels from {labels_file}")
    df = pd.read_csv(labels_file, sep='\t')
    
    # Create variant_id in the same format as VCF (1-based coordinates, no chr prefix)
    chrom_no_prefix = df['chrom'].str.replace('chr', '', regex=False)
    df['vcf_variant_id'] = chrom_no_prefix + ':' + (df['pos_0based'] + 1).astype(str) + '_' + df['ref'] + '/' + df['alt']
    
    # Convert label to binary
    df['binary_label'] = (df['label'] == 'splice-altering').astype(int)
    
    logging.info(f"Loaded {len(df)} variants")
    logging.info(f"Splice-altering: {df['binary_label'].sum()}")
    logging.info(f"Normal: {len(df) - df['binary_label'].sum()}")
    
    return df[['vcf_variant_id', 'binary_label']].set_index('vcf_variant_id')

def load_pangolin_vcf(vcf_file):
    """Load Pangolin results from VCF"""
    logging.info(f"Loading Pangolin results from {vcf_file}")
    results = []
    
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            chrom, pos, _, ref, alt, _, _, info = fields[:8]
            
            # Create variant ID
            variant_id = f"{chrom}:{pos}_{ref}/{alt}"
            
            # Extract Pangolin score
            pangolin_score = 0.0
            if 'Pangolin=' in info:
                pangolin_field = info.split('Pangolin=')[1].split(';')[0]
                pangolin_score = parse_pangolin_score(pangolin_field)
            
            results.append({
                'variant_id': variant_id,
                'pangolin_score': pangolin_score
            })
    
    df = pd.DataFrame(results).set_index('variant_id')
    logging.info(f"Loaded {len(df)} Pangolin predictions")
    logging.info(f"Non-zero scores: {(df['pangolin_score'] > 0).sum()}")
    
    return df

def compute_detailed_metrics(y_true, y_scores, thresholds=[0.001, 0.01, 0.05, 0.1, 0.2]):
    """Compute detailed classification metrics at multiple thresholds"""
    
    # Threshold-independent metrics
    auroc = roc_auc_score(y_true, y_scores)
    auprc = average_precision_score(y_true, y_scores)
    
    # Threshold-dependent metrics
    metrics_by_threshold = {}
    
    for threshold in thresholds:
        y_pred = (y_scores >= threshold).astype(int)
        
        # Basic metrics
        accuracy = accuracy_score(y_true, y_pred)
        precision = precision_score(y_true, y_pred, zero_division=0)
        recall = recall_score(y_true, y_pred, zero_division=0)
        f1 = f1_score(y_true, y_pred, zero_division=0)
        
        # Confusion matrix
        tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
        specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
        
        metrics_by_threshold[f'threshold_{threshold}'] = {
            'threshold': threshold,
            'accuracy': accuracy,
            'precision': precision,
            'recall': recall,
            'f1_score': f1,
            'specificity': specificity,
            'tp': int(tp), 'fp': int(fp), 'tn': int(tn), 'fn': int(fn)
        }
    
    return {
        'auroc': auroc,
        'auprc': auprc,
        'metrics_by_threshold': metrics_by_threshold
    }

def plot_confusion_matrices(y_true, y_scores, thresholds, output_dir):
    """Generate confusion matrix plots for multiple thresholds"""
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.ravel()
    
    for i, threshold in enumerate(thresholds):
        if i >= 6:  # Only plot first 6 thresholds
            break
            
        y_pred = (y_scores >= threshold).astype(int)
        cm = confusion_matrix(y_true, y_pred)
        
        # Plot confusion matrix
        sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', 
                   xticklabels=['Normal', 'Splice-altering'],
                   yticklabels=['Normal', 'Splice-altering'],
                   ax=axes[i])
        axes[i].set_title(f'Threshold = {threshold}')
        axes[i].set_xlabel('Predicted')
        axes[i].set_ylabel('True')
    
    # Remove empty subplots
    for j in range(len(thresholds), 6):
        fig.delaxes(axes[j])
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'confusion_matrices.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'confusion_matrices.pdf'), bbox_inches='tight')
    plt.close()
    
    logging.info("✅ Confusion matrices saved")

def plot_roc_curve(y_true, y_scores, output_dir):
    """Generate ROC curve"""
    
    fpr, tpr, _ = roc_curve(y_true, y_scores)
    auroc = roc_auc_score(y_true, y_scores)
    
    plt.figure(figsize=(8, 6))
    plt.plot(fpr, tpr, linewidth=2, label=f'Pangolin (AUROC = {auroc:.4f})')
    plt.plot([0, 1], [0, 1], 'k--', alpha=0.5, label='Random')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curve - Pangolin on SpliceVarDB')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.savefig(os.path.join(output_dir, 'roc_curve.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'roc_curve.pdf'), bbox_inches='tight')
    plt.close()
    
    logging.info("✅ ROC curve saved")

def plot_precision_recall_curve(y_true, y_scores, output_dir):
    """Generate Precision-Recall curve"""
    
    precision, recall, thresholds = precision_recall_curve(y_true, y_scores)
    auprc = average_precision_score(y_true, y_scores)
    baseline = np.sum(y_true) / len(y_true)  # Random baseline
    
    plt.figure(figsize=(8, 6))
    plt.plot(recall, precision, linewidth=2, label=f'Pangolin (AUPRC = {auprc:.4f})')
    plt.axhline(y=baseline, color='k', linestyle='--', alpha=0.5, 
                label=f'Random (Baseline = {baseline:.4f})')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Precision-Recall Curve - Pangolin on SpliceVarDB')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.savefig(os.path.join(output_dir, 'precision_recall_curve.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'precision_recall_curve.pdf'), bbox_inches='tight')
    plt.close()
    
    logging.info("✅ Precision-Recall curve saved")

def plot_score_distribution(y_true, y_scores, output_dir):
    """Plot score distribution by class"""
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Score distribution histogram
    positive_scores = y_scores[y_true == 1]
    negative_scores = y_scores[y_true == 0]
    
    ax1.hist(negative_scores, bins=50, alpha=0.7, label='Normal', density=True, color='blue')
    ax1.hist(positive_scores, bins=50, alpha=0.7, label='Splice-altering', density=True, color='red')
    ax1.set_xlabel('Pangolin Score')
    ax1.set_ylabel('Density')
    ax1.set_title('Score Distribution by Class')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Box plot
    data_to_plot = [negative_scores, positive_scores]
    ax2.boxplot(data_to_plot, labels=['Normal', 'Splice-altering'])
    ax2.set_ylabel('Pangolin Score')
    ax2.set_title('Score Distribution Box Plot')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'score_distribution.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'score_distribution.pdf'), bbox_inches='tight')
    plt.close()
    
    logging.info("✅ Score distribution plots saved")

def plot_threshold_performance(metrics_data, output_dir):
    """Plot performance metrics across different thresholds"""
    
    thresholds = []
    accuracies = []
    precisions = []
    recalls = []
    f1_scores = []
    specificities = []
    
    for key, metrics in metrics_data['metrics_by_threshold'].items():
        thresholds.append(metrics['threshold'])
        accuracies.append(metrics['accuracy'])
        precisions.append(metrics['precision'])
        recalls.append(metrics['recall'])
        f1_scores.append(metrics['f1_score'])
        specificities.append(metrics['specificity'])
    
    plt.figure(figsize=(12, 8))
    plt.plot(thresholds, accuracies, 'o-', label='Accuracy', linewidth=2)
    plt.plot(thresholds, precisions, 's-', label='Precision', linewidth=2)
    plt.plot(thresholds, recalls, '^-', label='Recall', linewidth=2)
    plt.plot(thresholds, f1_scores, 'd-', label='F1-Score', linewidth=2)
    plt.plot(thresholds, specificities, 'v-', label='Specificity', linewidth=2)
    
    plt.xlabel('Threshold')
    plt.ylabel('Score')
    plt.title('Performance Metrics vs Threshold - Pangolin on SpliceVarDB')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.xscale('log')
    
    plt.savefig(os.path.join(output_dir, 'threshold_performance.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'threshold_performance.pdf'), bbox_inches='tight')
    plt.close()
    
    logging.info("✅ Threshold performance plot saved")

def generate_classification_report(y_true, y_scores, threshold, output_dir):
    """Generate detailed classification report"""
    
    y_pred = (y_scores >= threshold).astype(int)
    
    # Generate classification report
    report = classification_report(y_true, y_pred, 
                                 target_names=['Normal', 'Splice-altering'],
                                 digits=4)
    
    # Save to file
    report_file = os.path.join(output_dir, f'classification_report_threshold_{threshold}.txt')
    with open(report_file, 'w') as f:
        f.write(f"Classification Report - Pangolin on SpliceVarDB\n")
        f.write(f"Threshold: {threshold}\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("=" * 60 + "\n\n")
        f.write(report)
    
    logging.info(f"✅ Classification report saved (threshold={threshold})")

def main():
    parser = argparse.ArgumentParser(
        description="Generate detailed visualizations for Pangolin performance on SpliceVarDB"
    )
    
    parser.add_argument(
        '--labels',
        type=str,
        default='/Users/byronsun/Desktop/AS_复现模型/BIB_review/data/processed_data/splicevardb/evo2_splicevardb_dataset_dedup.tsv',
        help='SpliceVarDB labels TSV file'
    )
    
    parser.add_argument(
        '--pangolin_vcf',
        type=str,
        default='results/pangolin_splicevardb_full.vcf',
        help='Pangolin output VCF file'
    )
    
    parser.add_argument(
        '--output_dir',
        type=str,
        default='results/detailed_analysis',
        help='Output directory for visualizations'
    )
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    logging.info("=== Pangolin Detailed Visualization Analysis ===")
    
    # Load data
    labels_df = load_splicevardb_labels(args.labels)
    pangolin_df = load_pangolin_vcf(args.pangolin_vcf)
    
    # Merge data
    logging.info("Merging data...")
    merged = labels_df.join(pangolin_df, how='inner')
    
    logging.info(f"Matched variants: {len(merged)}")
    logging.info(f"Positive labels: {merged['binary_label'].sum()}")
    logging.info(f"Negative labels: {len(merged) - merged['binary_label'].sum()}")
    logging.info(f"Non-zero scores: {(merged['pangolin_score'] > 0).sum()}")
    
    # Extract arrays
    y_true = merged['binary_label'].values
    y_scores = merged['pangolin_score'].values
    
    # Compute detailed metrics
    logging.info("Computing detailed metrics...")
    thresholds = [0.001, 0.01, 0.05, 0.1, 0.2]
    metrics_data = compute_detailed_metrics(y_true, y_scores, thresholds)
    
    # Generate visualizations
    logging.info("Generating visualizations...")
    
    # 1. Confusion matrices
    plot_confusion_matrices(y_true, y_scores, thresholds, args.output_dir)
    
    # 2. ROC curve
    plot_roc_curve(y_true, y_scores, args.output_dir)
    
    # 3. Precision-Recall curve
    plot_precision_recall_curve(y_true, y_scores, args.output_dir)
    
    # 4. Score distribution
    plot_score_distribution(y_true, y_scores, args.output_dir)
    
    # 5. Threshold performance
    plot_threshold_performance(metrics_data, args.output_dir)
    
    # 6. Classification reports
    for threshold in [0.01, 0.05, 0.1]:
        generate_classification_report(y_true, y_scores, threshold, args.output_dir)
    
    # Save comprehensive results
    comprehensive_results = {
        'model': 'Pangolin',
        'dataset': 'SpliceVarDB_Human_Gold_Standard',
        'analysis_timestamp': datetime.now().strftime("%Y%m%d_%H%M%S"),
        'total_variants': len(merged),
        'positive_variants': int(merged['binary_label'].sum()),
        'negative_variants': int(len(merged) - merged['binary_label'].sum()),
        'non_zero_scores': int((merged['pangolin_score'] > 0).sum()),
        'zero_score_percentage': float((merged['pangolin_score'] == 0).sum() / len(merged) * 100),
        'overall_metrics': {
            'auroc': metrics_data['auroc'],
            'auprc': metrics_data['auprc']
        },
        'threshold_metrics': metrics_data['metrics_by_threshold'],
        'score_statistics': {
            'min': float(y_scores.min()),
            'max': float(y_scores.max()),
            'mean': float(y_scores.mean()),
            'std': float(y_scores.std()),
            'median': float(np.median(y_scores)),
            'q25': float(np.percentile(y_scores, 25)),
            'q75': float(np.percentile(y_scores, 75))
        }
    }
    
    # Save comprehensive results
    results_file = os.path.join(args.output_dir, 'comprehensive_analysis.json')
    with open(results_file, 'w') as f:
        json.dump(comprehensive_results, f, indent=2)
    
    logging.info(f"✅ Comprehensive results saved to: {results_file}")
    
    # Print summary
    print(f"\n{'='*80}")
    print(f"PANGOLIN DETAILED ANALYSIS SUMMARY")
    print(f"{'='*80}")
    print(f"Dataset: {len(merged):,} variants")
    print(f"Positive: {comprehensive_results['positive_variants']:,} ({comprehensive_results['positive_variants']/len(merged)*100:.1f}%)")
    print(f"Negative: {comprehensive_results['negative_variants']:,} ({comprehensive_results['negative_variants']/len(merged)*100:.1f}%)")
    print(f"Non-zero scores: {comprehensive_results['non_zero_scores']:,} ({100-comprehensive_results['zero_score_percentage']:.1f}%)")
    print(f"\nOverall Performance:")
    print(f"  AUROC: {metrics_data['auroc']:.4f}")
    print(f"  AUPRC: {metrics_data['auprc']:.4f}")
    print(f"\nRecommended Threshold (0.01):")
    rec_metrics = metrics_data['metrics_by_threshold']['threshold_0.01']
    print(f"  Accuracy: {rec_metrics['accuracy']:.4f}")
    print(f"  Precision: {rec_metrics['precision']:.4f}")
    print(f"  Recall: {rec_metrics['recall']:.4f}")
    print(f"  F1-Score: {rec_metrics['f1_score']:.4f}")
    print(f"  Specificity: {rec_metrics['specificity']:.4f}")
    print(f"\nAll visualizations saved to: {args.output_dir}/")
    print(f"{'='*80}")

if __name__ == '__main__':
    main()
