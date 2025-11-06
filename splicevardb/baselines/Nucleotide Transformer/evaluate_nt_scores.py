#!/usr/bin/env python3
import os
import json
import argparse
import logging
import numpy as np
import pandas as pd
from sklearn.metrics import (
    roc_auc_score, average_precision_score, accuracy_score,
    precision_score, recall_score, f1_score, confusion_matrix
)
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime

logging.basicConfig(level=logging.INFO, format='%(asctime)s - [%(levelname)s] - %(message)s')

def compute_metrics(y_true, y_scores, threshold=0.01):
    """Compute classification metrics at a given threshold"""
    y_pred = (y_scores >= threshold).astype(int)
    
    # Basic metrics
    auroc = roc_auc_score(y_true, y_scores)
    auprc = average_precision_score(y_true, y_scores)
    acc = accuracy_score(y_true, y_pred)
    prec = precision_score(y_true, y_pred, zero_division=0)
    rec = recall_score(y_true, y_pred, zero_division=0)
    f1 = f1_score(y_true, y_pred, zero_division=0)
    
    # Confusion matrix
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    spec = tn / (tn + fp) if (tn + fp) > 0 else 0
    
    return {
        'threshold': threshold,
        'auroc': float(auroc),
        'auprc': float(auprc),
        'accuracy': float(acc),
        'precision': float(prec),
        'recall': float(rec),
        'f1_score': float(f1),
        'specificity': float(spec),
        'tp': int(tp), 'fp': int(fp), 'tn': int(tn), 'fn': int(fn)
    }

def plot_score_distribution(y_true, y_scores, output_dir):
    """Plot score distribution by class"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Histogram
    pos_scores = y_scores[y_true == 1]
    neg_scores = y_scores[y_true == 0]
    
    ax1.hist(neg_scores, bins=50, alpha=0.7, label='Normal', density=True, color='blue')
    ax1.hist(pos_scores, bins=50, alpha=0.7, label='Splice-altering', density=True, color='red')
    ax1.set_xlabel('NT Score (cosine distance)')
    ax1.set_ylabel('Density')
    ax1.set_title('Score Distribution by Class')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Box plot
    data_to_plot = [neg_scores, pos_scores]
    ax2.boxplot(data_to_plot, labels=['Normal', 'Splice-altering'])
    ax2.set_ylabel('NT Score (cosine distance)')
    ax2.set_title('Score Distribution Box Plot')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'nt_score_distribution.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    logging.info("✅ Score distribution plot saved")

def main():
    parser = argparse.ArgumentParser(
        description="Evaluate Nucleotide Transformer variant scores"
    )
    parser.add_argument('--labels', required=True, 
                       help='TSV file with variant labels')
    parser.add_argument('--scores', required=True, 
                       help='CSV file with variant scores from NT')
    parser.add_argument('--label_column', default='label',
                       help='Label column name (default: label)')
    parser.add_argument('--out_dir', required=True,
                       help='Output directory for results')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.out_dir, exist_ok=True)
    
    logging.info("=== Nucleotide Transformer Evaluation ===")
    logging.info(f"Labels: {args.labels}")
    logging.info(f"Scores: {args.scores}")
    
    # Load labels
    logging.info("Loading labels...")
    df_lab = pd.read_csv(args.labels, sep='\t')
    
    # Convert labels to binary
    if 'binary_label' in df_lab.columns:
        df_lab['y'] = df_lab['binary_label']
    else:
        df_lab['y'] = (df_lab[args.label_column] == 'splice-altering').astype(int)
    
    # Create variant ID for matching - use original variant_id if available
    if 'variant_id' in df_lab.columns:
        df_lab['vid'] = df_lab['variant_id']
    elif {'chrom', 'pos_0based', 'ref', 'alt'}.issubset(df_lab.columns):
        chrom = df_lab['chrom'].astype(str).str.replace('^chr', '', regex=True)
        df_lab['vid'] = chrom + ':' + (df_lab['pos_0based'] + 1).astype(str) + '_' + df_lab['ref'] + '/' + df_lab['alt']
    else:
        raise ValueError("No suitable columns found for creating variant IDs")
    
    df_lab = df_lab[['vid', 'y']].drop_duplicates()
    
    # Load scores
    logging.info("Loading scores...")
    df_sc = pd.read_csv(args.scores)
    df_sc['vid'] = df_sc['variant_id']
    
    # Use cosine distance as primary score
    score_col = 'score_cosine' if 'score_cosine' in df_sc.columns else 'score_l2'
    logging.info(f"Using score column: {score_col}")
    
    # Merge data
    merged = df_lab.merge(df_sc[['vid', score_col]], on='vid', how='inner')
    y_true = merged['y'].values
    y_scores = merged[score_col].values
    
    logging.info(f"Merged: {len(merged)} variants")
    logging.info(f"Positive labels: {merged['y'].sum()} ({merged['y'].sum()/len(merged)*100:.1f}%)")
    logging.info(f"Non-zero scores: {(y_scores > 0).sum()} ({(y_scores > 0).sum()/len(y_scores)*100:.1f}%)")
    
    # Compute metrics
    logging.info("Computing metrics...")
    metrics = compute_metrics(y_true, y_scores, threshold=0.01)
    
    # Save results
    results = {
        'model': 'Nucleotide_Transformer_v2_500M',
        'dataset': 'SpliceVarDB_Human_Gold_Standard',
        'evaluation_timestamp': datetime.now().strftime("%Y%m%d_%H%M%S"),
        'total_variants': len(merged),
        'positive_variants': int(merged['y'].sum()),
        'negative_variants': int(len(merged) - merged['y'].sum()),
        'non_zero_scores': int((y_scores > 0).sum()),
        'zero_score_percentage': float((y_scores == 0).sum() / len(y_scores) * 100),
        'metrics': metrics,
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
    
    # Save JSON results
    results_file = os.path.join(args.out_dir, 'nt_evaluation_results.json')
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    # Generate plots
    plot_score_distribution(y_true, y_scores, args.out_dir)
    
    # Print summary
    print(f"\n{'='*60}")
    print("NUCLEOTIDE TRANSFORMER EVALUATION RESULTS")
    print(f"{'='*60}")
    print(f"Dataset: {len(merged):,} variants")
    print(f"Positive: {results['positive_variants']:,} ({results['positive_variants']/len(merged)*100:.1f}%)")
    print(f"Negative: {results['negative_variants']:,} ({results['negative_variants']/len(merged)*100:.1f}%)")
    print(f"Non-zero scores: {results['non_zero_scores']:,} ({100-results['zero_score_percentage']:.1f}%)")
    print(f"\nPerformance (threshold=0.01):")
    print(f"  AUROC: {metrics['auroc']:.4f}")
    print(f"  AUPRC: {metrics['auprc']:.4f}")
    print(f"  Accuracy: {metrics['accuracy']:.4f}")
    print(f"  Precision: {metrics['precision']:.4f}")
    print(f"  Recall: {metrics['recall']:.4f}")
    print(f"  F1-Score: {metrics['f1_score']:.4f}")
    print(f"  Specificity: {metrics['specificity']:.4f}")
    print(f"\nScore Statistics:")
    print(f"  Range: [{results['score_statistics']['min']:.4f}, {results['score_statistics']['max']:.4f}]")
    print(f"  Mean ± Std: {results['score_statistics']['mean']:.4f} ± {results['score_statistics']['std']:.4f}")
    print(f"  Median: {results['score_statistics']['median']:.4f}")
    print(f"\nResults saved to: {results_file}")
    print(f"{'='*60}")
    
    logging.info(f"✅ Evaluation completed! Results saved to: {results_file}")

if __name__ == '__main__':
    main()