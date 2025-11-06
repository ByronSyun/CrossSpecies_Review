#!/usr/bin/env python3
"""
Evaluate Pangolin performance on SpliceVarDB dataset.
"""

import pandas as pd
import numpy as np
import json
import argparse
import re
import logging
from sklearn.metrics import roc_auc_score, average_precision_score, accuracy_score
from sklearn.metrics import precision_score, recall_score, f1_score, confusion_matrix
from datetime import datetime
import os

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
    # Remove 'chr' prefix and use / separator to match VCF format
    chrom_no_prefix = df['chrom'].str.replace('chr', '', regex=False)
    df['vcf_variant_id'] = chrom_no_prefix + ':' + (df['pos_0based'] + 1).astype(str) + '_' + df['ref'] + '/' + df['alt']
    
    # Convert label to binary
    df['binary_label'] = (df['label'] == 'splice-altering').astype(int)
    
    logging.info(f"Loaded {len(df)} variants")
    logging.info(f"Splice-altering: {df['binary_label'].sum()}")
    logging.info(f"Normal: {len(df) - df['binary_label'].sum()}")
    logging.info(f"Sample variant IDs: {df['vcf_variant_id'].head(3).tolist()}")
    
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
    logging.info(f"Sample variant IDs: {df.index[:3].tolist()}")
    
    return df

def compute_metrics(y_true, y_scores, threshold=0.01):
    """Compute classification metrics"""
    y_pred = (y_scores >= threshold).astype(int)
    
    # Basic metrics
    auroc = roc_auc_score(y_true, y_scores)
    auprc = average_precision_score(y_true, y_scores)
    accuracy = accuracy_score(y_true, y_pred)
    precision = precision_score(y_true, y_pred, zero_division=0)
    recall = recall_score(y_true, y_pred, zero_division=0)
    f1 = f1_score(y_true, y_pred, zero_division=0)
    
    # Confusion matrix
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
    
    return {
        'threshold': threshold,
        'auroc': auroc,
        'auprc': auprc,
        'accuracy': accuracy,
        'precision': precision,
        'recall': recall,
        'f1_score': f1,
        'specificity': specificity,
        'tp': int(tp), 'fp': int(fp), 'tn': int(tn), 'fn': int(fn)
    }

def main():
    parser = argparse.ArgumentParser(description="Evaluate Pangolin on SpliceVarDB")
    parser.add_argument('--labels', default='/mnt/userdata4/splicing/Evo2Splicevardb/splicevardb_data/evo2_splicevardb_dataset_dedup.tsv',
                       help='SpliceVarDB labels TSV file')
    parser.add_argument('--pangolin_vcf', default='results/pangolin_splicevardb_full.vcf',
                       help='Pangolin output VCF file')
    parser.add_argument('--output_dir', default='results/',
                       help='Output directory')
    
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    
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
    
    # Compute metrics
    y_true = merged['binary_label'].values
    y_scores = merged['pangolin_score'].values
    
    # Multiple thresholds
    thresholds = [0.001, 0.01, 0.05, 0.1, 0.2]
    all_metrics = {}
    
    for thresh in thresholds:
        all_metrics[f'threshold_{thresh}'] = compute_metrics(y_true, y_scores, threshold=thresh)
    
    # Results summary
    results = {
        'model': 'Pangolin',
        'dataset': 'SpliceVarDB_Human_Gold_Standard',
        'evaluation_timestamp': datetime.now().strftime("%Y%m%d_%H%M%S"),
        'total_variants': len(merged),
        'positive_variants': int(merged['binary_label'].sum()),
        'negative_variants': int(len(merged) - merged['binary_label'].sum()),
        'non_zero_scores': int((merged['pangolin_score'] > 0).sum()),
        'metrics_by_threshold': all_metrics,
        'score_statistics': {
            'min': float(y_scores.min()),
            'max': float(y_scores.max()),
            'mean': float(y_scores.mean()),
            'std': float(y_scores.std()),
            'zero_count': int((y_scores == 0).sum()),
            'zero_percentage': float((y_scores == 0).sum() / len(y_scores) * 100)
        }
    }
    
    # Save results
    output_file = os.path.join(args.output_dir, 'pangolin_splicevardb_evaluation.json')
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    # Print summary
    print(f"\n{'='*60}")
    print(f"PANGOLIN PERFORMANCE ON SPLICEVARDB (HUMAN GOLD STANDARD)")
    print(f"{'='*60}")
    print(f"Dataset: {len(merged):,} variants")
    print(f"Positive: {results['positive_variants']:,} ({results['positive_variants']/len(merged)*100:.1f}%)")
    print(f"Negative: {results['negative_variants']:,} ({results['negative_variants']/len(merged)*100:.1f}%)")
    print(f"Non-zero scores: {results['non_zero_scores']:,} ({results['non_zero_scores']/len(merged)*100:.1f}%)")
    print(f"Zero scores: {results['score_statistics']['zero_count']:,} ({results['score_statistics']['zero_percentage']:.1f}%)")
    
    print(f"\nSCORE STATISTICS:")
    print(f"  Min: {results['score_statistics']['min']:.4f}")
    print(f"  Max: {results['score_statistics']['max']:.4f}")
    print(f"  Mean: {results['score_statistics']['mean']:.4f}")
    print(f"  Std: {results['score_statistics']['std']:.4f}")
    
    print(f"\nPERFORMANCE BY THRESHOLD:")
    print(f"{'Threshold':<10} {'AUROC':<8} {'AUPRC':<8} {'Acc':<8} {'Prec':<8} {'Rec':<8} {'F1':<8} {'Spec':<8}")
    print("-" * 70)
    
    for thresh_name, metrics in all_metrics.items():
        thresh_val = metrics['threshold']
        print(f"{thresh_val:<10.3f} {metrics['auroc']:<8.4f} {metrics['auprc']:<8.4f} "
              f"{metrics['accuracy']:<8.4f} {metrics['precision']:<8.4f} {metrics['recall']:<8.4f} "
              f"{metrics['f1_score']:<8.4f} {metrics['specificity']:<8.4f}")
    
    # Recommended metrics (threshold=0.01)
    rec_metrics = all_metrics['threshold_0.01']
    print(f"\n{'='*40}")
    print(f"RECOMMENDED PERFORMANCE (threshold=0.01)")
    print(f"{'='*40}")
    print(f"AUROC: {rec_metrics['auroc']:.4f}")
    print(f"AUPRC: {rec_metrics['auprc']:.4f}")
    print(f"Accuracy: {rec_metrics['accuracy']:.4f}")
    print(f"Precision: {rec_metrics['precision']:.4f}")
    print(f"Recall: {rec_metrics['recall']:.4f}")
    print(f"F1-score: {rec_metrics['f1_score']:.4f}")
    print(f"Specificity: {rec_metrics['specificity']:.4f}")
    
    print(f"\nResults saved to: {output_file}")
    logging.info(f"✅ Evaluation completed successfully!")

if __name__ == '__main__':
    main()
