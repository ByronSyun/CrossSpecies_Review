#!/usr/bin/env python3
"""
Evaluate SpliceTransformer Cross-Species Performance on RatGTEx

This script evaluates the performance of SpliceTransformer on the RatGTEx 
silver-standard dataset to quantify cross-species prediction capabilities.

Key Features:
- Computes comprehensive binary classification metrics (AUROC, AUPRC, etc.)
- Generates ROC and Precision-Recall curves
- Calculates optimal thresholds for different objectives
- Creates publication-ready performance plots
- Saves results in formats compatible with cross-model comparison

Author: Generated for AS 复现模型 project
Date: 2024
"""

import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import (
    roc_auc_score, average_precision_score, roc_curve, precision_recall_curve,
    confusion_matrix, accuracy_score, precision_score, recall_score, f1_score
)
import argparse
import logging
import os
from pathlib import Path

# Setup logging and plotting
logging.basicConfig(level=logging.INFO, format='%(asctime)s - [%(levelname)s] - %(message)s')
plt.style.use('default')
sns.set_palette("husl")

def load_splicetransformer_csv_results(predictions_csv, labels_csv):
    """
    Load SpliceTransformer results from CSV predictions and CSV labels.
    
    Args:
        predictions_csv (str): Path to SpliceTransformer CSV output
        labels_csv (str): Path to CSV file with labels
        
    Returns:
        tuple: (scores, labels, metadata) or (None, None, None) if failed
    """
    try:
        # Load labels
        logging.info(f"Loading labels from: {labels_csv}")
        df_labels = pd.read_csv(labels_csv)
        logging.info(f"Loaded {len(df_labels):,} labeled variants")
        
        # Load SpliceTransformer CSV predictions
        logging.info(f"Loading SpliceTransformer predictions from: {predictions_csv}")
        
        # Read CSV file
        df_predictions = pd.read_csv(predictions_csv)
        logging.info(f"Loaded {len(df_predictions):,} predictions from CSV")
        
        if len(df_predictions) == 0:
            logging.error("No predictions found in CSV file")
            return None, None, None
        
        # Create variant_id for matching
        df_predictions['variant_id'] = (df_predictions['#CHROM'].astype(str) + ':' + 
                                       df_predictions['POS'].astype(str) + '_' + 
                                       df_predictions['REF'] + '/' + 
                                       df_predictions['ALT'])
        
        # Merge predictions with labels using variant_id
        merged = df_predictions.merge(df_labels[['variant_id', 'label']], on='variant_id', how='inner')
        
        if len(merged) == 0:
            logging.error("No matching variants found between predictions and labels")
            return None, None, None
        
        logging.info(f"Successfully matched {len(merged):,} variant records")
        
        # Handle duplicates: same variant may appear multiple times due to multiple tissues
        # Take the first occurrence of each variant_id (they should have the same label)
        merged_dedup = merged.drop_duplicates(subset=['variant_id'], keep='first')
        logging.info(f"After deduplication: {len(merged_dedup):,} unique variants")
        
        # Verify labels are consistent for duplicated variants
        label_check = merged.groupby('variant_id')['label'].nunique()
        inconsistent_labels = (label_check > 1).sum()
        if inconsistent_labels > 0:
            logging.warning(f"Found {inconsistent_labels} variants with inconsistent labels across tissues")
        
        # Extract final scores and labels
        scores = merged_dedup['score'].values
        labels = merged_dedup['label'].values
        
        # Calculate metadata
        metadata = {
            'total_variants': len(scores),
            'positive_samples': np.sum(labels),
            'negative_samples': len(labels) - np.sum(labels),
            'score_range': (np.min(scores), np.max(scores)),
            'score_stats': {
                'mean': np.mean(scores),
                'std': np.std(scores),
                'median': np.median(scores)
            }
        }
        
        logging.info(f"Final dataset: {len(scores):,} predictions")
        logging.info(f"Positive samples: {metadata['positive_samples']:,} ({metadata['positive_samples']/len(labels)*100:.1f}%)")
        logging.info(f"Score range: [{metadata['score_range'][0]:.3f}, {metadata['score_range'][1]:.3f}]")
        
        return scores, labels, metadata
        
    except Exception as e:
        logging.error(f"Failed to load SpliceTransformer results: {e}")
        import traceback
        logging.error(traceback.format_exc())
        return None, None, None

def load_splicetransformer_vcf_results(predictions_vcf, labels_csv):
    """
    Load SpliceTransformer results from VCF predictions and CSV labels.
    
    Args:
        predictions_vcf (str): Path to SpliceTransformer VCF output
        labels_csv (str): Path to CSV file with labels
        
    Returns:
        tuple: (scores, labels, metadata) or (None, None, None) if failed
    """
    try:
        # Load labels
        logging.info(f"Loading labels from: {labels_csv}")
        df_labels = pd.read_csv(labels_csv)
        logging.info(f"Loaded {len(df_labels):,} labeled variants")
        
        # Load SpliceTransformer VCF predictions
        logging.info(f"Loading SpliceTransformer predictions from: {predictions_vcf}")
        
        # Parse VCF file
        vcf_data = []
        with open(predictions_vcf, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue  # Skip header lines
                
                fields = line.strip().split('\t')
                if len(fields) >= 8:  # Standard VCF format
                    chrom, pos, variant_id, ref, alt, qual, filter_field, info = fields[:8]
                    
                    # Extract SpliceTransformer score from INFO field
                    score = None
                    if 'SCORE=' in info:
                        # Extract score value
                        for item in info.split(';'):
                            if item.startswith('SCORE='):
                                score = float(item.split('=')[1])
                                break
                    elif info != '.' and info != '':
                        # Sometimes the score is directly in the INFO field
                        try:
                            score = float(info)
                        except:
                            pass
                    
                    if score is not None:
                        vcf_data.append({
                            'CHROM': chrom,
                            'POS': int(pos),
                            'variant_id': variant_id,
                            'REF': ref,
                            'ALT': alt,
                            'score': score
                        })
        
        df_predictions = pd.DataFrame(vcf_data)
        logging.info(f"Parsed {len(df_predictions):,} predictions from VCF")
        
        if len(df_predictions) == 0:
            logging.error("No predictions found in VCF file")
            return None, None, None
        
        # Merge predictions with labels using variant_id
        merged = df_predictions.merge(df_labels[['variant_id', 'label']], on='variant_id', how='inner')
        
        if len(merged) == 0:
            logging.error("No matching variants found between predictions and labels")
            return None, None, None
        
        logging.info(f"Successfully matched {len(merged):,} variants")
        
        # Extract final scores and labels
        scores = merged['score'].values
        labels = merged['label'].values
        
        # Calculate metadata
        metadata = {
            'total_variants': len(scores),
            'positive_samples': np.sum(labels),
            'negative_samples': len(labels) - np.sum(labels),
            'score_range': (np.min(scores), np.max(scores)),
            'score_stats': {
                'mean': np.mean(scores),
                'std': np.std(scores),
                'median': np.median(scores)
            }
        }
        
        logging.info(f"Final dataset: {len(scores):,} predictions")
        logging.info(f"Positive samples: {metadata['positive_samples']:,} ({metadata['positive_samples']/len(labels)*100:.1f}%)")
        logging.info(f"Score range: [{metadata['score_range'][0]:.3f}, {metadata['score_range'][1]:.3f}]")
        
        return scores, labels, metadata
        
    except Exception as e:
        logging.error(f"Failed to load SpliceTransformer results: {e}")
        import traceback
        logging.error(traceback.format_exc())
        return None, None, None

def load_splicetransformer_results(json_file):
    """
    Load SpliceTransformer results from JSON file.
    
    Args:
        json_file (str): Path to JSON results file
        
    Returns:
        tuple: (scores, labels, metadata) or (None, None, None) if failed
    """
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
        
        results = data['results']
        scores = [r['splicetransformer_score'] for r in results]
        labels = [r['label'] for r in results]
        metadata = data['metadata']
        
        logging.info(f"Loaded {len(scores):,} predictions from {json_file}")
        logging.info(f"Score range: {min(scores):.3f} to {max(scores):.3f}")
        logging.info(f"Class distribution: {sum(labels)} positive, {len(labels) - sum(labels)} negative")
        
        return np.array(scores), np.array(labels), metadata
        
    except Exception as e:
        logging.error(f"Failed to load results from {json_file}: {e}")
        return None, None, None

def compute_binary_classification_metrics(scores, labels):
    """
    Compute comprehensive binary classification metrics.
    
    Args:
        scores (np.array): Prediction scores
        labels (np.array): True binary labels
        
    Returns:
        dict: Dictionary containing all metrics
    """
    
    # Basic metrics
    auroc = roc_auc_score(labels, scores)
    auprc = average_precision_score(labels, scores)
    
    # ROC curve data
    fpr, tpr, roc_thresholds = roc_curve(labels, scores)
    
    # Precision-Recall curve data  
    precision, recall, pr_thresholds = precision_recall_curve(labels, scores)
    
    # Find optimal thresholds
    # 1. Youden's J statistic (sensitivity + specificity - 1)
    j_scores = tpr - fpr
    optimal_idx = np.argmax(j_scores)
    optimal_threshold_youden = roc_thresholds[optimal_idx]
    
    # 2. F1-optimal threshold
    f1_scores = 2 * precision * recall / (precision + recall + 1e-10)
    f1_optimal_idx = np.argmax(f1_scores)
    optimal_threshold_f1 = pr_thresholds[f1_optimal_idx] if f1_optimal_idx < len(pr_thresholds) else pr_thresholds[-1]
    
    # Performance at optimal thresholds
    pred_youden = (scores >= optimal_threshold_youden).astype(int)
    pred_f1 = (scores >= optimal_threshold_f1).astype(int)
    
    # Confusion matrices
    tn_y, fp_y, fn_y, tp_y = confusion_matrix(labels, pred_youden).ravel()
    tn_f1, fp_f1, fn_f1, tp_f1 = confusion_matrix(labels, pred_f1).ravel()
    
    # Compile metrics
    metrics = {
        'auroc': auroc,
        'auprc': auprc,
        'optimal_thresholds': {
            'youden': optimal_threshold_youden,
            'f1': optimal_threshold_f1
        },
        'youden_metrics': {
            'threshold': optimal_threshold_youden,
            'accuracy': accuracy_score(labels, pred_youden),
            'precision': precision_score(labels, pred_youden, zero_division=0),
            'recall': recall_score(labels, pred_youden, zero_division=0),
            'f1': f1_score(labels, pred_youden, zero_division=0),
            'specificity': tn_y / (tn_y + fp_y) if (tn_y + fp_y) > 0 else 0,
            'confusion_matrix': [[tn_y, fp_y], [fn_y, tp_y]]
        },
        'f1_optimal_metrics': {
            'threshold': optimal_threshold_f1,
            'accuracy': accuracy_score(labels, pred_f1),
            'precision': precision_score(labels, pred_f1, zero_division=0),
            'recall': recall_score(labels, pred_f1, zero_division=0),
            'f1': f1_score(labels, pred_f1, zero_division=0),
            'specificity': tn_f1 / (tn_f1 + fp_f1) if (tn_f1 + fp_f1) > 0 else 0,
            'confusion_matrix': [[tn_f1, fp_f1], [fn_f1, tp_f1]]
        },
        'curve_data': {
            'roc': {'fpr': fpr.tolist(), 'tpr': tpr.tolist(), 'thresholds': roc_thresholds.tolist()},
            'pr': {'precision': precision.tolist(), 'recall': recall.tolist(), 'thresholds': pr_thresholds.tolist()}
        }
    }
    
    return metrics

def create_performance_plots(scores, labels, metrics, output_dir, model_name="SpliceTransformer"):
    """
    Create ROC and Precision-Recall curves.
    
    Args:
        scores (np.array): Prediction scores
        labels (np.array): True binary labels  
        metrics (dict): Computed metrics
        output_dir (str): Directory to save plots
        model_name (str): Model name for plot titles
    """
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # ROC Curve
    roc_data = metrics['curve_data']['roc']
    ax1.plot(roc_data['fpr'], roc_data['tpr'], linewidth=2, 
            label=f'{model_name} (AUROC = {metrics["auroc"]:.3f})')
    ax1.plot([0, 1], [0, 1], 'k--', alpha=0.5, label='Random')
    ax1.set_xlabel('False Positive Rate')
    ax1.set_ylabel('True Positive Rate')
    ax1.set_title(f'{model_name} ROC Curve\nRatGTEx Silver Standard')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Precision-Recall Curve
    pr_data = metrics['curve_data']['pr']
    baseline_precision = np.sum(labels) / len(labels)
    ax2.plot(pr_data['recall'], pr_data['precision'], linewidth=2,
            label=f'{model_name} (AUPRC = {metrics["auprc"]:.3f})')
    ax2.axhline(y=baseline_precision, color='k', linestyle='--', alpha=0.5, 
               label=f'Random (Precision = {baseline_precision:.3f})')
    ax2.set_xlabel('Recall')
    ax2.set_ylabel('Precision') 
    ax2.set_title(f'{model_name} Precision-Recall Curve\nRatGTEx Silver Standard')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save plot
    plot_path = os.path.join(output_dir, f'{model_name.lower()}_performance_curves.png')
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logging.info(f"Performance plots saved to: {plot_path}")

def generate_performance_summary(metrics, metadata, output_file):
    """
    Generate a summary report of model performance.
    
    Args:
        metrics (dict): Computed metrics
        metadata (dict): Model metadata
        output_file (str): Path to save summary
    """
    
    summary = {
        'model_info': {
            'model_name': 'SpliceTransformer',
            'dataset': 'RatGTEx Silver Standard',
            'total_variants': metadata.get('total_processed', 'Unknown'),
            'timestamp': metadata.get('timestamp', 'Unknown')
        },
        'performance_metrics': {
            'auroc': round(metrics['auroc'], 4),
            'auprc': round(metrics['auprc'], 4)
        },
        'optimal_performance': {
            'youden_threshold': {
                'threshold': round(metrics['optimal_thresholds']['youden'], 4),
                'accuracy': round(metrics['youden_metrics']['accuracy'], 4),
                'precision': round(metrics['youden_metrics']['precision'], 4),
                'recall': round(metrics['youden_metrics']['recall'], 4),
                'f1_score': round(metrics['youden_metrics']['f1'], 4),
                'specificity': round(metrics['youden_metrics']['specificity'], 4)
            },
            'f1_optimal_threshold': {
                'threshold': round(metrics['optimal_thresholds']['f1'], 4),
                'accuracy': round(metrics['f1_optimal_metrics']['accuracy'], 4),
                'precision': round(metrics['f1_optimal_metrics']['precision'], 4),
                'recall': round(metrics['f1_optimal_metrics']['recall'], 4),
                'f1_score': round(metrics['f1_optimal_metrics']['f1'], 4),
                'specificity': round(metrics['f1_optimal_metrics']['specificity'], 4)
            }
        },
        'confusion_matrices': {
            'youden_threshold': metrics['youden_metrics']['confusion_matrix'].tolist() if hasattr(metrics['youden_metrics']['confusion_matrix'], 'tolist') else metrics['youden_metrics']['confusion_matrix'],
            'f1_optimal_threshold': metrics['f1_optimal_metrics']['confusion_matrix'].tolist() if hasattr(metrics['f1_optimal_metrics']['confusion_matrix'], 'tolist') else metrics['f1_optimal_metrics']['confusion_matrix']
        }
    }
    
    # Save summary
    with open(output_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    # Print summary to console
    logging.info("\n" + "="*60)
    logging.info("SPLICETRANSFORMER CROSS-SPECIES PERFORMANCE SUMMARY")
    logging.info("="*60)
    logging.info(f"Dataset: {summary['model_info']['dataset']}")
    logging.info(f"Total Variants: {summary['model_info']['total_variants']:,}")
    logging.info("")
    logging.info("MAIN PERFORMANCE METRICS:")
    logging.info(f"  AUROC: {summary['performance_metrics']['auroc']:.4f}")
    logging.info(f"  AUPRC: {summary['performance_metrics']['auprc']:.4f}")
    logging.info("")
    logging.info("OPTIMAL THRESHOLD PERFORMANCE:")
    logging.info("  Youden's J Statistic:")
    youden = summary['optimal_performance']['youden_threshold']
    logging.info(f"    Threshold: {youden['threshold']:.4f}")
    logging.info(f"    Accuracy: {youden['accuracy']:.4f}")
    logging.info(f"    Precision: {youden['precision']:.4f}")
    logging.info(f"    Recall: {youden['recall']:.4f}")
    logging.info(f"    F1-Score: {youden['f1_score']:.4f}")
    logging.info(f"    Specificity: {youden['specificity']:.4f}")
    logging.info("")
    logging.info("  F1-Optimal Threshold:")
    f1_opt = summary['optimal_performance']['f1_optimal_threshold']
    logging.info(f"    Threshold: {f1_opt['threshold']:.4f}")
    logging.info(f"    Accuracy: {f1_opt['accuracy']:.4f}")
    logging.info(f"    Precision: {f1_opt['precision']:.4f}")
    logging.info(f"    Recall: {f1_opt['recall']:.4f}")
    logging.info(f"    F1-Score: {f1_opt['f1_score']:.4f}")
    logging.info(f"    Specificity: {f1_opt['specificity']:.4f}")
    logging.info("="*60)
    
    return summary

def main():
    parser = argparse.ArgumentParser(
        description="Evaluate SpliceTransformer cross-species performance on RatGTEx"
    )
    
    # Default paths
    default_predictions = "/Users/byronsun/Desktop/AS_复现模型/BIB_review/code/ratGTEx/baselines/spliceTransformer(need_hg38)/SpliceTransformer_ratgtex_predictions.vcf"
    default_labels = "/Users/byronsun/Desktop/AS_复现模型/BIB_review/code/ratGTEx/baselines/spliceTransformer(need_hg38)/ratgtex_labels.csv"
    default_output_dir = "/Users/byronsun/Desktop/AS_复现模型/BIB_review/code/ratGTEx/baselines/spliceTransformer(need_hg38)"
    
    parser.add_argument(
        '--predictions',
        type=str,
        default=default_predictions,
        help=f'Path to SpliceTransformer CSV/VCF predictions file (default: {default_predictions})'
    )
    
    parser.add_argument(
        '--labels',
        type=str,
        default=default_labels,
        help=f'Path to CSV file with labels (default: {default_labels})'
    )
    
    parser.add_argument(
        '--output_dir',
        type=str,
        default=default_output_dir,
        help=f'Directory to save evaluation outputs (default: {default_output_dir})'
    )
    
    args = parser.parse_args()
    
    logging.info("=== SpliceTransformer Cross-Species Performance Evaluation ===")
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load results (try CSV first, then VCF)
    if args.predictions.endswith('.csv'):
        scores, labels, metadata = load_splicetransformer_csv_results(args.predictions, args.labels)
    else:
        # Try CSV format first (SpliceTransformer often outputs CSV despite .vcf extension)
        scores, labels, metadata = load_splicetransformer_csv_results(args.predictions, args.labels)
        if scores is None:
            # Fallback to VCF format
            scores, labels, metadata = load_splicetransformer_vcf_results(args.predictions, args.labels)
    
    if scores is None:
        logging.error("Failed to load results!")
        return 1
    
    # Compute metrics
    logging.info("Computing performance metrics...")
    metrics = compute_binary_classification_metrics(scores, labels)
    
    # Create plots
    logging.info("Generating performance plots...")
    create_performance_plots(scores, labels, metrics, args.output_dir)
    
    # Generate summary
    summary_file = os.path.join(args.output_dir, 'splicetransformer_evaluation_results.json')
    summary = generate_performance_summary(metrics, metadata, summary_file)
    
    logging.info(f"✅ Evaluation completed! Summary saved to: {summary_file}")
    
    return 0

if __name__ == '__main__':
    exit(main())
