#!/usr/bin/env python3
"""
Evaluate SpliceTransformer Cross-Species Performance on ChickenGTEx
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

logging.basicConfig(level=logging.INFO, format='%(message)s')
plt.style.use('default')
sns.set_palette("husl")

def load_splicetransformer_vcf_results(predictions_vcf, labels_csv):
    try:
        logging.info(f"Loading labels from: {labels_csv}")
        df_labels = pd.read_csv(labels_csv)
        logging.info(f"Loaded {len(df_labels):,} labeled variants")

        logging.info(f"Loading SpliceTransformer predictions from: {predictions_vcf}")

        # Attempt 1: CSV-like file
        try:
            df_predictions = pd.read_csv(predictions_vcf, sep=None, engine='python')
            cols = {c.lower(): c for c in df_predictions.columns}
            chrom_col = cols.get('#chrom') or cols.get('chrom')
            pos_col = cols.get('pos')
            ref_col = cols.get('ref')
            alt_col = cols.get('alt')
            score_col = cols.get('score') or cols.get('spliceai_score') or cols.get('splicetransformer_score')
            
            if chrom_col and pos_col and ref_col and alt_col and score_col:
                tmp = df_predictions[[chrom_col, pos_col, ref_col, alt_col, score_col]].copy()
                tmp.columns = ['CHROM', 'POS', 'REF', 'ALT', 'score']
                tmp['variant_id'] = (tmp['CHROM'].astype(str) + ':' + tmp['POS'].astype(str)
                                     + '_' + tmp['REF'].astype(str) + '/' + tmp['ALT'].astype(str))
                df_predictions = tmp
            else:
                raise ValueError('CSV-like parse did not find expected columns')
        except Exception:
            # Attempt 2: True VCF with SCORE= in INFO
            vcf_data = []
            with open(predictions_vcf, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    fields = line.strip().split('\t')
                    if len(fields) >= 8:
                        chrom, pos, variant_id, ref, alt, qual, filter_field, info = fields[:8]
                        score = None
                        if 'SCORE=' in info:
                            for item in info.split(';'):
                                if item.startswith('SCORE='):
                                    try:
                                        score = float(item.split('=')[1])
                                    except Exception:
                                        score = None
                                    break
                        elif info != '.' and info != '':
                            try:
                                score = float(info)
                            except Exception:
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

        logging.info(f"Parsed {len(df_predictions):,} predictions")
        if len(df_predictions) == 0:
            logging.error("No predictions found (CSV/VCF parse failed)")
            return None, None, None
        
        merged = df_predictions.merge(df_labels[['variant_id', 'label']], on='variant_id', how='inner')
        if len(merged) == 0:
            logging.error("No matching variants found between predictions and labels")
            return None, None, None
        
        scores = merged['score'].values
        labels = merged['label'].values
        metadata = {
            'total_variants': len(scores),
            'positive_samples': int(np.sum(labels)),
            'negative_samples': int(len(labels) - np.sum(labels))
        }
        return scores, labels, metadata
    except Exception as e:
        logging.error(f"Failed to load SpliceTransformer results: {e}")
        return None, None, None

def compute_metrics(scores, labels):
    auroc = roc_auc_score(labels, scores)
    auprc = average_precision_score(labels, scores)
    fpr, tpr, roc_thresholds = roc_curve(labels, scores)
    precision, recall, pr_thresholds = precision_recall_curve(labels, scores)
    
    j_scores = tpr - fpr
    optimal_idx = np.argmax(j_scores)
    optimal_threshold_youden = roc_thresholds[optimal_idx]
    
    f1_scores = 2 * precision * recall / (precision + recall + 1e-10)
    f1_optimal_idx = np.argmax(f1_scores)
    optimal_threshold_f1 = pr_thresholds[f1_optimal_idx] if f1_optimal_idx < len(pr_thresholds) else pr_thresholds[-1]
    
    pred_youden = (scores >= optimal_threshold_youden).astype(int)
    pred_f1 = (scores >= optimal_threshold_f1).astype(int)
    
    tn_y, fp_y, fn_y, tp_y = confusion_matrix(labels, pred_youden).ravel()
    tn_f1, fp_f1, fn_f1, tp_f1 = confusion_matrix(labels, pred_f1).ravel()
    
    return {
        'auroc': auroc,
        'auprc': auprc,
        'optimal_thresholds': {'youden': optimal_threshold_youden, 'f1': optimal_threshold_f1},
        'youden_metrics': {
            'threshold': optimal_threshold_youden,
            'accuracy': accuracy_score(labels, pred_youden),
            'precision': precision_score(labels, pred_youden, zero_division=0),
            'recall': recall_score(labels, pred_youden, zero_division=0),
            'f1': f1_score(labels, pred_youden, zero_division=0),
        },
        'f1_optimal_metrics': {
            'threshold': optimal_threshold_f1,
            'accuracy': accuracy_score(labels, pred_f1),
            'precision': precision_score(labels, pred_f1, zero_division=0),
            'recall': recall_score(labels, pred_f1, zero_division=0),
            'f1': f1_score(labels, pred_f1, zero_division=0),
        },
        'curve_data': {
            'roc': {'fpr': fpr.tolist(), 'tpr': tpr.tolist(), 'thresholds': roc_thresholds.tolist()},
            'pr': {'precision': precision.tolist(), 'recall': recall.tolist(), 'thresholds': pr_thresholds.tolist()}
        }
    }

def main():
    parser = argparse.ArgumentParser(description="Evaluate SpliceTransformer on ChickenGTEx")
    default_predictions = "/Users/byronsun/Desktop/AS_复现模型/BIB_review/code/chickenGTEx/baselines/spliceTransformer/SpliceTransformer_chickengtex_predictions.vcf"
    default_labels = "/Users/byronsun/Desktop/AS_复现模型/BIB_review/code/chickenGTEx/baselines/spliceTransformer/chickengtex_labels.csv"
    default_output_dir = "/Users/byronsun/Desktop/AS_复现模型/BIB_review/code/chickenGTEx/baselines/spliceTransformer"
    
    parser.add_argument('--predictions', type=str, default=default_predictions)
    parser.add_argument('--labels', type=str, default=default_labels)
    parser.add_argument('--output_dir', type=str, default=default_output_dir)
    
    args = parser.parse_args()
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    scores, labels, metadata = load_splicetransformer_vcf_results(args.predictions, args.labels)
    if scores is None:
        return 1
    
    metrics = compute_metrics(scores, labels)
    
    summary_file = os.path.join(args.output_dir, 'splicetransformer_evaluation_results.json')
    with open(summary_file, 'w') as f:
        json.dump({'auroc': metrics['auroc'], 'auprc': metrics['auprc']}, f, indent=2)
    
    print(f"Saved: {summary_file}")
    return 0

if __name__ == '__main__':
    exit(main())

