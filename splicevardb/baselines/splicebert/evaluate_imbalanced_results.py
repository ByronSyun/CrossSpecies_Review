#!/usr/bin/env python3
"""
Quick evaluation of SpliceBERT performance on imbalanced SpliceVarDB dataset
"""

import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score, accuracy_score, precision_score, recall_score, f1_score
from pathlib import Path

def evaluate_imbalanced_results():
    """Evaluate SpliceBERT performance on imbalanced dataset"""
    
    # Paths
    base_dir = Path("/Users/byronsun/Desktop/AS_复现模型/BIB_review")
    
    # Load ground truth (imbalanced dataset)
    gt_file = base_dir / "data/processed_data/splicevardb/splicevardb_imbalanced.tsv"
    df_gt = pd.read_csv(gt_file, sep='\t')
    
    # Load SpliceBERT scores
    scores_file = base_dir / "code/splicevardb/baselines/splicebert/results/splicebert_imbalanced_scores.csv"
    df_scores = pd.read_csv(scores_file)
    
    print("=== SpliceBERT Imbalanced Dataset Evaluation ===")
    print(f"Ground truth samples: {len(df_gt)}")
    print(f"Score samples: {len(df_scores)}")
    
    # Merge data
    merged = df_gt[['variant_id', 'label_binary']].merge(
        df_scores[['variant_id', 'kl_context_score']], 
        on='variant_id', how='inner'
    )
    
    print(f"Merged samples: {len(merged)}")
    
    # Check class distribution
    n_positive = merged['label_binary'].sum()
    n_negative = len(merged) - n_positive
    pos_rate = n_positive / len(merged)
    
    print(f"\nClass distribution:")
    print(f"Positive (splice-altering): {n_positive} ({pos_rate:.3f})")
    print(f"Negative (normal): {n_negative} ({1-pos_rate:.3f})")
    
    # Prepare data for evaluation
    y_true = merged['label_binary'].values
    y_scores = merged['kl_context_score'].values
    
    # SpliceBERT uses negative KL-divergence, so higher (less negative) scores indicate more splice-altering
    # We need to flip the scores for proper evaluation
    y_scores_flipped = -y_scores
    
    # Calculate metrics
    auroc = roc_auc_score(y_true, y_scores_flipped)
    auprc = average_precision_score(y_true, y_scores_flipped)
    
    # Find optimal threshold using F1 score
    from sklearn.metrics import precision_recall_curve
    precision, recall, thresholds = precision_recall_curve(y_true, y_scores_flipped)
    f1_scores = 2 * (precision * recall) / (precision + recall + 1e-8)
    optimal_idx = np.argmax(f1_scores)
    optimal_threshold = thresholds[optimal_idx] if optimal_idx < len(thresholds) else thresholds[-1]
    
    # Calculate threshold-based metrics
    y_pred = (y_scores_flipped >= optimal_threshold).astype(int)
    accuracy = accuracy_score(y_true, y_pred)
    precision_score_val = precision_score(y_true, y_pred, zero_division=0)
    recall_score_val = recall_score(y_true, y_pred, zero_division=0)
    f1_score_val = f1_score(y_true, y_pred, zero_division=0)
    
    print(f"\n=== Performance Metrics ===")
    print(f"AUROC: {auroc:.4f}")
    print(f"AUPRC: {auprc:.4f}")
    print(f"Optimal threshold: {optimal_threshold:.4f}")
    print(f"Accuracy: {accuracy:.4f}")
    print(f"Precision: {precision_score_val:.4f}")
    print(f"Recall: {recall_score_val:.4f}")
    print(f"F1-score: {f1_score_val:.4f}")
    
    # Compare with previous balanced dataset results
    print(f"\n=== Comparison with Balanced Dataset ===")
    print(f"Balanced SpliceVarDB AUROC: 0.5163")
    print(f"Imbalanced SpliceVarDB AUROC: {auroc:.4f}")
    print(f"Improvement: {auroc/0.5163:.2f}x")
    
    if auroc > 0.5163:
        print("✅ HYPOTHESIS CONFIRMED: SpliceBERT performs better on imbalanced data!")
    else:
        print("❌ HYPOTHESIS NOT CONFIRMED: No improvement on imbalanced data")
    
    # Score distribution analysis
    print(f"\n=== Score Distribution ===")
    print(f"Score range: [{y_scores.min():.2f}, {y_scores.max():.2f}]")
    print(f"Mean score: {y_scores.mean():.2f}")
    print(f"Std score: {y_scores.std():.2f}")
    
    pos_scores = y_scores[y_true == 1]
    neg_scores = y_scores[y_true == 0]
    print(f"Positive samples mean score: {pos_scores.mean():.2f} ± {pos_scores.std():.2f}")
    print(f"Negative samples mean score: {neg_scores.mean():.2f} ± {neg_scores.std():.2f}")
    
    return {
        'auroc': auroc,
        'auprc': auprc,
        'accuracy': accuracy,
        'precision': precision_score_val,
        'recall': recall_score_val,
        'f1': f1_score_val,
        'n_samples': len(merged),
        'pos_rate': pos_rate
    }

if __name__ == "__main__":
    results = evaluate_imbalanced_results()
