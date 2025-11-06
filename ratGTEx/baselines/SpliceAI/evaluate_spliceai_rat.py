#!/usr/bin/env python3
"""
Evaluate SpliceAI performance on rat data.
"""

import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score, accuracy_score, precision_score, recall_score, f1_score, precision_recall_curve
import argparse
import os

def standardize_variant_id(variant_id):
    """Convert variant IDs to standard format chr:pos:ref>alt"""
    if isinstance(variant_id, str):
        if '_' in variant_id and '/' in variant_id:
            # Convert 15:55290237_A/G to 15:55290237:A>G
            if ':' in variant_id:
                chrom_pos, ref_alt = variant_id.split('_')
                ref, alt = ref_alt.split('/')
                return f"{chrom_pos}:{ref}>{alt}"
        elif ':' in variant_id and '>' in variant_id:
            # Already in standard format
            return variant_id
    return variant_id

def load_ground_truth(tsv_path):
    """Load rat ground truth data"""
    print(f"Loading ground truth from: {tsv_path}")
    
    # Rat TSV format: variant_id, ref_seq, alt_seq, label, tissue_id
    df = pd.read_csv(tsv_path, sep='\t', header=None, 
                     names=['variant_id', 'ref_seq', 'alt_seq', 'label', 'tissue_id'])
    
    # Standardize variant IDs
    df['variant_id'] = df['variant_id'].apply(standardize_variant_id)
    
    print(f"Loaded {len(df)} ground truth variants")
    print(f"Positive samples: {sum(df['label'] == 1)}")
    print(f"Negative samples: {sum(df['label'] == 0)}")
    
    return df[['variant_id', 'label']].rename(columns={'label': 'ground_truth'})

def load_spliceai_predictions(tsv_path):
    """Load SpliceAI predictions"""
    print(f"Loading SpliceAI predictions from: {tsv_path}")
    
    df = pd.read_csv(tsv_path, sep='\t')
    print(f"Loaded {len(df)} SpliceAI predictions")
    
    # Create standardized variant ID
    df['variant_id'] = df['CHROM'].astype(str) + ':' + df['POS'].astype(str) + ':' + df['REF'] + '>' + df['ALT']
    
    # Use MAX_DS as the prediction score
    df = df[['variant_id', 'MAX_DS']].rename(columns={'MAX_DS': 'spliceai_score'})
    
    print(f"SpliceAI score range: {df['spliceai_score'].min():.4f} to {df['spliceai_score'].max():.4f}")
    print(f"Non-zero scores: {sum(df['spliceai_score'] > 0)}")
    
    return df

def evaluate_performance(ground_truth_df, predictions_df):
    """Evaluate SpliceAI performance"""
    print("\n=== Evaluating SpliceAI Performance ===")
    
    # Merge ground truth with predictions
    merged_df = ground_truth_df.merge(predictions_df, on='variant_id', how='inner')
    print(f"Merged {len(merged_df)} variants (ground truth + predictions)")
    
    if len(merged_df) == 0:
        print("ERROR: No overlapping variants between ground truth and predictions!")
        return None
    
    # Calculate coverage
    coverage = len(merged_df) / len(ground_truth_df) * 100
    print(f"Coverage: {coverage:.2f}% ({len(merged_df)}/{len(ground_truth_df)})")
    
    y_true = merged_df['ground_truth'].values
    y_scores = merged_df['spliceai_score'].values
    
    # Remove any NaN scores
    valid_mask = ~np.isnan(y_scores)
    y_true = y_true[valid_mask]
    y_scores = y_scores[valid_mask]
    
    if len(y_true) == 0:
        print("ERROR: No valid scores after removing NaNs!")
        return None
    
    print(f"Valid predictions: {len(y_true)}")
    
    # Calculate threshold-independent metrics
    auroc = roc_auc_score(y_true, y_scores)
    auprc = average_precision_score(y_true, y_scores)
    
    # Find optimal threshold using F1-score
    precision_curve, recall_curve, thresholds = precision_recall_curve(y_true, y_scores)
    f1_scores = 2 * precision_curve[:-1] * recall_curve[:-1] / (precision_curve[:-1] + recall_curve[:-1] + 1e-8)
    
    if len(f1_scores) > 0:
        optimal_idx = np.argmax(f1_scores)
        optimal_threshold = thresholds[optimal_idx]
    else:
        optimal_threshold = 0.5
    
    # Calculate threshold-dependent metrics
    y_pred = (y_scores >= optimal_threshold).astype(int)
    
    accuracy = accuracy_score(y_true, y_pred)
    precision = precision_score(y_true, y_pred, zero_division=0)
    recall = recall_score(y_true, y_pred, zero_division=0)
    f1 = f1_score(y_true, y_pred, zero_division=0)
    
    # Print results
    print(f"\n=== Results ===")
    print(f"Coverage: {coverage:.2f}%")
    print(f"AUROC: {auroc:.4f}")
    print(f"AUPRC: {auprc:.4f}")
    print(f"Optimal threshold: {optimal_threshold:.4f}")
    print(f"Accuracy: {accuracy:.4f}")
    print(f"Precision: {precision:.4f}")
    print(f"Recall: {recall:.4f}")
    print(f"F1-score: {f1:.4f}")
    
    # Return results dictionary
    results = {
        'model': 'SpliceAI',
        'n_variants': len(merged_df),
        'coverage': coverage,
        'auroc': auroc,
        'auprc': auprc,
        'threshold': optimal_threshold,
        'accuracy': accuracy,
        'precision': precision,
        'recall': recall,
        'f1': f1
    }
    
    return results, merged_df

def main():
    parser = argparse.ArgumentParser(description="Evaluate SpliceAI performance on rat data")
    parser.add_argument('--ground_truth', required=True, help='Path to rat ground truth TSV file')
    parser.add_argument('--predictions', required=True, help='Path to SpliceAI predictions TSV file')
    parser.add_argument('--output', help='Path to save results CSV (optional)')
    
    args = parser.parse_args()
    
    # Load data
    ground_truth_df = load_ground_truth(args.ground_truth)
    predictions_df = load_spliceai_predictions(args.predictions)
    
    # Evaluate performance
    results, merged_df = evaluate_performance(ground_truth_df, predictions_df)
    
    if results is None:
        print("Evaluation failed!")
        return 1
    
    # Save results if output path provided
    if args.output:
        results_df = pd.DataFrame([results])
        results_df.to_csv(args.output, index=False)
        print(f"\nResults saved to: {args.output}")
        
        # Also save the merged data for further analysis
        merged_output = args.output.replace('.csv', '_merged_data.csv')
        merged_df.to_csv(merged_output, index=False)
        print(f"Merged data saved to: {merged_output}")
    
    print("\nEvaluation completed successfully!")
    return 0

if __name__ == "__main__":
    exit(main())
