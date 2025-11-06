#!/usr/bin/env python3
"""
Create balanced benchmark for PigGTEx zero-shot evaluation.
Combines positive samples with equal number of randomly sampled negatives.
Output format matches ratGTEx: variant_id, ref_seq, alt_seq, label, tissue_id (no header)
"""

import pandas as pd
import numpy as np
import logging
import argparse

logging.basicConfig(level=logging.INFO, format='%(asctime)s - [%(levelname)s] - %(message)s')

def create_balanced_benchmark(positive_file, negative_pool_file, output_file, random_seed=42):
    """
    Create balanced benchmark by combining all positives with equal number of negatives.
    """
    logging.info("=" * 80)
    logging.info("Creating PigGTEx Balanced Benchmark for Zero-Shot Evaluation")
    logging.info("=" * 80)
    
    # Load positive samples
    logging.info(f"Loading positive samples from: {positive_file}")
    positive_df = pd.read_csv(
        positive_file, 
        sep='\t', 
        header=None,
        names=['variant_id', 'ref_sequence', 'alt_sequence', 'label', 'tissue_id']
    )
    n_positive = len(positive_df)
    logging.info(f"Loaded {n_positive:,} positive samples")
    
    # Load negative pool
    logging.info(f"Loading negative pool from: {negative_pool_file}")
    negative_pool_df = pd.read_csv(
        negative_pool_file,
        sep='\t',
        header=None,
        names=['variant_id', 'ref_sequence', 'alt_sequence', 'label', 'tissue_id']
    )
    n_negative_pool = len(negative_pool_df)
    logging.info(f"Loaded {n_negative_pool:,} negative samples in pool")
    
    # Check if we have enough negatives
    if n_negative_pool < n_positive:
        logging.error(f"Insufficient negatives: need {n_positive:,}, but only have {n_negative_pool:,}")
        return
    
    # Randomly sample negatives to match positive count
    logging.info(f"Randomly sampling {n_positive:,} negatives from pool (seed={random_seed})")
    np.random.seed(random_seed)
    negative_df = negative_pool_df.sample(n=n_positive, random_state=random_seed)
    
    # Combine positive and negative
    logging.info("Combining positive and negative samples")
    combined_df = pd.concat([positive_df, negative_df], ignore_index=True)
    
    # Shuffle the combined dataset
    logging.info(f"Shuffling combined dataset (seed={random_seed})")
    combined_df = combined_df.sample(frac=1.0, random_state=random_seed).reset_index(drop=True)
    
    # Final statistics
    n_total = len(combined_df)
    n_pos_final = (combined_df['label'] == 1).sum()
    n_neg_final = (combined_df['label'] == 0).sum()
    
    logging.info("=" * 80)
    logging.info("Final Benchmark Statistics")
    logging.info("=" * 80)
    logging.info(f"Total samples: {n_total:,}")
    logging.info(f"Positive (label=1): {n_pos_final:,} ({n_pos_final/n_total:.1%})")
    logging.info(f"Negative (label=0): {n_neg_final:,} ({n_neg_final/n_total:.1%})")
    logging.info(f"Balance ratio: {n_pos_final}/{n_neg_final} = {n_pos_final/n_neg_final:.3f}")
    
    # Save to file (no header, same format as ratGTEx)
    logging.info(f"Saving balanced benchmark to: {output_file}")
    combined_df.to_csv(output_file, sep='\t', index=False, header=False)
    
    logging.info("✅ Balanced benchmark created successfully")
    logging.info("=" * 80)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create balanced PigGTEx benchmark')
    parser.add_argument(
        '--positive', 
        type=str,
        default='./processed_data/pig_positive_samples_HIGH_QUALITY.tsv',
        help='Path to positive samples file'
    )
    parser.add_argument(
        '--negative_pool',
        type=str, 
        default='./processed_data/pig_negative_pool_HIGH_QUALITY.tsv',
        help='Path to negative pool file'
    )
    parser.add_argument(
        '--output',
        type=str,
        default='./processed_data/piggtex_silver_benchmark_balanced.tsv',
        help='Path to output balanced benchmark file'
    )
    parser.add_argument(
        '--seed',
        type=int,
        default=42,
        help='Random seed for sampling and shuffling'
    )
    
    args = parser.parse_args()
    
    create_balanced_benchmark(
        positive_file=args.positive,
        negative_pool_file=args.negative_pool,
        output_file=args.output,
        random_seed=args.seed
    )

