#!/usr/bin/env python3
"""
Create imbalanced versions of SpliceVarDB dataset for testing SpliceBERT performance hypothesis.
Strategy: Use all negative samples, subsample positive samples to achieve target imbalance.
"""

import pandas as pd
from pathlib import Path
import numpy as np

def create_imbalanced_datasets(input_file: Path, output_dir: Path) -> pd.DataFrame:
    """
    Create multiple imbalanced versions of SpliceVarDB with different positive rates
    to match the conditions where SpliceBERT was originally validated.
    Strategy: Use all negative samples, subsample positive samples.
    """
    
    # Load original balanced dataset
    df = pd.read_csv(input_file, sep='\t')
    print(f"Original dataset: {len(df)} variants")
    
    # Convert label to binary (splice-altering = 1, normal = 0)
    df['label_binary'] = (df['label'] == 'splice-altering').astype(int)
    print(f"Positive samples: {df['label_binary'].sum()} ({df['label_binary'].mean():.1%})")
    print(f"Negative samples: {len(df) - df['label_binary'].sum()} ({(1-df['label_binary'].mean()):.1%})")
    
    # Separate positive and negative samples
    positives = df[df['label_binary'] == 1].copy()
    negatives = df[df['label_binary'] == 0].copy()
    
    print(f"\nSeparated: {len(positives)} positives, {len(negatives)} negatives")
    print("Strategy: Use all negatives, subsample positives to achieve target imbalance")
    
    # Create single imbalanced dataset: all negatives + 5% of positives
    print(f"\n=== Creating imbalanced dataset ===")
    print("Strategy: Use all negatives + 5% of positives")
    
    # Use all negatives
    n_negatives = len(negatives)
    
    # Use 5% of positives
    positive_fraction = 0.05
    n_positives_subset = int(len(positives) * positive_fraction)
    
    print(f"Using all {n_negatives} negative samples")
    print(f"Using {n_positives_subset} positive samples ({positive_fraction:.1%} of {len(positives)} available)")
    
    # Randomly subsample positives
    positives_subset = positives.sample(n=n_positives_subset, random_state=42)
    
    # Combine positives and all negatives
    imbalanced_df = pd.concat([positives_subset, negatives], ignore_index=True)
    
    # Shuffle the dataset
    imbalanced_df = imbalanced_df.sample(frac=1, random_state=42).reset_index(drop=True)
    
    # Verify the class distribution
    actual_pos_rate = imbalanced_df['label_binary'].mean()
    print(f"Actual positive rate: {actual_pos_rate:.3f}")
    
    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Save the imbalanced dataset
    output_file = output_dir / "splicevardb_imbalanced.tsv"
    imbalanced_df.to_csv(output_file, sep='\t', index=False)
    print(f"Saved: {output_file}")
    
    # Create summary
    summary_data = {
        'description': f'Imbalanced SpliceVarDB: {positive_fraction:.1%} positives',
        'positive_fraction': positive_fraction,
        'actual_pos_rate': actual_pos_rate,
        'total_samples': len(imbalanced_df),
        'positive_samples': imbalanced_df['label_binary'].sum(),
        'negative_samples': len(imbalanced_df) - imbalanced_df['label_binary'].sum(),
        'filename': output_file.name
    }
    
    # Save summary
    summary_df = pd.DataFrame([summary_data])
    summary_file = output_dir / "splicevardb_imbalanced_summary.csv"
    summary_df.to_csv(summary_file, index=False)
    
    print(f"\n=== Summary ===")
    print(f"Description: {summary_data['description']}")
    print(f"Total samples: {summary_data['total_samples']:,}")
    print(f"Positive samples: {summary_data['positive_samples']:,} ({summary_data['actual_pos_rate']:.3f})")
    print(f"Negative samples: {summary_data['negative_samples']:,}")
    print(f"Filename: {summary_data['filename']}")
    print(f"\nSummary saved to: {summary_file}")
    
    return summary_df

def main():
    # Paths
    base_dir = Path("/Users/byronsun/Desktop/AS_复现模型/BIB_review")
    input_file = base_dir / "data/processed_data/splicevardb/evo2_splicevardb_dataset_dedup.tsv"
    output_dir = base_dir / "data/processed_data/splicevardb"
    
    print("Creating imbalanced SpliceVarDB datasets for SpliceBERT validation...")
    print(f"Input: {input_file}")
    print(f"Output directory: {output_dir}")
    
    if not input_file.exists():
        print(f"Error: Input file not found: {input_file}")
        return
    
    summary = create_imbalanced_datasets(input_file, output_dir)
    print("\nDone!")

if __name__ == "__main__":
    main()