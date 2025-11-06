#!/usr/bin/env python3
"""
Create labels file for SpliceTransformer evaluation from ChickenGTEx TSV.
"""

import pandas as pd
import argparse
from pathlib import Path

def create_labels_file(input_tsv, output_csv):
    """
    Extract variant_id and labels from ChickenGTEx TSV for SpliceTransformer evaluation.
    """
    df = pd.read_csv(input_tsv, sep='\t', header=None,
                     names=['variant_id', 'ref_seq', 'alt_seq', 'label', 'tissue_id'])
    
    print(f"Loaded {len(df)} variants from {input_tsv}")
    
    labels_df = df[['variant_id', 'label']].copy()
    labels_df['label'] = labels_df['label'].astype(int)
    
    print(f"Label distribution:")
    print(f"  Positive (1): {labels_df['label'].sum()} ({labels_df['label'].mean():.1%})")
    print(f"  Negative (0): {len(labels_df) - labels_df['label'].sum()} ({1-labels_df['label'].mean():.1%})")
    
    labels_df.to_csv(output_csv, index=False)
    print(f"Saved labels file: {output_csv}")
    
    print("\nSample labels:")
    print(labels_df.head())
    
    return labels_df

def main():
    parser = argparse.ArgumentParser(description="Create labels file for SpliceTransformer evaluation (ChickenGTEx)")
    parser.add_argument('--input_tsv', type=str,
                        default='/Users/byronsun/Desktop/AS_复现模型/BIB_review/data/processed_data/chickenGTEx/chickengtex_silver_benchmark_balanced.tsv',
                        help='Input ChickenGTEx TSV file')
    parser.add_argument('--output_csv', type=str,
                        default='/Users/byronsun/Desktop/AS_复现模型/BIB_review/code/chickenGTEx/baselines/spliceTransformer/chickengtex_labels.csv',
                        help='Output labels CSV file')
    
    args = parser.parse_args()
    
    if not Path(args.input_tsv).exists():
        print(f"Error: Input file not found: {args.input_tsv}")
        return 1
    
    create_labels_file(args.input_tsv, args.output_csv)
    print("\nDone!")
    return 0

if __name__ == "__main__":
    exit(main())

