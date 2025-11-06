#!/usr/bin/env python3
"""
Analyze model prediction coverage on SpliceVarDB
Check for missing predictions and understand the reasons
"""

import pandas as pd
import numpy as np

def analyze_coverage(df, model_name):
    """Analyze prediction coverage for a specific model"""
    total_variants = len(df)
    
    if model_name not in df.columns:
        return None
    
    scores = df[model_name]
    
    # Count different types of missing/invalid predictions
    n_predictions = scores.notna().sum()
    n_missing = scores.isna().sum()
    n_zero = (scores == 0).sum()
    n_valid_nonzero = ((scores != 0) & scores.notna()).sum()
    
    coverage_pct = (n_predictions / total_variants) * 100
    
    result = {
        'model': model_name,
        'total_variants': total_variants,
        'predictions': n_predictions,
        'missing': n_missing,
        'zero_scores': n_zero,
        'valid_nonzero': n_valid_nonzero,
        'coverage_pct': coverage_pct
    }
    
    return result

def main():
    # Load data
    data_file = "/Users/byronsun/Desktop/AS_复现模型/BIB_review/code/splicevardb/baselines/splicevardb_all_model_predictions.csv"
    df = pd.read_csv(data_file)
    
    print("🔍 Model Prediction Coverage Analysis")
    print("=" * 80)
    print(f"Total variants in SpliceVarDB: {len(df)}")
    print()
    
    models = ['AlphaGenome', 'Evo2', 'Nucleotide_Transformer', 'Pangolin', 
              'SpliceAI', 'SpliceBERT', 'SpliceTransformer']
    
    coverage_results = []
    
    for model in models:
        if model in df.columns:
            result = analyze_coverage(df, model)
            if result:
                coverage_results.append(result)
                
                print(f"📊 {model}")
                print(f"  Total predictions: {result['predictions']:,}/{result['total_variants']:,} ({result['coverage_pct']:.2f}%)")
                print(f"  Missing predictions: {result['missing']:,}")
                print(f"  Zero scores: {result['zero_scores']:,}")
                print(f"  Valid non-zero: {result['valid_nonzero']:,}")
                
                # Check for specific issues
                if result['missing'] > 0:
                    print(f"  ⚠️  {result['missing']:,} variants have no predictions")
                    
                    # Sample missing variants
                    missing_mask = df[model].isna()
                    if missing_mask.sum() > 0:
                        missing_examples = df[missing_mask]['variant_id'].head(5).tolist()
                        print(f"     Examples: {', '.join(missing_examples[:3])}...")
                
                if result['zero_scores'] > (result['total_variants'] * 0.1):  # More than 10% zeros
                    print(f"  ⚠️  High proportion of zero scores ({result['zero_scores']:,})")
                
                print()
    
    # Create summary table
    if coverage_results:
        coverage_df = pd.DataFrame(coverage_results)
        
        print("📈 Coverage Summary Table")
        print("=" * 80)
        print(coverage_df[['model', 'predictions', 'missing', 'coverage_pct']].to_string(index=False))
        print()
        
        # Identify models with incomplete coverage
        incomplete_models = coverage_df[coverage_df['coverage_pct'] < 100]
        
        if len(incomplete_models) > 0:
            print("⚠️  Models with Incomplete Coverage:")
            print("=" * 50)
            for _, row in incomplete_models.iterrows():
                print(f"{row['model']}: {row['predictions']:,}/{row['total_variants']:,} ({row['coverage_pct']:.2f}%)")
                print(f"  Missing: {row['missing']:,} variants")
            print()
            
            # Analyze missing patterns
            print("🔍 Missing Prediction Analysis:")
            print("=" * 40)
            
            for model in incomplete_models['model'].values:
                print(f"\n{model} missing predictions:")
                missing_mask = df[model].isna()
                missing_variants = df[missing_mask]
                
                # Check if missing patterns correlate with specific characteristics
                if len(missing_variants) > 0:
                    # By chromosome
                    missing_by_chr = missing_variants['variant_id'].str.extract(r'(chr\w+):').value_counts().head()
                    print(f"  Top chromosomes with missing predictions:")
                    for chr_name, count in missing_by_chr.items():
                        print(f"    {chr_name}: {count} variants")
                    
                    # By ground truth label
                    missing_by_label = missing_variants['ground_truth'].value_counts()
                    print(f"  Missing by ground truth:")
                    for label, count in missing_by_label.items():
                        pct = (count / len(missing_variants)) * 100
                        print(f"    {label}: {count} ({pct:.1f}%)")
        else:
            print("✅ All models have complete coverage (100%)")
    
    # Special analysis for zero scores
    print("\n🎯 Zero Score Analysis:")
    print("=" * 30)
    
    for model in models:
        if model in df.columns:
            scores = df[model].dropna()
            if len(scores) > 0:
                n_zeros = (scores == 0).sum()
                zero_pct = (n_zeros / len(scores)) * 100
                
                if zero_pct > 5:  # Report if >5% zeros
                    print(f"{model}: {n_zeros:,} zero scores ({zero_pct:.1f}%)")
                    
                    # Check if zeros correlate with ground truth
                    zero_mask = (df[model] == 0) & df[model].notna()
                    zero_variants = df[zero_mask]
                    if len(zero_variants) > 0:
                        zero_by_label = zero_variants['ground_truth'].value_counts()
                        print(f"  Zero scores by ground truth:")
                        for label, count in zero_by_label.items():
                            pct = (count / len(zero_variants)) * 100
                            print(f"    {label}: {count} ({pct:.1f}%)")
    
    print("\n📝 Summary for Discussion:")
    print("=" * 30)
    print("Coverage Issues to Discuss:")
    
    for _, row in coverage_df.iterrows():
        if row['coverage_pct'] < 100:
            missing_pct = (row['missing'] / row['total_variants']) * 100
            print(f"- {row['model']}: {missing_pct:.1f}% variants lack predictions ({row['missing']:,}/{row['total_variants']:,})")
            
            # Suggest reasons based on model type
            if row['model'] == 'Pangolin':
                print(f"  Likely reason: Gene region filtering (only processes variants in annotated genes)")
            elif row['model'] == 'SpliceAI':
                print(f"  Likely reason: Processing limitations or coordinate issues")
            elif 'Transformer' in row['model']:
                print(f"  Likely reason: Sequence extraction or tokenization issues")
        else:
            print(f"- {row['model']}: Complete coverage (100%)")

if __name__ == "__main__":
    main()
