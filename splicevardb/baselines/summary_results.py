#!/usr/bin/env python3
"""
Summary script to compile all model predictions on SpliceVarDB
Extracts variant IDs, ground truth labels, and model prediction scores
into a unified table for analysis.
"""

import pandas as pd
import json
import re
import os
from pathlib import Path

def extract_spliceai_scores(info_field):
    """Extract SpliceAI scores from VCF INFO field"""
    if info_field == '.' or not info_field.startswith('SpliceAI='):
        return None
    # SpliceAI format: ALT|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL
    parts = info_field.replace('SpliceAI=', '').split('|')
    if len(parts) >= 6:
        # Max score from DS_AG, DS_AL, DS_DG, DS_DL
        scores = [float(parts[2]), float(parts[3]), float(parts[4]), float(parts[5])]
        return max(scores)
    return None

def extract_pangolin_scores(info_field):
    """Extract maximum absolute Pangolin score from VCF INFO field"""
    if info_field == '.' or not info_field.startswith('Pangolin='):
        return None
    # Pangolin format: GENE|pos:score|pos:score|...
    pangolin_info = info_field.replace('Pangolin=', '')
    scores = []
    # Extract all position:score pairs
    score_pattern = r'[-+]?\d+:[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?'
    matches = re.findall(score_pattern, pangolin_info)
    for match in matches:
        try:
            score = float(match.split(':')[1])
            scores.append(abs(score))
        except (ValueError, IndexError):
            continue
    
    if scores:
        return max(scores)
    return 0.0

def standardize_variant_id(variant_id):
    """Convert variant IDs to standard format chr:pos:ref>alt"""
    if isinstance(variant_id, str):
        if '-' in variant_id:
            # Convert chr1-100110484-C-A to chr1:100110484:C>A
            parts = variant_id.split('-')
            if len(parts) == 4:
                return f"{parts[0]}:{parts[1]}:{parts[2]}>{parts[3]}"
        elif ':' in variant_id and '>' in variant_id:
            # Already in standard format
            return variant_id
    return variant_id

def main():
    # Define file paths
    base_dir = Path("/Users/byronsun/Desktop/AS_复现模型/BIB_review/code/splicevardb/baselines")
    ground_truth_file = "/Users/byronsun/Desktop/AS_复现模型/BIB_review/data/processed_data/splicevardb/splicevardb_filter_benchmark.tsv"
    
    # Model result files (9 zero-shot models: 4 task-specific + 5 GFMs)
    # MLP models (Evo2_MLP, DNABERT2_MLP) excluded - evaluated on test sets only
    model_files = {
        # Task-specific models
        'Pangolin': base_dir / "pangolin/results/pangolin_splicevardb_full.vcf",
        'SpliceAI': base_dir / "spliceAI/spliceai_predictions_splicevardb.vcf",
        'MMSplice': base_dir / "mmsplice/results/mmsplice_splicevardb_scores.csv",
        'SpliceTransformer': base_dir / "spliceTransformer/SpliceTransformer_splicevardb_predictions.vcf",
        # Genomic Foundation Models
        'AlphaGenome': base_dir / "alphagenome/alphagenome_COMPLETE_results_20250904_094904.json",
        'Evo2_zeroshot': base_dir / "Evo2/results/evo2_predictions_evo2_splicevardb_dataset_dedup.tsv",
        'Nucleotide_Transformer': base_dir / "Nucleotide Transformer/results/nt_splicevardb_scores.csv",
        'SpliceBERT': base_dir / "splicebert/results/splicebert_splicevardb_fixed_scores.csv",
        'DNABERT2_Logistic': base_dir / "DNABert_2/results/dnabert2_variant_scores.tsv"
    }
    
    print("Loading ground truth labels...")
    # Load ground truth
    gt_df = pd.read_csv(ground_truth_file, sep='\t')
    print(f"Loaded {len(gt_df)} ground truth variants")
    
    # Create result dataframe with variant ID and ground truth
    result_df = gt_df[['hg38', 'classification']].copy()
    result_df.columns = ['variant_id', 'ground_truth']
    
    # Standardize variant IDs
    result_df['variant_id'] = result_df['variant_id'].apply(standardize_variant_id)
    print("Standardized variant ID formats")
    
    print("\nProcessing model predictions...")
    
    # 1. AlphaGenome (JSON format)
    print("Processing AlphaGenome...")
    if model_files['AlphaGenome'].exists():
        with open(model_files['AlphaGenome'], 'r') as f:
            alphagenome_data = json.load(f)
        
        ag_scores = {}
        for result in alphagenome_data['results']:
            if result['status'] == 'success':
                variant_id = result['variant_id']
                max_abs_change = result['analysis']['differential']['max_abs_change']
                ag_scores[variant_id] = max_abs_change
        
        result_df['AlphaGenome'] = result_df['variant_id'].map(ag_scores)
        print(f"  - Processed {len(ag_scores)} variants")
    else:
        print(f"  - File not found: {model_files['AlphaGenome']}")
        result_df['AlphaGenome'] = None
    
    # 2. Evo2 zero-shot (TSV format)
    print("Processing Evo2 (zero-shot)...")
    if model_files['Evo2_zeroshot'].exists():
        evo2_df = pd.read_csv(model_files['Evo2_zeroshot'], sep='\t')
        evo2_scores = dict(zip(evo2_df['variant_id'], evo2_df['delta_logp'].abs()))
        result_df['Evo2_zeroshot'] = result_df['variant_id'].map(evo2_scores)
        print(f"  - Processed {len(evo2_scores)} variants")
    else:
        print(f"  - File not found: {model_files['Evo2_zeroshot']}")
        result_df['Evo2_zeroshot'] = None
    
    # 3. Nucleotide Transformer (CSV format)
    print("Processing Nucleotide Transformer...")
    if model_files['Nucleotide_Transformer'].exists():
        nt_df = pd.read_csv(model_files['Nucleotide_Transformer'])
        # Use cosine distance as primary score
        nt_scores = dict(zip(nt_df['variant_id'], nt_df['score_cosine']))
        result_df['Nucleotide_Transformer'] = result_df['variant_id'].map(nt_scores)
        print(f"  - Processed {len(nt_scores)} variants")
    else:
        print(f"  - File not found: {model_files['Nucleotide_Transformer']}")
        result_df['Nucleotide_Transformer'] = None
    
    # 4. Pangolin (VCF format)
    print("Processing Pangolin...")
    if model_files['Pangolin'].exists():
        pangolin_scores = {}
        with open(model_files['Pangolin'], 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 8:
                    chrom, pos, _, ref, alt = parts[:5]
                    info_field = parts[7]
                    variant_id = f"chr{chrom}:{pos}:{ref}>{alt}"
                    score = extract_pangolin_scores(info_field)
                    if score is not None:
                        pangolin_scores[variant_id] = score
        
        result_df['Pangolin'] = result_df['variant_id'].map(pangolin_scores)
        print(f"  - Processed {len(pangolin_scores)} variants")
    else:
        print(f"  - File not found: {model_files['Pangolin']}")
        result_df['Pangolin'] = None
    
    # 5. SpliceAI (VCF format)
    print("Processing SpliceAI...")
    if model_files['SpliceAI'].exists():
        spliceai_scores = {}
        with open(model_files['SpliceAI'], 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 8:
                    chrom, pos, _, ref, alt = parts[:5]
                    info_field = parts[7]
                    variant_id = f"chr{chrom}:{pos}:{ref}>{alt}"
                    score = extract_spliceai_scores(info_field)
                    if score is not None:
                        spliceai_scores[variant_id] = score
        
        result_df['SpliceAI'] = result_df['variant_id'].map(spliceai_scores)
        print(f"  - Processed {len(spliceai_scores)} variants")
    else:
        print(f"  - File not found: {model_files['SpliceAI']}")
        result_df['SpliceAI'] = None
    
    # 6. SpliceBERT (CSV format)
    print("Processing SpliceBERT...")
    if model_files['SpliceBERT'].exists():
        splicebert_df = pd.read_csv(model_files['SpliceBERT'])
        # Use absolute value of KL context score
        splicebert_scores = dict(zip(splicebert_df['variant_id'], splicebert_df['kl_context_score'].abs()))
        result_df['SpliceBERT'] = result_df['variant_id'].map(splicebert_scores)
        print(f"  - Processed {len(splicebert_scores)} variants")
    else:
        print(f"  - File not found: {model_files['SpliceBERT']}")
        result_df['SpliceBERT'] = None
    
    # 7. SpliceTransformer (VCF/CSV hybrid format - check)
    print("Processing SpliceTransformer...")
    if model_files['SpliceTransformer'].exists():
        # Check if file has scores
        with open(model_files['SpliceTransformer'], 'r') as f:
            sample_lines = [line for line in f if not line.startswith('#')][:5]
        
        if sample_lines:
            # Extract scores from the 'score' column (column index 6)
            splicetr_scores = {}
            with open(model_files['SpliceTransformer'], 'r') as f:
                for line in f:
                    if line.startswith('#') or line.startswith(',#CHROM'):
                        continue
                    parts = line.strip().split(',')
                    if len(parts) >= 7:
                        try:
                            chrom = parts[1]
                            pos = parts[2]
                            ref = parts[4]
                            alt = parts[5]
                            score = float(parts[6])
                            variant_id = f"chr{chrom}:{pos}:{ref}>{alt}"
                            splicetr_scores[variant_id] = score
                        except (ValueError, IndexError):
                            continue
            
            result_df['SpliceTransformer'] = result_df['variant_id'].map(splicetr_scores)
            print(f"  - Processed {len(splicetr_scores)} variants")
        else:
            print("  - No data found in SpliceTransformer file")
            result_df['SpliceTransformer'] = None
    else:
        print(f"  - File not found: {model_files['SpliceTransformer']}")
        result_df['SpliceTransformer'] = None

    # 8. DNABERT-2 Logistic (TSV format)
    print("Processing DNABERT-2 (Logistic + diff)...")
    if model_files['DNABERT2_Logistic'].exists():
        try:
            dnabert2_lr_df = pd.read_csv(model_files['DNABERT2_Logistic'], sep='\t')
            # Expect columns: variant_id, prob_splice_altering, label
            lr_scores = dict(zip(dnabert2_lr_df['variant_id'], dnabert2_lr_df['prob_splice_altering']))
            result_df['DNABERT2_Logistic'] = result_df['variant_id'].map(lr_scores)
            print(f"  - Processed {len(lr_scores)} variants")
        except Exception as e:
            print(f"  - Failed to parse DNABERT-2 Logistic file: {e}")
            result_df['DNABERT2_Logistic'] = None
    else:
        print(f"  - File not found: {model_files['DNABERT2_Logistic']}")
        result_df['DNABERT2_Logistic'] = None

    # 9. MMSplice (CSV format) — optional
    print("Processing MMSplice (optional)...")
    if model_files['MMSplice'].exists():
        try:
            mmsplice_df = pd.read_csv(model_files['MMSplice'])
            # Use pathogenicity as primary scalar; aggregate by variant ID taking max absolute
            # Source CSV may include multiple transcript entries per variant ID
            grouped = mmsplice_df.groupby('ID')['pathogenicity'].max()
            mmsplice_scores = grouped.to_dict()
            result_df['MMSplice_pathogenicity'] = result_df['variant_id'].map(mmsplice_scores)
            print(f"  - Processed {len(mmsplice_scores)} variants")
        except Exception as e:
            print(f"  - Failed to parse MMSplice file: {e}")
            result_df['MMSplice_pathogenicity'] = None
    else:
        print(f"  - File not found: {model_files['MMSplice']}")
        result_df['MMSplice_pathogenicity'] = None
    
    # Convert ground truth to binary labels
    result_df['ground_truth_binary'] = (result_df['ground_truth'] == 'Splice-altering').astype(int)
    
    # Save results
    output_file = base_dir / "splicevardb_all_model_predictions.csv"
    result_df.to_csv(output_file, index=False)
    
    print(f"\n✅ Summary complete!")
    print(f"📊 Total variants: {len(result_df)}")
    print(f"📁 Results saved to: {output_file}")
    
    # Print data availability summary (9 zero-shot models: 4 task-specific + 5 GFMs)
    print(f"\n📈 Model prediction coverage (9 zero-shot models):")
    print(f"  Task-Specific Models:")
    task_specific_models = ['Pangolin', 'SpliceAI', 'MMSplice_pathogenicity', 'SpliceTransformer']
    for col in task_specific_models:
        if col in result_df.columns:
            non_null_count = result_df[col].notna().sum()
            coverage = (non_null_count / len(result_df)) * 100
            print(f"    {col}: {non_null_count}/{len(result_df)} ({coverage:.1f}%)")
    
    print(f"  Genomic Foundation Models:")
    gfm_models = ['AlphaGenome', 'Evo2_zeroshot', 'Nucleotide_Transformer', 'SpliceBERT', 'DNABERT2_Logistic']
    for col in gfm_models:
        if col in result_df.columns:
            non_null_count = result_df[col].notna().sum()
            coverage = (non_null_count / len(result_df)) * 100
            print(f"    {col}: {non_null_count}/{len(result_df)} ({coverage:.1f}%)")
    
    print(f"\n  Note: Evo2_MLP and DNABERT2_MLP evaluated on test sets only, reported separately")
    
    # Show ground truth distribution
    gt_dist = result_df['ground_truth'].value_counts()
    print(f"\n🎯 Ground truth distribution:")
    for label, count in gt_dist.items():
        pct = (count / len(result_df)) * 100
        print(f"  - {label}: {count} ({pct:.1f}%)")
    
    # Preview first few rows
    print(f"\n👀 Preview of results:")
    print(result_df.head())

if __name__ == "__main__":
    main()
