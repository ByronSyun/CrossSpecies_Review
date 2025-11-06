#!/usr/bin/env python3
"""
Merge chickenGTEx ground truth with model predictions (zero-shot models only).
Outputs a consolidated CSV for downstream evaluation.

Inputs:
- Ground truth (balanced): chickengtex_silver_benchmark_balanced.tsv (no header)
- AlphaGenome JSON: { results: [ { variant_id, alphagenome_score, ... } ] }
- Pangolin JSON: { metadata: {...}, results: [ { variant_id, pangolin_score, ... } ] }
- SpliceAI TSV: CHROM, POS, ID, REF, ALT, GENE, DS_AG, DS_AL, DS_DG, DS_DL, MAX_DS
- SpliceBERT CSV: variant_id, kl_context_score
- Evo2 TSV: variant_id, chrom, pos_0based, ref, alt, label, logp_ref, logp_alt, delta_logp, sequence
- Nucleotide Transformer CSV: variant_id, true_label, cosine_distance
- MMSplice CSV: variant_id, efficiency, delta_logit_psi, pathogenicity
- SpliceTransformer VCF: VCF format with SpliceTransformer scores
- DNABERT2 Logistic TSV: variant_id, prob_splice_altering

Output:
- chickengtex_all_model_predictions.csv with columns:
  variant_id, ground_truth, Pangolin, SpliceAI, MMSplice_pathogenicity, SpliceTransformer, AlphaGenome, Evo2_zeroshot, Nucleotide_Transformer, SpliceBERT, DNABERT2_Logistic
  
Note: Evo2_MLP and DNABERT2_MLP are excluded from this CSV because they are evaluated on test sets only (not all variants)
"""

import json
from pathlib import Path
import pandas as pd
import numpy as np


def normalize_variant_id(vid: str) -> str:
    """
    Normalize variant ID to format 'chrom:pos_ref/alt'
    Handles both 'chrom_pos_ref_alt' and 'chrom:pos_ref/alt' formats
    """
    if pd.isna(vid):
        return vid
    
    # If already in correct format (contains both ':' and '/'), return as is
    if ':' in vid and '/' in vid:
        return vid
    
    # If in format 'chrom_pos_ref_alt', convert to 'chrom:pos_ref/alt'
    if vid.count('_') == 3:
        parts = vid.split('_')
        chrom, pos, ref, alt = parts[0], parts[1], parts[2], parts[3]
        return f"{chrom}:{pos}_{ref}/{alt}"
    
    return vid


def load_ground_truth(path: Path) -> pd.DataFrame:
    # File has 5 columns, no header: variant_id, ref_seq(8192), alt_seq(8192), label, tissue_id
    df = pd.read_csv(path, sep='\t', header=None,
                     names=['variant_id', 'ref_seq', 'alt_seq', 'label', 'tissue_id'],
                     dtype={'variant_id': str})
    df['ground_truth'] = df['label'].astype(int)
    return df[['variant_id', 'ground_truth']]


def load_alphagenome(path: Path) -> pd.DataFrame:
    """Load AlphaGenome results from JSON format (Chicken uses JSON, not TSV)"""
    with open(path, 'r') as f:
        data = json.load(f)
    results = data.get('results', [])
    df = pd.DataFrame(results)
    if 'variant_id' not in df.columns or 'alphagenome_score' not in df.columns:
        print(f"Warning: AlphaGenome JSON missing expected columns")
        return pd.DataFrame(columns=['variant_id', 'AlphaGenome'])
    out = df[['variant_id', 'alphagenome_score']].copy()
    out['variant_id'] = out['variant_id'].apply(normalize_variant_id)  # Normalize format
    out = out.rename(columns={'alphagenome_score': 'AlphaGenome'})
    out = out.drop_duplicates(subset=['variant_id'], keep='first')
    return out


def load_pangolin(path: Path) -> pd.DataFrame:
    data = json.loads(Path(path).read_text())
    results = data.get('results', [])
    df = pd.DataFrame(results)
    assert 'variant_id' in df.columns and 'pangolin_score' in df.columns, 'Pangolin JSON missing columns'
    out = df[['variant_id', 'pangolin_score']].copy().rename(columns={'pangolin_score': 'Pangolin'})
    # Remove duplicates, keep first occurrence
    out = out.drop_duplicates(subset=['variant_id'], keep='first')
    return out


def load_spliceai(path: Path) -> pd.DataFrame:
    """Load SpliceAI results and convert to variant_id format"""
    df = pd.read_csv(path, sep='\t')
    # Create variant_id in the format chr:pos_ref/alt
    df['variant_id'] = df['CHROM'].astype(str) + ':' + df['POS'].astype(str) + '_' + df['REF'] + '/' + df['ALT']
    # Use MAX_DS as the score
    out = df[['variant_id', 'MAX_DS']].copy().rename(columns={'MAX_DS': 'SpliceAI'})
    out = out.drop_duplicates(subset=['variant_id'], keep='first')
    return out


def load_splicebert(path: Path) -> pd.DataFrame:
    """Load SpliceBERT results"""
    df = pd.read_csv(path)
    assert {'variant_id', 'kl_context_score'}.issubset(df.columns), 'SpliceBERT CSV missing columns'
    # Invert score: higher kl_context_score should mean more splice-altering
    out = df[['variant_id', 'kl_context_score']].copy()
    out['SpliceBERT'] = -out['kl_context_score']  # Negate for correct directionality
    out = out[['variant_id', 'SpliceBERT']].drop_duplicates(subset=['variant_id'], keep='first')
    return out


def load_evo2(path: Path) -> pd.DataFrame:
    """Load Evo2 zero-shot results (delta_logp)"""
    df = pd.read_csv(path, sep='\t')
    assert {'variant_id', 'delta_logp'}.issubset(df.columns), 'Evo2 TSV missing columns'
    # Use -delta_logp as the score (higher = more splice-altering)
    out = df[['variant_id', 'delta_logp']].copy()
    out['variant_id'] = out['variant_id'].apply(normalize_variant_id)  # Normalize format
    out['Evo2_zeroshot'] = -out['delta_logp']
    out = out[['variant_id', 'Evo2_zeroshot']].drop_duplicates(subset=['variant_id'], keep='first')
    return out


def load_evo2_mlp(path: Path) -> pd.DataFrame:
    """Load Evo2 MLP results (embedding+MLP)"""
    df = pd.read_csv(path, sep='\t')
    assert {'variant_id', 'prediction_probability'}.issubset(df.columns), 'Evo2 MLP TSV missing columns'
    out = df[['variant_id', 'prediction_probability']].copy().rename(columns={'prediction_probability': 'Evo2_MLP'})
    out = out.drop_duplicates(subset=['variant_id'], keep='first')
    return out


def load_nucleotide_transformer(path: Path) -> pd.DataFrame:
    """Load Nucleotide Transformer results"""
    df = pd.read_csv(path)
    assert {'variant_id', 'score_cosine'}.issubset(df.columns), 'NT CSV missing columns'
    out = df[['variant_id', 'score_cosine']].copy().rename(columns={'score_cosine': 'Nucleotide_Transformer'})
    out = out.drop_duplicates(subset=['variant_id'], keep='first')
    return out


def load_mmsplice(path: Path) -> pd.DataFrame:
    """Load MMSplice results - use pathogenicity score"""
    df = pd.read_csv(path)
    # MMSplice uses 'ID' as variant identifier in format '1:30766214:T>C'
    # Need to convert to '1:30766214_T/C' to match ground truth
    if 'ID' in df.columns:
        # Convert format: '1:30766214:T>C' -> '1:30766214_T/C'
        # Replace last ':' with '_' and '>' with '/'
        def convert_id(id_str):
            parts = id_str.rsplit(':', 1)  # Split from right, only 1 time
            if len(parts) == 2:
                return parts[0] + '_' + parts[1].replace('>', '/')
            return id_str
        df['variant_id'] = df['ID'].apply(convert_id)
    # Group by variant_id and take max absolute pathogenicity (automatically removes duplicates)
    if 'pathogenicity' in df.columns:
        mmsplice_scores = df.groupby('variant_id')['pathogenicity'].max().reset_index()
        mmsplice_scores = mmsplice_scores.rename(columns={'pathogenicity': 'MMSplice_pathogenicity'})
        return mmsplice_scores
    else:
        print('Warning: pathogenicity column not found in MMSplice results')
        return pd.DataFrame(columns=['variant_id', 'MMSplice_pathogenicity'])


def load_dnabert2_logistic(path: Path) -> pd.DataFrame:
    """Load DNABERT-2 Logistic results for chicken data"""
    df = pd.read_csv(path, sep='\t')
    out = df[['variant_id', 'prob_splice_altering']].copy()
    out['variant_id'] = out['variant_id'].apply(normalize_variant_id)  # Normalize format
    out = out.rename(columns={'prob_splice_altering': 'DNABERT2_Logistic'})
    out = out.drop_duplicates(subset=['variant_id'], keep='first')
    return out


def load_dnabert2_mlp(path: Path) -> pd.DataFrame:
    """Load DNABERT-2 MLP results for rat data"""
    df = pd.read_csv(path, sep='\t')
    assert {'variant_id', 'prob_splice_altering'}.issubset(df.columns), 'DNABERT-2 MLP TSV missing columns'
    out = df[['variant_id', 'prob_splice_altering']].copy().rename(columns={'prob_splice_altering': 'DNABERT2_MLP'})
    out = out.drop_duplicates(subset=['variant_id'], keep='first')
    return out


def load_splicetransformer(path: Path) -> pd.DataFrame:
    """Load SpliceTransformer results from VCF/CSV file"""
    df = pd.read_csv(path)
    # Create variant_id in format 'CHROM:POS_REF/ALT'
    df['variant_id'] = df['#CHROM'].astype(str) + ':' + df['POS'].astype(str) + '_' + df['REF'] + '/' + df['ALT']
    out = df[['variant_id', 'score']].copy().rename(columns={'score': 'SpliceTransformer'})
    out = out.drop_duplicates(subset=['variant_id'], keep='first')
    return out


def main():
    base_dir = Path('/Users/byronsun/Desktop/AS_复现模型/BIB_review')

    gt_file = base_dir / 'data/processed_data/chickenGTEx/chickengtex_silver_benchmark_balanced.tsv'
    ag_file = base_dir / 'code/chickenGTEx/baselines/alphagenome/results/chickengtex_alphagenome_complete_results_20251018_205155.json'
    pg_file = base_dir / 'code/chickenGTEx/baselines/pangolin/results/pangolin_chickengtex_results.json'
    sa_file = base_dir / 'code/chickenGTEx/baselines/SpliceAI/results/spliceai_parsed_scores_chicken.tsv'
    evo2_file = base_dir / 'code/chickenGTEx/baselines/Evo2/results/evo2_chicken_predictions.tsv'
    mmsplice_file = base_dir / 'code/chickenGTEx/baselines/mmsplice/results/mmsplice_chickengtex_scores.csv'
    dnabert2_log_file = base_dir / 'code/chickenGTEx/baselines/DNABert_2/results/dnabert2_chicken_logistic_scores.tsv'
    st_file = base_dir / 'code/chickenGTEx/baselines/spliceTransformer/results/SpliceTransformer_chickengtex_predictions.vcf'
    out_csv = base_dir / 'code/chickenGTEx/baselines/chickengtex_all_model_predictions.csv'

    print('Loading ground truth...')
    gt = load_ground_truth(gt_file)
    print(f'  Ground truth variants: {len(gt):,}')

    # Load all models with error handling
    models_data = {}
    
    print('Loading AlphaGenome predictions...')
    try:
        models_data['AlphaGenome'] = load_alphagenome(ag_file)
        print(f'  AlphaGenome predictions: {len(models_data["AlphaGenome"]):,}')
    except Exception as e:
        print(f'  AlphaGenome error: {e}')
        models_data['AlphaGenome'] = pd.DataFrame(columns=['variant_id', 'AlphaGenome'])

    print('Loading Pangolin predictions...')
    try:
        models_data['Pangolin'] = load_pangolin(pg_file)
        print(f'  Pangolin predictions: {len(models_data["Pangolin"]):,}')
    except Exception as e:
        print(f'  Pangolin error: {e}')
        models_data['Pangolin'] = pd.DataFrame(columns=['variant_id', 'Pangolin'])

    print('Loading SpliceAI predictions...')
    try:
        models_data['SpliceAI'] = load_spliceai(sa_file)
        print(f'  SpliceAI predictions: {len(models_data["SpliceAI"]):,}')
    except Exception as e:
        print(f'  SpliceAI error: {e}')
        models_data['SpliceAI'] = pd.DataFrame(columns=['variant_id', 'SpliceAI'])

    # SpliceBERT skipped - poor performance on human benchmark (AUROC 0.5163)
    
    print('Loading Evo2 zero-shot predictions...')
    try:
        models_data['Evo2_zeroshot'] = load_evo2(evo2_file)
        print(f'  Evo2 zero-shot predictions: {len(models_data["Evo2_zeroshot"]):,}')
    except Exception as e:
        print(f'  Evo2 zero-shot error: {e}')
        models_data['Evo2_zeroshot'] = pd.DataFrame(columns=['variant_id', 'Evo2_zeroshot'])

    # Note: Evo2_MLP excluded - evaluated on test set only, reported separately
    
    # Nucleotide_Transformer skipped - poor performance on human benchmark (AUROC 0.5073)

    print('Loading MMSplice predictions...')
    try:
        models_data['MMSplice_pathogenicity'] = load_mmsplice(mmsplice_file)
        print(f'  MMSplice predictions: {len(models_data["MMSplice_pathogenicity"]):,}')
    except Exception as e:
        print(f'  MMSplice error: {e}')
        models_data['MMSplice_pathogenicity'] = pd.DataFrame(columns=['variant_id', 'MMSplice_pathogenicity'])

    print('Loading DNABERT-2 Logistic predictions...')
    try:
        models_data['DNABERT2_Logistic'] = load_dnabert2_logistic(dnabert2_log_file)
        print(f'  DNABERT-2 Logistic predictions: {len(models_data["DNABERT2_Logistic"]):,}')
    except Exception as e:
        print(f'  DNABERT-2 Logistic error: {e}')
        models_data['DNABERT2_Logistic'] = pd.DataFrame(columns=['variant_id', 'DNABERT2_Logistic'])

    # Note: DNABERT2_MLP excluded - evaluated on test set only, reported separately

    print('Loading SpliceTransformer predictions...')
    try:
        models_data['SpliceTransformer'] = load_splicetransformer(st_file)
        print(f'  SpliceTransformer predictions: {len(models_data["SpliceTransformer"]):,}')
    except Exception as e:
        print(f'  SpliceTransformer error: {e}')
        models_data['SpliceTransformer'] = pd.DataFrame(columns=['variant_id', 'SpliceTransformer'])

    print('Merging...')
    df = gt.copy()
    for model_name, model_df in models_data.items():
        df = df.merge(model_df, on='variant_id', how='left')

    # Coverage report (7 zero-shot models: 4 task-specific + 3 GFMs)
    # Note: Nucleotide_Transformer and SpliceBERT skipped due to poor human performance
    total = len(df)
    print(f'\nCoverage Report (7 zero-shot models):')
    print(f'  Task-Specific Models:')
    for model in ['Pangolin', 'SpliceAI', 'MMSplice_pathogenicity', 'SpliceTransformer']:
        if model in df.columns:
            cov = df[model].notna().sum() / total * 100
            print(f'    {model:<28}: {cov:>6.2f}%')
    
    print(f'  Genomic Foundation Models:')
    for model in ['AlphaGenome', 'Evo2_zeroshot', 'DNABERT2_Logistic']:
        if model in df.columns:
            cov = df[model].notna().sum() / total * 100
            print(f'    {model:<28}: {cov:>6.2f}%')
    
    print(f'\n  Note: Nucleotide_Transformer and SpliceBERT skipped (poor human performance)')
    print(f'        Evo2_MLP and DNABERT2_MLP evaluated on test sets only, reported separately')

    print(f'\nSaving consolidated table -> {out_csv}')
    df.to_csv(out_csv, index=False)

    # Small preview
    print('\nPreview:')
    print(df.head(3).to_string(index=False))


if __name__ == '__main__':
    main()


