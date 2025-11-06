#!/usr/bin/env python3
"""
Plot PigGTEx performance results for all models
Generate publication-ready figures mirroring RatGTEx and SpliceVarDB styles
Colors are matched with SpliceVarDB plots for consistency
"""

import os
import json
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import (
    roc_curve,
    precision_recall_curve,
    roc_auc_score,
    average_precision_score,
)

warnings.filterwarnings("ignore")

# Style
plt.style.use("seaborn-v0_8-whitegrid")
sns.set_palette("husl")

# MODEL COLORS - must match human/rat plots
MODEL_COLORS = {
    'AlphaGenome': '#1f77b4',          # Blue
    'Evo2_zeroshot': '#ff7f0e',        # Orange
    'Nucleotide_Transformer': '#2ca02c', # Green
    'Pangolin': '#d62728',              # Red
    'SpliceAI': '#9467bd',              # Purple
    'SpliceBERT': '#8c564b',            # Brown
    'SpliceTransformer': '#e377c2',     # Pink
    'DNABERT2_Logistic': '#17becf',     # Cyan
    'MMSplice_pathogenicity': '#bcbd22' # Olive
}

MODEL_LABELS = {
    'AlphaGenome': 'AlphaGenome',
    'Evo2_zeroshot': 'Evo2',
    'Nucleotide_Transformer': 'Nucleotide Transformer',
    'Pangolin': 'Pangolin',
    'SpliceAI': 'SpliceAI',
    'SpliceBERT': 'SpliceBERT',
    'SpliceTransformer': 'SpliceTransformer',
    'DNABERT2_Logistic': 'DNABERT-2 (Logistic)',
    'MMSplice_pathogenicity': 'MMSplice'
}


def read_labels(labels_path: str) -> pd.DataFrame:
    df = pd.read_csv(
        labels_path,
        sep='\t', header=None,
        names=['variant_id', 'ref_seq', 'alt_seq', 'label', 'chrom']
    )
    return df[['variant_id', 'label']].drop_duplicates('variant_id')


def load_alphagenome(json_path: str) -> pd.DataFrame:
    if not os.path.exists(json_path):
        return pd.DataFrame(columns=['variant_id', 'AlphaGenome'])
    with open(json_path) as f:
        data = json.load(f)
    rows = data.get('results', [])
    df = pd.DataFrame(rows)
    # Expect columns: variant_id, alphagenome_score
    score_col = 'alphagenome_score' if 'alphagenome_score' in df.columns else 'score'
    out = df[['variant_id', score_col]].rename(columns={score_col: 'AlphaGenome'})
    return out.drop_duplicates('variant_id')


def load_pangolin(json_path: str) -> pd.DataFrame:
    if not os.path.exists(json_path):
        return pd.DataFrame(columns=['variant_id', 'Pangolin'])
    with open(json_path) as f:
        data = json.load(f)
    rows = data.get('results', [])
    df = pd.DataFrame(rows)
    score_col = 'pangolin_score' if 'pangolin_score' in df.columns else 'score'
    out = df[['variant_id', score_col]].rename(columns={score_col: 'Pangolin'})
    return out.drop_duplicates('variant_id')


def load_spliceai(tsv_path: str) -> pd.DataFrame:
    if not os.path.exists(tsv_path):
        return pd.DataFrame(columns=['variant_id', 'SpliceAI'])
    df = pd.read_csv(tsv_path, sep='\t')
    # Build variant_id if needed
    if 'variant_id' not in df.columns and {
        'CHROM', 'POS', 'REF', 'ALT'
    }.issubset(df.columns):
        df['variant_id'] = (
            df['CHROM'].astype(str) + ':' + df['POS'].astype(str) + '_' +
            df['REF'].astype(str) + '/' + df['ALT'].astype(str)
        )
    # Choose score
    if 'MAX_DS' in df.columns:
        df['SpliceAI'] = df['MAX_DS'].astype(float)
        return df[['variant_id', 'SpliceAI']].drop_duplicates('variant_id')
    for col in ['score', 'max_ds', 'delta_score', 'DS_MAX']:
        if col in df.columns:
            return df[['variant_id', col]].rename(columns={col: 'SpliceAI'}).drop_duplicates('variant_id')
    # Fallback: take max of DS_* columns
    ds_cols = [c for c in df.columns if c.upper() in {'DS_AG', 'DS_AL', 'DS_DG', 'DS_DL'}]
    if ds_cols:
        df['SpliceAI'] = df[ds_cols].abs().max(axis=1)
        return df[['variant_id', 'SpliceAI']].drop_duplicates('variant_id')
    return pd.DataFrame(columns=['variant_id', 'SpliceAI'])


def load_nt(csv_path: str) -> pd.DataFrame:
    if not os.path.exists(csv_path):
        return pd.DataFrame(columns=['variant_id', 'Nucleotide_Transformer'])
    df = pd.read_csv(csv_path)
    score_col = 'score_cosine' if 'score_cosine' in df.columns else 'score'
    out = df[['variant_id', score_col]].rename(columns={score_col: 'Nucleotide_Transformer'})
    return out.drop_duplicates('variant_id')


def load_splicebert(csv_path: str) -> pd.DataFrame:
    if not os.path.exists(csv_path):
        return pd.DataFrame(columns=['variant_id', 'SpliceBERT'])
    df = pd.read_csv(csv_path)
    # Expect kl_context_score; keep sign as-is, inversion during metric calc/plot
    for col in ['kl_context_score', 'score']:
        if col in df.columns:
            return df[['variant_id', col]].rename(columns={col: 'SpliceBERT'}).drop_duplicates('variant_id')
    return pd.DataFrame(columns=['variant_id', 'SpliceBERT'])


def load_evo2_zeroshot(tsv_path: str) -> pd.DataFrame:
    if not os.path.exists(tsv_path):
        return pd.DataFrame(columns=['variant_id', 'Evo2_zeroshot'])
    df = pd.read_csv(tsv_path, sep='\t')
    # Expect delta_logp; invert sign so higher = more splice-altering
    col = 'delta_logp' if 'delta_logp' in df.columns else 'score'
    df['Evo2_zeroshot'] = -df[col]
    return df[['variant_id', 'Evo2_zeroshot']].drop_duplicates('variant_id')


def load_evo2_mlp(tsv_path: str) -> pd.DataFrame:
    # Optional: treat as alternative Evo2 usage; plotted as Evo2 if present? We keep as Evo2_MLP (not in head plots)
    if not os.path.exists(tsv_path):
        return pd.DataFrame(columns=['variant_id', 'Evo2_MLP'])
    df = pd.read_csv(tsv_path, sep='\t')
    for col in ['score', 'prob', 'proba', 'mlp_score']:
        if col in df.columns:
            out = df[['variant_id', col]].rename(columns={col: 'Evo2_MLP'})
            return out.drop_duplicates('variant_id')
    return pd.DataFrame(columns=['variant_id', 'Evo2_MLP'])


def load_dnabert2_logistic(tsv_path: str) -> pd.DataFrame:
    if not os.path.exists(tsv_path):
        return pd.DataFrame(columns=['variant_id', 'DNABERT2_Logistic'])
    df = pd.read_csv(tsv_path, sep='\t|,', engine='python')
    for col in ['proba', 'prob', 'score']:
        if col in df.columns:
            return df[['variant_id', col]].rename(columns={col: 'DNABERT2_Logistic'}).drop_duplicates('variant_id')
    # If only logits exist, map to 0.5 default
    if 'variant_id' in df.columns and 'DNABERT2_Logistic' not in df.columns:
        df['DNABERT2_Logistic'] = 0.5
        return df[['variant_id', 'DNABERT2_Logistic']].drop_duplicates('variant_id')
    return pd.DataFrame(columns=['variant_id', 'DNABERT2_Logistic'])


def load_splicetransformer(pred_path: str) -> pd.DataFrame:
    if not os.path.exists(pred_path):
        return pd.DataFrame(columns=['variant_id', 'SpliceTransformer'])
    # First try CSV-like
    try:
        df = pd.read_csv(pred_path, sep='\t|,', engine='python')
        # Try to infer columns
        if 'score' in df.columns and {'#CHROM','POS','REF','ALT'}.issubset(df.columns):
            df['variant_id'] = df['#CHROM'].astype(str) + ':' + df['POS'].astype(str) + '_' + df['REF'].astype(str) + '/' + df['ALT'].astype(str)
            return df[['variant_id', 'score']].rename(columns={'score': 'SpliceTransformer'}).drop_duplicates('variant_id')
    except Exception:
        pass
    # Fallback to VCF parsing
    records = []
    with open(pred_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue
            chrom, pos, _id, ref, alt, _qual, _filter, info = parts[:8]
            score = None
            if 'score=' in info:
                try:
                    score = float(info.split('score=')[1].split(';')[0])
                except Exception:
                    score = None
            if score is not None:
                vid = f"{chrom}:{pos}_{ref}/{alt}"
                records.append((vid, score))
    if not records:
        return pd.DataFrame(columns=['variant_id', 'SpliceTransformer'])
    df = pd.DataFrame(records, columns=['variant_id', 'SpliceTransformer'])
    return df.drop_duplicates('variant_id')


def load_mmsplice(csv_path: str) -> pd.DataFrame:
    if not os.path.exists(csv_path):
        return pd.DataFrame(columns=['variant_id', 'MMSplice_pathogenicity'])
    df = pd.read_csv(csv_path)
    # ID like 1:30766214:T>C -> 1:30766214_T/C
    if 'ID' in df.columns:
        ids = df['ID'].astype(str)
        ids = ids.apply(lambda s: s.rsplit(':', 1)[0].replace('>', '/').replace(':', '_', 1) if ':' in s else s)
        # Correction: rsplit then replace last ':' with '_' properly
        def mm_id(s: str) -> str:
            # split into chrom:pos and ref>alt
            if s.count(':') >= 2 and '>' in s:
                chrom_pos, refalt = s.rsplit(':', 1)
                return chrom_pos.replace(':', ':') + '_' + refalt.replace('>', '/')
            return s
        ids = df['ID'].astype(str).apply(mm_id)
        df['variant_id'] = ids
    if 'pathogenicity' in df.columns:
        out = df[['variant_id', 'pathogenicity']].rename(columns={'pathogenicity': 'MMSplice_pathogenicity'})
        return out.drop_duplicates('variant_id')
    # Fallback: delta_logit_psi absolute
    if 'delta_logit_psi' in df.columns:
        df['MMSplice_pathogenicity'] = df['delta_logit_psi'].abs()
        return df[['variant_id', 'MMSplice_pathogenicity']].drop_duplicates('variant_id')
    return pd.DataFrame(columns=['variant_id', 'MMSplice_pathogenicity'])


def build_merged_dataframe(labels_path: str) -> pd.DataFrame:
    labels = read_labels(labels_path)

    base = labels.rename(columns={'label': 'ground_truth'}).copy()

    # Paths
    ROOT = "/Users/byronsun/Desktop/AS_复现模型/BIB_review/code/pigGTEx/baselines"

    alphagenome = load_alphagenome(f"{ROOT}/alphagenome/results/piggtex_alphagenome_complete_results_20251019_073042.json")
    pangolin = load_pangolin(f"{ROOT}/Pangolin/results/pangolin_piggtex_results.json")
    spliceai = load_spliceai(f"{ROOT}/SpliceAI/results/spliceai_parsed_scores_pig.tsv")
    nt = load_nt(f"{ROOT}/Nucleotide_Transformer/results/nt_pig_scores.csv")
    splicebert = load_splicebert(f"{ROOT}/SpliceBERT/results/splicebert_pig_scores.csv")
    evo2_zero = load_evo2_zeroshot(f"{ROOT}/Evo2/results/evo2_pig_predictions.tsv")
    evo2_mlp = load_evo2_mlp(f"{ROOT}/Evo2/results/evo2_mlp_pig_predictions.tsv")
    dnabert2_log = load_dnabert2_logistic(f"{ROOT}/DNABert_2/results/dnabert2_pig_logistic_scores.tsv")
    splicetransformer = load_splicetransformer(f"{ROOT}/spliceTransformer/results/SpliceTransformer_piggtex_predictions.vcf")
    mmsplice = load_mmsplice(f"{ROOT}/mmsplice/results/mmsplice_piggtex_scores.csv")

    # Merge
    dfs = [
        alphagenome, evo2_zero, nt, pangolin, spliceai, splicebert,
        splicetransformer, dnabert2_log, mmsplice
    ]
    merged = base
    for d in dfs:
        if d is None or d.empty:
            continue
        merged = merged.merge(d, on='variant_id', how='left')

    # Drop duplicates just in case
    merged = merged.drop_duplicates('variant_id')
    return merged


def plot_roc_curves(df: pd.DataFrame, output_path: str) -> None:
    """Plot ROC curves for all models (excluding SpliceBERT, NT, DNABERT2-Logistic which performed poorly on human, and MMSplice due to catastrophic coverage collapse)"""
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    
    # Exclude SpliceBERT, Nucleotide_Transformer, DNABERT2_Logistic (poor human performance)
    # Exclude MMSplice_pathogenicity (coverage collapse 5-12% makes metrics uninterpretable)
    models = ['AlphaGenome', 'Evo2_zeroshot', 'Pangolin', 
              'SpliceAI', 'SpliceTransformer']
    
    y_true = df['ground_truth'].values
    for model in models:
        if model in df.columns:
            y_scores = df[model].values
            mask = ~np.isnan(y_scores)
            y_true_clean = y_true[mask]
            y_scores_clean = y_scores[mask]
            if len(y_true_clean) > 0 and len(np.unique(y_true_clean)) > 1:
                try:
                    fpr, tpr, _ = roc_curve(y_true_clean, y_scores_clean)
                    auroc = roc_auc_score(y_true_clean, y_scores_clean)
                    ax.plot(
                        fpr, tpr, color=MODEL_COLORS[model], linewidth=4.0,
                        label=f"{MODEL_LABELS[model]} (AUROC = {auroc:.3f})"
                    )
                except Exception:
                    pass
    ax.plot([0, 1], [0, 1], 'k--', alpha=0.5, linewidth=2.5)
    ax.set_xlabel('False Positive Rate', fontsize=22, fontweight='bold')
    ax.set_ylabel('True Positive Rate', fontsize=22, fontweight='bold')
    ax.set_title('ROC Curves - PigGTEx Cross-Species Benchmark', fontsize=24, fontweight='bold')
    ax.legend(loc='lower right', fontsize=20, frameon=True, fancybox=True, shadow=True, framealpha=0.95)
    ax.tick_params(axis='both', which='major', labelsize=18)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{output_path}/pig_roc_curves.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_path}/pig_roc_curves.pdf", bbox_inches='tight')
    plt.close()


def plot_precision_recall_curves(df: pd.DataFrame, output_path: str) -> None:
    """Plot Precision-Recall curves for all models (excluding SpliceBERT, NT, DNABERT2-Logistic which performed poorly on human, and MMSplice due to catastrophic coverage collapse)"""
    fig, ax = plt.subplots(1, 1, figsize=(14, 11))
    
    # Exclude SpliceBERT, Nucleotide_Transformer, DNABERT2_Logistic (poor human performance)
    # Exclude MMSplice_pathogenicity (coverage collapse 5-12% makes metrics uninterpretable)
    models = ['AlphaGenome', 'Evo2_zeroshot', 'Pangolin', 
              'SpliceAI', 'SpliceTransformer']
    
    y_true = df['ground_truth'].values
    baseline_precision = y_true.mean()
    for model in models:
        if model in df.columns:
            y_scores = df[model].values
            mask = ~np.isnan(y_scores)
            y_true_clean = y_true[mask]
            y_scores_clean = y_scores[mask]
            if len(y_true_clean) > 0 and len(np.unique(y_true_clean)) > 1:
                try:
                    precision, recall, _ = precision_recall_curve(y_true_clean, y_scores_clean)
                    auprc = average_precision_score(y_true_clean, y_scores_clean)
                    ax.plot(
                        recall, precision, color=MODEL_COLORS[model], linewidth=4.0,
                        label=f"{MODEL_LABELS[model]} (AUPRC = {auprc:.3f})"
                    )
                except Exception:
                    pass
    # Plot baseline (random classifier) - no label in legend
    ax.axhline(y=baseline_precision, color='k', linestyle='--', alpha=0.5, linewidth=2.5)
    
    ax.set_xlabel('Recall', fontsize=26, fontweight='bold')
    ax.set_ylabel('Precision', fontsize=26, fontweight='bold')
    ax.set_title('Precision-Recall Curves - PigGTEx Cross-Species Benchmark', fontsize=28, fontweight='bold')
    ax.legend(loc='lower right', fontsize=28, frameon=True, fancybox=True, shadow=True, framealpha=0.95)
    ax.tick_params(axis='both', which='major', labelsize=22, width=2)
    ax.grid(True, alpha=0.3)
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1.05])
    plt.tight_layout()
    plt.savefig(f"{output_path}/pig_precision_recall_curves.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_path}/pig_precision_recall_curves.pdf", bbox_inches='tight')
    plt.close()


def plot_performance_bar_chart(df: pd.DataFrame, output_path: str) -> None:
    # Exclude SpliceBERT, Nucleotide_Transformer, DNABERT2_Logistic (poor human performance)
    models = ['AlphaGenome', 'Evo2_zeroshot', 'Pangolin', 
              'SpliceAI', 'SpliceTransformer', 'MMSplice_pathogenicity']
    y_true = df['ground_truth'].values
    auroc_scores, auprc_scores, names = [], [], []
    for model in models:
        if model in df.columns:
            y_scores = df[model].values
            mask = ~np.isnan(y_scores)
            y_true_clean = y_true[mask]
            y_scores_clean = y_scores[mask]
            if len(y_true_clean) > 0 and len(np.unique(y_true_clean)) > 1:
                try:
                    auroc_scores.append(roc_auc_score(y_true_clean, y_scores_clean))
                    auprc_scores.append(average_precision_score(y_true_clean, y_scores_clean))
                    names.append(MODEL_LABELS[model])
                except Exception:
                    pass
    fig, ax = plt.subplots(1, 1, figsize=(14, 8))
    x = np.arange(len(names))
    width = 0.35
    bars1 = ax.bar(x - width/2, auroc_scores, width, label='AUROC', alpha=0.8, color='#2E86AB')
    bars2 = ax.bar(x + width/2, auprc_scores, width, label='AUPRC', alpha=0.8, color='#A23B72')
    ax.set_xlabel('Models', fontsize=14, fontweight='bold')
    ax.set_ylabel('Performance Score', fontsize=14, fontweight='bold')
    ax.set_title('Model Performance Comparison - PigGTEx Cross-Species Benchmark', fontsize=16, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=45, ha='right', fontsize=12)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_ylim(0, 1.05)
    for b in list(bars1) + list(bars2):
        h = b.get_height()
        ax.annotate(f'{h:.3f}', (b.get_x() + b.get_width()/2, h), xytext=(0, 3),
                    textcoords='offset points', ha='center', va='bottom', fontsize=10)
    plt.tight_layout()
    plt.savefig(f"{output_path}/pig_performance_comparison.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_path}/pig_performance_comparison.pdf", bbox_inches='tight')
    plt.close()


def plot_coverage(df: pd.DataFrame, output_path: str) -> None:
    # Exclude SpliceBERT, Nucleotide_Transformer, DNABERT2_Logistic (poor human performance)
    models = ['AlphaGenome', 'Evo2_zeroshot', 'Pangolin', 
              'SpliceAI', 'SpliceTransformer', 'MMSplice_pathogenicity']
    totals = len(df)
    data, names, colors = [], [], []
    for m in models:
        if m in df.columns:
            cov = df[m].notna().sum() / totals * 100
            data.append(cov)
            names.append(MODEL_LABELS[m])
            colors.append(MODEL_COLORS[m])
    fig, ax = plt.subplots(1, 1, figsize=(14, 6))
    bars = ax.barh(names, data, alpha=0.8, color=colors, edgecolor='black', linewidth=1.2)
    ax.set_xlabel('Coverage (%)', fontsize=14, fontweight='bold')
    ax.set_title('Model Coverage - PigGTEx Cross-Species Benchmark', fontsize=16, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='x')
    ax.set_xlim(0, 105)
    for bar, val in zip(bars, data):
        ax.text(val + 1, bar.get_y() + bar.get_height()/2, f'{val:.1f}%', va='center', ha='left', fontsize=11, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f"{output_path}/pig_coverage_analysis.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_path}/pig_coverage_analysis.pdf", bbox_inches='tight')
    plt.close()


def main():
    labels_path = "/Users/byronsun/Desktop/AS_复现模型/BIB_review/data/processed_data/pigGTEx/piggtex_silver_benchmark_balanced.tsv"
    output_path = "/Users/byronsun/Desktop/AS_复现模型/BIB_review/code/pigGTEx/baselines/visualisation"
    os.makedirs(output_path, exist_ok=True)

    print("📊 Building merged dataframe from individual model outputs (PigGTEx)...")
    df = build_merged_dataframe(labels_path)
    print(f"Loaded labels: {len(df)} variants; positive rate: {df['ground_truth'].mean()*100:.2f}%")

    # Save merged CSV similar to ratGTEx summary
    merged_csv = "/Users/byronsun/Desktop/AS_复现模型/BIB_review/code/pigGTEx/baselines/piggtex_all_model_predictions.csv"
    df.to_csv(merged_csv, index=False)
    print(f"💾 Saved merged predictions CSV: {merged_csv}")

    # Generate plots
    print("1) ROC curves...")
    plot_roc_curves(df, output_path)
    print("2) Precision-Recall curves...")
    plot_precision_recall_curves(df, output_path)
    print("3) Performance comparison bar chart...")
    plot_performance_bar_chart(df, output_path)
    print("4) Coverage analysis...")
    plot_coverage(df, output_path)
    print("✅ PigGTEx plots saved to:", output_path)


if __name__ == "__main__":
    main()


