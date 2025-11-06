#!/usr/bin/env python3
"""
Generate three publication-ready visualisations from SpliceVarDB predictions:
1) Score distributions and separability: KDE/violin-style small multiples + optional waterfall for positives.
2) Score correlation heatmap (Spearman) across models to reveal redundancy/complementarity.
3) Bootstrap 95% CI forest plot for AUROC and AUPRC (resampling variants; no model re-runs).

Inputs
- splicevardb_all_model_predictions.csv (created by summary_results.py)

Outputs
- visualisation/distributions_ridge.{png,pdf}
- visualisation/waterfall_selected.{png,pdf}
- visualisation/correlation_heatmap.{png,pdf}
- visualisation/forest_bootstrap_ci.{png,pdf}

Notes
- English comments as required. Plots styled for journal use.
"""

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_auc_score, average_precision_score


DEFAULT_MODELS = [
    'AlphaGenome',
    'Evo2',
    'Nucleotide_Transformer',
    'Pangolin',
    'SpliceAI',
    'SpliceBERT',
    'SpliceTransformer',
    'DNABERT2_Logistic',
    'MMSplice_pathogenicity',
]


def load_data(csv_path: str, models: list[str]) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    # Ensure binary labels available
    if 'ground_truth_binary' not in df.columns and 'ground_truth' in df.columns:
        df['ground_truth_binary'] = (df['ground_truth'] == 'Splice-altering').astype(int)
    # Keep only required columns
    keep_cols = ['variant_id', 'ground_truth_binary'] + [m for m in models if m in df.columns]
    return df[keep_cols].copy()


def plot_score_distributions(df: pd.DataFrame, models: list[str], out_dir: str):
    """Small multiples: per-model KDE for positive vs negative; also violin overlay."""
    sns.set_style('whitegrid')
    y = df['ground_truth_binary'].values

    n = len(models)
    ncols = 3
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 5, nrows * 3.3), sharex=False)
    axes = axes.flatten()

    for i, model in enumerate(models):
        ax = axes[i]
        if model not in df.columns:
            ax.axis('off')
            continue
        scores = df[model].values
        mask = ~np.isnan(scores)
        s = scores[mask]
        yt = y[mask]
        if s.size == 0:
            ax.axis('off')
            continue
        # Two KDEs: positives and negatives
        sns.kdeplot(x=s[yt == 1], ax=ax, color='#d62728', label='Splice-altering', linewidth=2,
                    warn_singular=False, common_norm=False)
        sns.kdeplot(x=s[yt == 0], ax=ax, color='#1f77b4', label='Normal', linewidth=2,
                    warn_singular=False, common_norm=False)
        ax.set_title(model.replace('_', ' '), fontsize=11, fontweight='bold')
        ax.set_xlabel('Score', fontsize=10)
        ax.set_ylabel('Density', fontsize=10)
        # Show AUROC/AUPRC in subtitle if feasible
        try:
            auroc = roc_auc_score(yt, s)
            auprc = average_precision_score(yt, s)
            ax.text(0.98, 0.95, f"AUROC {auroc:.3f}\nAUPRC {auprc:.3f}", transform=ax.transAxes,
                    ha='right', va='top', fontsize=9,
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.6, edgecolor='0.7'))
        except Exception:
            pass
        if i == 0:
            ax.legend(frameon=True, fontsize=9)

    # Hide unused axes
    for j in range(i + 1, len(axes)):
        axes[j].axis('off')

    plt.tight_layout()
    os.makedirs(out_dir, exist_ok=True)
    plt.savefig(os.path.join(out_dir, 'distributions_ridge.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(out_dir, 'distributions_ridge.pdf'), bbox_inches='tight')
    plt.close(fig)


def plot_waterfall_for_selected(df: pd.DataFrame, selected: list[str], out_dir: str, top_n: int = 1000):
    """Waterfall bars for positives sorted descending to illustrate thresholding difficulty."""
    pos_df = df[df['ground_truth_binary'] == 1]
    n = len(selected)
    fig, axes = plt.subplots(1, n, figsize=(n * 5, 3.5), sharey=True)
    if n == 1:
        axes = [axes]
    for ax, model in zip(axes, selected):
        if model not in pos_df.columns:
            ax.axis('off')
            continue
        s = pos_df[model].dropna().values
        s_sorted = np.sort(s)[::-1]
        if top_n is not None:
            s_sorted = s_sorted[:top_n]
        ax.bar(np.arange(len(s_sorted)), s_sorted, width=1.0, color='#8c564b', alpha=0.9)
        ax.set_title(f'{model.replace("_", " ")}: positive variants', fontsize=11)
        ax.set_xlabel('Rank (desc)')
        ax.set_ylabel('Score')
    plt.tight_layout()
    os.makedirs(out_dir, exist_ok=True)
    plt.savefig(os.path.join(out_dir, 'waterfall_selected.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(out_dir, 'waterfall_selected.pdf'), bbox_inches='tight')
    plt.close(fig)


def plot_correlation_heatmap(df: pd.DataFrame, models: list[str], out_dir: str):
    """Spearman correlation across model scores (pairwise, complete observations)."""
    scores_df = df[[m for m in models if m in df.columns]].copy()
    corr = scores_df.corr(method='spearman', min_periods=100)
    plt.figure(figsize=(8, 6))
    sns.heatmap(corr, cmap='vlag', vmin=-1, vmax=1, annot=True, fmt='.2f',
                cbar_kws={'label': 'Spearman'})
    plt.title('Model score correlation (Spearman)')
    plt.tight_layout()
    os.makedirs(out_dir, exist_ok=True)
    plt.savefig(os.path.join(out_dir, 'correlation_heatmap.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(out_dir, 'correlation_heatmap.pdf'), bbox_inches='tight')
    plt.close()


def bootstrap_ci(y_true: np.ndarray, y_score: np.ndarray, n_boot: int = 200, rng=None):
    """Return point estimates and 95% CI for AUROC and AUPRC using bootstrap on indices."""
    if rng is None:
        rng = np.random.RandomState(42)
    # Keep only finite
    mask = np.isfinite(y_score)
    y = y_true[mask]
    s = y_score[mask]
    if y.size == 0:
        return np.nan, (np.nan, np.nan), np.nan, (np.nan, np.nan)
    try:
        auroc_full = roc_auc_score(y, s)
        auprc_full = average_precision_score(y, s)
    except Exception:
        auroc_full = np.nan
        auprc_full = np.nan
    auc_boot = []
    apr_boot = []
    n = len(y)
    for _ in range(n_boot):
        idx = rng.randint(0, n, size=n)
        try:
            auc_boot.append(roc_auc_score(y[idx], s[idx]))
            apr_boot.append(average_precision_score(y[idx], s[idx]))
        except Exception:
            continue
    if len(auc_boot) == 0:
        return auroc_full, (np.nan, np.nan), auprc_full, (np.nan, np.nan)
    lo_auroc, hi_auroc = np.percentile(auc_boot, [2.5, 97.5])
    lo_auprc, hi_auprc = np.percentile(apr_boot, [2.5, 97.5])
    return auroc_full, (lo_auroc, hi_auroc), auprc_full, (lo_auprc, hi_auprc)


def plot_forest_ci(df: pd.DataFrame, models: list[str], out_dir: str, n_boot: int = 200):
    """Forest plot with AUROC/AUPRC and their 95% CIs."""
    y = df['ground_truth_binary'].values
    rows = []
    for m in models:
        if m not in df.columns:
            continue
        s = df[m].values
        auroc, (lo_auc, hi_auc), auprc, (lo_ap, hi_ap) = bootstrap_ci(y, s, n_boot=n_boot)
        rows.append({
            'model': m,
            'auroc': auroc, 'auc_lo': lo_auc, 'auc_hi': hi_auc,
            'auprc': auprc, 'ap_lo': lo_ap, 'ap_hi': hi_ap,
        })
    res = pd.DataFrame(rows)
    # Sort by AUROC desc
    res = res.sort_values('auroc', ascending=False)

    # Forest for AUROC
    fig, ax = plt.subplots(figsize=(7, 0.5 * len(res) + 1.5))
    y_pos = np.arange(len(res))
    ax.errorbar(res['auroc'], y_pos, xerr=[res['auroc'] - res['auc_lo'], res['auc_hi'] - res['auroc']],
                fmt='o', color='#1f77b4', ecolor='#1f77b4', elinewidth=2, capsize=3)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(res['model'])
    ax.set_xlabel('AUROC (95% CI)')
    ax.set_title('AUROC with 95% bootstrap CI')
    ax.set_xlim(0.0, 1.0)
    plt.tight_layout()
    os.makedirs(out_dir, exist_ok=True)
    plt.savefig(os.path.join(out_dir, 'forest_bootstrap_ci_auroc.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(out_dir, 'forest_bootstrap_ci_auroc.pdf'), bbox_inches='tight')
    plt.close(fig)

    # Forest for AUPRC
    fig, ax = plt.subplots(figsize=(7, 0.5 * len(res) + 1.5))
    ax.errorbar(res['auprc'], y_pos, xerr=[res['auprc'] - res['ap_lo'], res['ap_hi'] - res['auprc']],
                fmt='o', color='#d62728', ecolor='#d62728', elinewidth=2, capsize=3)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(res['model'])
    ax.set_xlabel('AUPRC (95% CI)')
    ax.set_title('AUPRC with 95% bootstrap CI')
    ax.set_xlim(0.0, 1.0)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'forest_bootstrap_ci_auprc.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(out_dir, 'forest_bootstrap_ci_auprc.pdf'), bbox_inches='tight')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description='Distributions, correlation and bootstrap CI plots for SpliceVarDB predictions')
    parser.add_argument('--csv', default='/Users/byronsun/Desktop/AS_复现模型/BIB_review/code/splicevardb/baselines/splicevardb_all_model_predictions.csv',
                        help='Path to merged predictions CSV')
    parser.add_argument('--out', default='/Users/byronsun/Desktop/AS_复现模型/BIB_review/code/splicevardb/baselines/visualisation',
                        help='Output directory for figures')
    parser.add_argument('--n_boot', type=int, default=200, help='Number of bootstrap resamples for CI')
    parser.add_argument('--waterfall_models', nargs='*', default=['Evo2', 'Nucleotide_Transformer', 'SpliceBERT'],
                        help='Models to include in waterfall plot (positives only)')
    parser.add_argument('--waterfall_top', type=int, default=1000, help='Top-N positives to display in each waterfall')
    args = parser.parse_args()

    models = DEFAULT_MODELS
    df = load_data(args.csv, models)

    # 1) Distributions/separability
    plot_score_distributions(df, models, args.out)
    plot_waterfall_for_selected(df, args.waterfall_models, args.out, top_n=args.waterfall_top)

    # 2) Correlation heatmap
    plot_correlation_heatmap(df, models, args.out)

    # 3) Bootstrap CI forest
    plot_forest_ci(df, models, args.out, n_boot=args.n_boot)

    print('Done. Figures saved to:', args.out)


if __name__ == '__main__':
    main()


