#!/usr/bin/env python3
"""
Plot within-model strategy comparison for Evo2 and DNABERT-2
Compares zero-shot methods vs embedding+MLP approaches on both human and rat data

This figure is for the auxiliary analysis section, showing:
- Evo2: ΔlogP vs embedding+MLP
- DNABERT-2: Logistic probe vs embedding+MLP
- Performance on both SpliceVarDB (human) and ratGTEx (rat)
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (required to enable 3D projection)
import json

# Set publication-quality style
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 11
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['ytick.major.width'] = 1.5

def load_chicken_mlp_aurocs() -> tuple[float, float]:
    """Return (evo2_chicken_mlp_auroc, dnabert2_chicken_mlp_auroc) from provided JSON files."""
    evo2_path = Path('/Users/byronsun/Desktop/AS_复现模型/BIB_review/code/chickenGTEx/baselines/Evo2/results/evo2_mlp_chicken_train_test_metrics.json')
    db2_path = Path('/Users/byronsun/Desktop/AS_复现模型/BIB_review/code/chickenGTEx/baselines/DNABert_2/results/dnabert2_chicken_mlp_metrics.json')
    with evo2_path.open('r') as f:
        evo2_json = json.load(f)
    with db2_path.open('r') as f:
        db2_json = json.load(f)
    return float(evo2_json['auroc']), float(db2_json['auroc'])


# Assemble data (keep human/rat/pig numbers; add chicken from JSON and zero-shot baselines)
evo2_ch_mlp, dnabert2_ch_mlp = load_chicken_mlp_aurocs()

data = {
    'Model': ['Evo2', 'Evo2', 'DNABERT-2', 'DNABERT-2'],
    'Method': ['Zero-shot\n(ΔlogP)', 'Embedding\n+MLP', 'Zero-shot\n(Logistic)', 'Embedding\n+MLP'],
    'Human_AUROC': [0.7089, 0.8790, 0.5000, 0.9603],
    'Rat_AUROC':   [0.4848, 0.7761, 0.5000, 0.7112],
    'Pig_AUROC':   [0.5042, 0.6851, 0.5000, 0.7008],
    'Chicken_AUROC': [0.4934, evo2_ch_mlp, 0.5000, dnabert2_ch_mlp],
}

df = pd.DataFrame(data)

# Color scheme
colors = {
    'Evo2_zeroshot': '#ff7f0e',  # Orange
    'Evo2_mlp': '#ff4500',       # Darker orange
    'DNABERT2_zeroshot': '#17becf',  # Cyan
    'DNABERT2_mlp': '#008080'    # Darker cyan
}


def plot_grouped_bar_comparison(df, output_path):
    """
    Grouped bar chart comparing zero-shot vs embedding+MLP for both models across 3 species
    """
    fig, axes = plt.subplots(1, 4, figsize=(26, 6))
    
    # Human SpliceVarDB
    ax = axes[0]
    x = np.arange(2)  # Evo2, DNABERT-2
    width = 0.35
    
    human_zeroshot = [df.loc[0, 'Human_AUROC'], df.loc[2, 'Human_AUROC']]
    human_mlp = [df.loc[1, 'Human_AUROC'], df.loc[3, 'Human_AUROC']]
    
    bars1 = ax.bar(x - width/2, human_zeroshot, width, 
                   label='Zero-shot', color=['#ff7f0e', '#17becf'],
                   edgecolor='black', linewidth=1.5, alpha=0.8)
    bars2 = ax.bar(x + width/2, human_mlp, width,
                   label='Embedding+MLP', color=['#ff4500', '#008080'],
                   edgecolor='black', linewidth=1.5, alpha=0.8)
    
    # Add value labels on bars
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{height:.3f}',
                   ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    # Add improvement arrows and percentages
    for i, (z, m) in enumerate(zip(human_zeroshot, human_mlp)):
        improvement = ((m - z) / z) * 100
        ax.annotate('', xy=(i + width/2, m - 0.02), xytext=(i - width/2, z + 0.02),
                   arrowprops=dict(arrowstyle='->', lw=2, color='green', alpha=0.6))
        ax.text(i, (z + m) / 2, f'+{improvement:.1f}%',
               ha='center', va='center', fontsize=9, 
               bbox=dict(boxstyle='round,pad=0.3', facecolor='lightgreen', alpha=0.7))
    
    ax.set_ylabel('AUROC', fontsize=13, fontweight='bold')
    ax.set_title('Human (SpliceVarDB)', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(['Evo2', 'DNABERT-2'], fontsize=12, fontweight='bold')
    ax.set_ylim(0, 1.05)
    ax.legend(loc='upper left', fontsize=11, frameon=True, fancybox=True, shadow=True)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.axhline(y=0.5, color='red', linestyle='--', linewidth=1, alpha=0.5, label='Random')
    
    # Rat ratGTEx
    ax = axes[1]
    
    rat_zeroshot = [df.loc[0, 'Rat_AUROC'], df.loc[2, 'Rat_AUROC']]
    rat_mlp = [df.loc[1, 'Rat_AUROC'], df.loc[3, 'Rat_AUROC']]
    
    bars1 = ax.bar(x - width/2, rat_zeroshot, width,
                   label='Zero-shot', color=['#ff7f0e', '#17becf'],
                   edgecolor='black', linewidth=1.5, alpha=0.8)
    bars2 = ax.bar(x + width/2, rat_mlp, width,
                   label='Embedding+MLP', color=['#ff4500', '#008080'],
                   edgecolor='black', linewidth=1.5, alpha=0.8)
    
    # Add value labels on bars
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{height:.3f}',
                   ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    # Add improvement arrows and percentages
    for i, (z, m) in enumerate(zip(rat_zeroshot, rat_mlp)):
        improvement = ((m - z) / z) * 100
        ax.annotate('', xy=(i + width/2, m - 0.02), xytext=(i - width/2, z + 0.02),
                   arrowprops=dict(arrowstyle='->', lw=2, color='green', alpha=0.6))
        ax.text(i, (z + m) / 2, f'+{improvement:.1f}%',
               ha='center', va='center', fontsize=9,
               bbox=dict(boxstyle='round,pad=0.3', facecolor='lightgreen', alpha=0.7))
    
    ax.set_ylabel('AUROC', fontsize=13, fontweight='bold')
    ax.set_title('Rat (ratGTEx)', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(['Evo2', 'DNABERT-2'], fontsize=12, fontweight='bold')
    ax.set_ylim(0, 1.05)
    ax.legend(loc='upper left', fontsize=11, frameon=True, fancybox=True, shadow=True)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.axhline(y=0.5, color='red', linestyle='--', linewidth=1, alpha=0.5, label='Random')
    
    # Pig PigGTEx
    ax = axes[2]
    
    pig_zeroshot = [df.loc[0, 'Pig_AUROC'], df.loc[2, 'Pig_AUROC']]
    pig_mlp = [df.loc[1, 'Pig_AUROC'], df.loc[3, 'Pig_AUROC']]
    
    bars1 = ax.bar(x - width/2, pig_zeroshot, width,
                   label='Zero-shot', color=['#ff7f0e', '#17becf'],
                   edgecolor='black', linewidth=1.5, alpha=0.8)
    bars2 = ax.bar(x + width/2, pig_mlp, width,
                   label='Embedding+MLP', color=['#ff4500', '#008080'],
                   edgecolor='black', linewidth=1.5, alpha=0.8)
    
    # Add value labels on bars
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{height:.3f}',
                   ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    # Add improvement arrows and percentages
    for i, (z, m) in enumerate(zip(pig_zeroshot, pig_mlp)):
        improvement = ((m - z) / z) * 100
        ax.annotate('', xy=(i + width/2, m - 0.02), xytext=(i - width/2, z + 0.02),
                   arrowprops=dict(arrowstyle='->', lw=2, color='green', alpha=0.6))
        ax.text(i, (z + m) / 2, f'+{improvement:.1f}%',
               ha='center', va='center', fontsize=9,
               bbox=dict(boxstyle='round,pad=0.3', facecolor='lightgreen', alpha=0.7))
    
    ax.set_ylabel('AUROC', fontsize=13, fontweight='bold')
    ax.set_title('Pig (PigGTEx)', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(['Evo2', 'DNABERT-2'], fontsize=12, fontweight='bold')
    ax.set_ylim(0, 1.05)
    ax.legend(loc='upper left', fontsize=11, frameon=True, fancybox=True, shadow=True)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.axhline(y=0.5, color='red', linestyle='--', linewidth=1, alpha=0.5, label='Random')

    # Chicken ChickenGTEx
    ax = axes[3]
    chicken_zeroshot = [df.loc[0, 'Chicken_AUROC'], df.loc[2, 'Chicken_AUROC']]
    chicken_mlp = [df.loc[1, 'Chicken_AUROC'], df.loc[3, 'Chicken_AUROC']]

    bars1 = ax.bar(x - width/2, chicken_zeroshot, width,
                   label='Zero-shot', color=['#ff7f0e', '#17becf'],
                   edgecolor='black', linewidth=1.5, alpha=0.8)
    bars2 = ax.bar(x + width/2, chicken_mlp, width,
                   label='Embedding+MLP', color=['#ff4500', '#008080'],
                   edgecolor='black', linewidth=1.5, alpha=0.8)

    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                    f'{height:.3f}',
                    ha='center', va='bottom', fontsize=10, fontweight='bold')

    for i, (z, m) in enumerate(zip(chicken_zeroshot, chicken_mlp)):
        improvement = ((m - z) / z) * 100 if z != 0 else 0.0
        ax.annotate('', xy=(i + width/2, m - 0.02), xytext=(i - width/2, z + 0.02),
                    arrowprops=dict(arrowstyle='->', lw=2, color='green', alpha=0.6))
        ax.text(i, (z + m) / 2, f'+{improvement:.1f}%',
                ha='center', va='center', fontsize=9,
                bbox=dict(boxstyle='round,pad=0.3', facecolor='lightgreen', alpha=0.7))

    ax.set_ylabel('AUROC', fontsize=13, fontweight='bold')
    ax.set_title('Chicken (ChickenGTEx)', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(['Evo2', 'DNABERT-2'], fontsize=12, fontweight='bold')
    ax.set_ylim(0, 1.05)
    ax.legend(loc='upper left', fontsize=11, frameon=True, fancybox=True, shadow=True)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.axhline(y=0.5, color='red', linestyle='--', linewidth=1, alpha=0.5, label='Random')
    
    plt.tight_layout()
    plt.savefig(f"{output_path}/within_model_strategy_comparison_grouped.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_path}/within_model_strategy_comparison_grouped.pdf", bbox_inches='tight')
    plt.close()


def plot_cross_species_degradation(df, output_path):
    """
    Line plot showing cross-species performance degradation for each strategy
    """
    fig, ax = plt.subplots(1, 1, figsize=(12, 7))
    
    species = ['Human', 'Rat', 'Pig', 'Chicken']
    
    # Evo2
    evo2_zs_values = [df.loc[0, 'Human_AUROC'], df.loc[0, 'Rat_AUROC'], df.loc[0, 'Pig_AUROC'], df.loc[0, 'Chicken_AUROC']]
    evo2_mlp_values = [df.loc[1, 'Human_AUROC'], df.loc[1, 'Rat_AUROC'], df.loc[1, 'Pig_AUROC'], df.loc[1, 'Chicken_AUROC']]
    
    ax.plot(species, evo2_zs_values,
           marker='o', markersize=12, linewidth=3, color='#ff7f0e', 
           label='Evo2 (Zero-shot ΔlogP)', linestyle='-')
    ax.plot(species, evo2_mlp_values,
           marker='s', markersize=12, linewidth=3, color='#ff4500',
           label='Evo2 (Embedding+MLP)', linestyle='--')
    
    # DNABERT-2
    db2_zs_values = [df.loc[2, 'Human_AUROC'], df.loc[2, 'Rat_AUROC'], df.loc[2, 'Pig_AUROC'], df.loc[2, 'Chicken_AUROC']]
    db2_mlp_values = [df.loc[3, 'Human_AUROC'], df.loc[3, 'Rat_AUROC'], df.loc[3, 'Pig_AUROC'], df.loc[3, 'Chicken_AUROC']]
    
    ax.plot(species, db2_zs_values,
           marker='o', markersize=12, linewidth=3, color='#17becf',
           label='DNABERT-2 (Zero-shot Logistic)', linestyle='-')
    ax.plot(species, db2_mlp_values,
           marker='s', markersize=12, linewidth=3, color='#008080',
           label='DNABERT-2 (Embedding+MLP)', linestyle='--')
    
    # Add value labels
    for i, sp in enumerate(species):
        for row_idx, values in enumerate([evo2_zs_values, evo2_mlp_values, db2_zs_values, db2_mlp_values]):
            val = values[i]
            ax.text(i, val, f"{val:.3f}", 
                   ha='center', va='bottom', fontsize=8, fontweight='bold')
    
    ax.axhline(y=0.5, color='red', linestyle='--', linewidth=2, alpha=0.5, label='Random')
    ax.set_ylabel('AUROC', fontsize=14, fontweight='bold')
    ax.set_xlabel('Species', fontsize=14, fontweight='bold')
    ax.set_title('Cross-Species Performance: Zero-shot vs Embedding+MLP', 
                fontsize=15, fontweight='bold')
    ax.set_ylim(0.4, 1.0)
    ax.legend(loc='upper right', fontsize=11, frameon=True, fancybox=True, shadow=True)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    plt.savefig(f"{output_path}/within_model_cross_species_degradation.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_path}/within_model_cross_species_degradation.pdf", bbox_inches='tight')
    plt.close()


def plot_improvement_heatmap(df, output_path):
    """
    Heatmap showing improvement from zero-shot to embedding+MLP
    """
    # Calculate improvements
    improvements = pd.DataFrame({
        'Human': [
            df.loc[1, 'Human_AUROC'] - df.loc[0, 'Human_AUROC'],  # Evo2
            df.loc[3, 'Human_AUROC'] - df.loc[2, 'Human_AUROC']   # DNABERT-2
        ],
        'Rat': [
            df.loc[1, 'Rat_AUROC'] - df.loc[0, 'Rat_AUROC'],      # Evo2
            df.loc[3, 'Rat_AUROC'] - df.loc[2, 'Rat_AUROC']       # DNABERT-2
        ],
        'Pig': [
            df.loc[1, 'Pig_AUROC'] - df.loc[0, 'Pig_AUROC'],      # Evo2
            df.loc[3, 'Pig_AUROC'] - df.loc[2, 'Pig_AUROC']       # DNABERT-2
        ],
        'Chicken': [
            df.loc[1, 'Chicken_AUROC'] - df.loc[0, 'Chicken_AUROC'],  # Evo2
            df.loc[3, 'Chicken_AUROC'] - df.loc[2, 'Chicken_AUROC']   # DNABERT-2
        ]
    }, index=['Evo2', 'DNABERT-2'])
    
    fig, ax = plt.subplots(1, 1, figsize=(10, 5))
    
    sns.heatmap(improvements, annot=True, fmt='.3f', cmap='RdYlGn', center=0,
               cbar_kws={'label': 'AUROC Improvement\n(Embedding+MLP - Zero-shot)'},
               linewidths=2, linecolor='black', ax=ax,
               vmin=-0.1, vmax=0.5, annot_kws={'fontsize': 14, 'fontweight': 'bold'})
    
    ax.set_xlabel('Species', fontsize=13, fontweight='bold')
    ax.set_ylabel('Foundation Model', fontsize=13, fontweight='bold')
    ax.set_title('Performance Gain from Supervised Embedding Interpretation', 
                fontsize=14, fontweight='bold', pad=20)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    
    plt.tight_layout()
    plt.savefig(f"{output_path}/within_model_improvement_heatmap.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_path}/within_model_improvement_heatmap.pdf", bbox_inches='tight')
    plt.close()


def plot_combined_figure(df, output_path):
    """
    Combined figure with 3 panels: (A) Human comparison, (B) Rat comparison, (C) Cross-species
    """
    fig = plt.figure(figsize=(18, 6))
    gs = fig.add_gridspec(1, 3, hspace=0.3, wspace=0.3)
    
    # Panel A: Human SpliceVarDB
    ax1 = fig.add_subplot(gs[0, 0])
    x = np.arange(2)
    width = 0.35
    
    human_zeroshot = [df.loc[0, 'Human_AUROC'], df.loc[2, 'Human_AUROC']]
    human_mlp = [df.loc[1, 'Human_AUROC'], df.loc[3, 'Human_AUROC']]
    
    bars1 = ax1.bar(x - width/2, human_zeroshot, width,
                    label='Zero-shot', color=['#ff7f0e', '#17becf'],
                    edgecolor='black', linewidth=1.5, alpha=0.8)
    bars2 = ax1.bar(x + width/2, human_mlp, width,
                    label='Embedding+MLP', color=['#ff4500', '#008080'],
                    edgecolor='black', linewidth=1.5, alpha=0.8)
    
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height,
                    f'{height:.3f}', ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    ax1.set_ylabel('AUROC', fontsize=12, fontweight='bold')
    ax1.set_title('(A) Human (SpliceVarDB)', fontsize=13, fontweight='bold', loc='left')
    ax1.set_xticks(x)
    ax1.set_xticklabels(['Evo2', 'DNABERT-2'], fontsize=11)
    ax1.set_ylim(0, 1.05)
    ax1.legend(loc='upper left', fontsize=10, frameon=True)
    ax1.grid(axis='y', alpha=0.3, linestyle='--')
    ax1.axhline(y=0.5, color='red', linestyle='--', linewidth=1, alpha=0.5)
    
    # Panel B: Rat ratGTEx
    ax2 = fig.add_subplot(gs[0, 1])
    
    rat_zeroshot = [df.loc[0, 'Rat_AUROC'], df.loc[2, 'Rat_AUROC']]
    rat_mlp = [df.loc[1, 'Rat_AUROC'], df.loc[3, 'Rat_AUROC']]
    
    bars1 = ax2.bar(x - width/2, rat_zeroshot, width,
                    label='Zero-shot', color=['#ff7f0e', '#17becf'],
                    edgecolor='black', linewidth=1.5, alpha=0.8)
    bars2 = ax2.bar(x + width/2, rat_mlp, width,
                    label='Embedding+MLP', color=['#ff4500', '#008080'],
                    edgecolor='black', linewidth=1.5, alpha=0.8)
    
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height,
                    f'{height:.3f}', ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    ax2.set_ylabel('AUROC', fontsize=12, fontweight='bold')
    ax2.set_title('(B) Rat (ratGTEx)', fontsize=13, fontweight='bold', loc='left')
    ax2.set_xticks(x)
    ax2.set_xticklabels(['Evo2', 'DNABERT-2'], fontsize=11)
    ax2.set_ylim(0, 1.05)
    ax2.legend(loc='upper left', fontsize=10, frameon=True)
    ax2.grid(axis='y', alpha=0.3, linestyle='--')
    ax2.axhline(y=0.5, color='red', linestyle='--', linewidth=1, alpha=0.5)
    
    # Panel C: Cross-species comparison
    ax3 = fig.add_subplot(gs[0, 2])
    
    species = ['Human', 'Rat', 'Pig', 'Chicken']
    
    ax3.plot(species, [df.loc[0, 'Human_AUROC'], df.loc[0, 'Rat_AUROC'], df.loc[0, 'Pig_AUROC'], df.loc[0, 'Chicken_AUROC']],
            marker='o', markersize=10, linewidth=2.5, color='#ff7f0e',
            label='Evo2 (Zero-shot)', linestyle='-')
    ax3.plot(species, [df.loc[1, 'Human_AUROC'], df.loc[1, 'Rat_AUROC'], df.loc[1, 'Pig_AUROC'], df.loc[1, 'Chicken_AUROC']],
            marker='s', markersize=10, linewidth=2.5, color='#ff4500',
            label='Evo2 (Emb+MLP)', linestyle='--')
    ax3.plot(species, [df.loc[2, 'Human_AUROC'], df.loc[2, 'Rat_AUROC'], df.loc[2, 'Pig_AUROC'], df.loc[2, 'Chicken_AUROC']],
            marker='o', markersize=10, linewidth=2.5, color='#17becf',
            label='DNABERT-2 (Zero-shot)', linestyle='-')
    ax3.plot(species, [df.loc[3, 'Human_AUROC'], df.loc[3, 'Rat_AUROC'], df.loc[3, 'Pig_AUROC'], df.loc[3, 'Chicken_AUROC']],
            marker='s', markersize=10, linewidth=2.5, color='#008080',
            label='DNABERT-2 (Emb+MLP)', linestyle='--')
    
    ax3.axhline(y=0.5, color='red', linestyle='--', linewidth=1.5, alpha=0.5, label='Random')
    ax3.set_ylabel('AUROC', fontsize=12, fontweight='bold')
    ax3.set_title('(C) Cross-Species Performance', fontsize=13, fontweight='bold', loc='left')
    ax3.set_ylim(0.4, 1.0)
    ax3.legend(loc='upper right', fontsize=9, frameon=True)
    ax3.grid(axis='y', alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    plt.savefig(f"{output_path}/within_model_combined_figure.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_path}/within_model_combined_figure.pdf", bbox_inches='tight')
    plt.close()


def plot_two_layer_strategy_figure(df, output_path):
    """
    Single figure, two visual layers per species:
    - Front layer: Evo2 (two bars per species: Zero-shot, Embedding+MLP)
    - Back layer:  DNABERT-2 (two bars per species: Zero-shot, Embedding+MLP)
    Legend encodes strategy only, using shared colours across models.
    """
    species = ['Human', 'Rat', 'Pig', 'Chicken']

    # Values per species
    evo2_zs = [df.loc[0, f'{sp}_AUROC'] for sp in species]
    evo2_mlp = [df.loc[1, f'{sp}_AUROC'] for sp in species]
    db2_zs  = [df.loc[2, f'{sp}_AUROC'] for sp in species]
    db2_mlp = [df.loc[3, f'{sp}_AUROC'] for sp in species]

    x = np.arange(len(species), dtype=float)
    width_back = 0.24      # slimmer back layer
    width_front = 0.16     # slimmer front layer
    delta = 0.10           # smaller separation within each species

    # Colours shared by both models for strategies
    color_zs = '#9E9E9E'    # grey
    color_mlp = '#1F77B4'   # blue

    fig, ax = plt.subplots(1, 1, figsize=(12.5, 6.0))

    # Back layer: DNABERT-2
    ax.bar(x - delta, db2_zs, width=width_back, color=color_zs, alpha=0.35,
           edgecolor='black', linewidth=1.2, label='Zero-shot', zorder=1)
    ax.bar(x + delta, db2_mlp, width=width_back, color=color_mlp, alpha=0.35,
           edgecolor='black', linewidth=1.2, label='Embedding+MLP', zorder=1)

    # Front layer: Evo2 (narrower so背层可见)
    bars_front_zs = ax.bar(x - delta, evo2_zs, width=width_front, color=color_zs,
                           edgecolor='black', linewidth=1.2, zorder=2)
    bars_front_mlp = ax.bar(x + delta, evo2_mlp, width=width_front, color=color_mlp,
                            edgecolor='black', linewidth=1.2, zorder=2)

    # Value labels (front layer)
    for bars in [bars_front_zs, bars_front_mlp]:
        for bar in bars:
            h = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2, h + 0.01, f"{h:.3f}",
                    ha='center', va='bottom', fontsize=9, fontweight='bold', zorder=3)

    # Axes formatting
    ax.set_xticks(x)
    ax.set_xticklabels(species, fontsize=12, fontweight='bold')
    ax.set_ylabel('AUROC', fontsize=13, fontweight='bold')
    ax.set_title('Within-Model Strategy Comparison across Species (Front: Evo2; Back: DNABERT-2)',
                 fontsize=14, fontweight='bold')
    ax.set_ylim(0.4, 1.05)
    ax.grid(axis='y', alpha=0.3, linestyle='--', zorder=0)
    ax.axhline(y=0.5, color='red', linestyle='--', linewidth=1.2, alpha=0.5)

    # Legend: strategies only
    handles = [
        Patch(facecolor=color_zs, edgecolor='black', label='Zero-shot'),
        Patch(facecolor=color_mlp, edgecolor='black', label='Embedding+MLP')
    ]
    ax.legend(handles=handles, loc='upper right', fontsize=11, frameon=True, fancybox=True,
              title='Strategy')

    # Sub-label for layers to indicate models
    ax.text(0.01, 0.97, 'Front layer: Evo2', transform=ax.transAxes, fontsize=11,
            fontweight='bold', ha='left', va='center')
    ax.text(0.01, 0.93, 'Back layer: DNABERT-2', transform=ax.transAxes, fontsize=11,
            ha='left', va='center')

    plt.tight_layout()
    plt.savefig(f"{output_path}/within_model_two_layer_comparison.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_path}/within_model_two_layer_comparison.pdf", bbox_inches='tight')
    plt.close()


def plot_two_layer_strategy_3d(df, output_path):
    """
    3D perspective version of the two-layer comparison requested:
    - X-axis: Species (Human, Rat, Pig, Chicken)
    - Y-axis (depth): Model layer (front=Evo2 at y=0, back=DNABERT-2 at y=1)
    - For each species and layer: two bars (Zero-shot, Embedding+MLP)
    Legend encodes strategy only; colouring shared between models.
    """
    species = ['Human', 'Rat', 'Pig', 'Chicken']

    evo2_zs = [df.loc[0, f'{sp}_AUROC'] for sp in species]
    evo2_mlp = [df.loc[1, f'{sp}_AUROC'] for sp in species]
    db2_zs  = [df.loc[2, f'{sp}_AUROC'] for sp in species]
    db2_mlp = [df.loc[3, f'{sp}_AUROC'] for sp in species]

    # Common colours for strategies
    color_zs = '#9E9E9E'
    color_mlp = '#6BAED6'  # lighter blue for better visibility

    fig = plt.figure(figsize=(12.5, 7.0))
    ax = fig.add_subplot(111, projection='3d')

    # Bar geometry
    # Geometry (front layer a bit narrower to expose the back layer)
    dx_front = 0.18              # front layer bar width
    dx_back  = 0.24              # back layer bar width
    dy = 0.10                    # thinner depth (front-to-back direction)
    x = np.arange(len(species), dtype=float)  # species positions on X axis (tick centers)
    sep = 0.12                   # smaller separation around species center

    # Layer Y positions (only two): front (Evo2) and back (DNABERT-2)
    y_front = 0.0
    y_back = 1.0

    # Back layer: DNABERT-2 — two strategies centered at (x -/+ sep)
    ax.bar3d(x - sep - dx_back/2, y_back - dy/2, np.zeros_like(x), dx_back, dy, db2_zs,
             color=color_zs, alpha=1.0, edgecolor='black', linewidth=1.0, zorder=1)
    ax.bar3d(x + sep - dx_back/2, y_back - dy/2, np.zeros_like(x), dx_back, dy, db2_mlp,
             color=color_mlp, alpha=1.0, edgecolor='black', linewidth=1.0, zorder=1)

    # Front layer: Evo2 — two strategies centered at (x -/+ sep)
    bars_front_zs = ax.bar3d(x - sep - dx_front/2, y_front - dy/2, np.zeros_like(x), dx_front, dy, evo2_zs,
                             color=color_zs, alpha=1.0, edgecolor='black', linewidth=1.0, zorder=2)
    bars_front_mlp = ax.bar3d(x + sep - dx_front/2, y_front - dy/2, np.zeros_like(x), dx_front, dy, evo2_mlp,
                              color=color_mlp, alpha=1.0, edgecolor='black', linewidth=1.0, zorder=2)

    # Axes labels (no label for X and Y, only Z)
    ax.set_zlabel('AUROC', labelpad=12, fontsize=12, fontweight='bold')
    ax.set_title('Two-Layer Strategy Comparison (3D) — Front: Evo2, Back: DNABERT-2',
                 fontsize=14, fontweight='bold', pad=18)

    # Ticks
    ax.set_xticks(x)
    ax.set_xticklabels(species, fontsize=11, fontweight='bold')
    ax.set_yticks([y_front, y_back])
    ax.set_yticklabels(['Evo2', 'DNABERT-2'], fontsize=11, fontweight='bold')
    ax.set_zlim(0.0, 1.05)

    # Set axis limits to align ticks exactly under grouped bars
    ax.set_xlim(-0.6, len(species)-0.4)
    ax.set_ylim(-0.2, 1.2)

    # Ground plane at Z=0 (red reference plane at bottom)
    xx, yy = np.meshgrid(np.linspace(-0.6, len(species)-0.4, 2), np.linspace(-0.2, 1.2, 2))
    zz = np.full_like(xx, 0.0)
    ax.plot_surface(xx, yy, zz, alpha=0.20, color='#FF6B6B', zorder=0)

    # Camera/view for perspective similar to the reference image
    ax.view_init(elev=23, azim=-60)

    # Legend: strategy only
    handles = [
        Patch(facecolor=color_zs, edgecolor='black', label='Zero-shot'),
        Patch(facecolor=color_mlp, edgecolor='black', label='Embedding+MLP')
    ]
    ax.legend(handles=handles, loc='upper right', bbox_to_anchor=(1.15, 1.0),
              fontsize=11, frameon=True, fancybox=True, title='Strategy')

    plt.tight_layout()
    plt.savefig(f"{output_path}/within_model_two_layer_comparison_3d.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_path}/within_model_two_layer_comparison_3d.pdf", bbox_inches='tight')
    plt.savefig(f"{output_path}/within_model_two_layer_comparison_3d.svg", format='svg', bbox_inches='tight')
    plt.close()

def main():
    output_dir = Path("/Users/byronsun/Desktop/AS_复现模型/BIB_review/code/splicevardb/baselines/visualisation")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("Generating within-model strategy comparison figures...")
    
    # 1. Grouped bar comparison (Human and Rat side by side)
    print("  [1/4] Grouped bar comparison...")
    plot_grouped_bar_comparison(df, output_dir)
    
    # 2. Cross-species degradation line plot
    print("  [2/4] Cross-species degradation plot...")
    plot_cross_species_degradation(df, output_dir)
    
    # 3. Improvement heatmap
    print("  [3/4] Improvement heatmap...")
    plot_improvement_heatmap(df, output_dir)
    
    # 4. Combined figure (recommended for publication)
    print("  [4/4] Combined figure (3-panel)...")
    plot_combined_figure(df, output_dir)
    
    # 5. Two-layer single figure (requested)
    print("  [+] Two-layer single figure...")
    plot_two_layer_strategy_figure(df, output_dir)
    
    # 6. Two-layer 3D figure (requested perspective)
    print("  [+] Two-layer 3D figure...")
    plot_two_layer_strategy_3d(df, output_dir)
    
    print(f"\n✓ All figures saved to: {output_dir}")
    print("\nGenerated files:")
    print("  - within_model_strategy_comparison_grouped.png/pdf")
    print("  - within_model_cross_species_degradation.png/pdf")
    print("  - within_model_improvement_heatmap.png/pdf")
    print("  - within_model_combined_figure.png/pdf (RECOMMENDED for paper)")
    print("  - within_model_two_layer_comparison.png/pdf")
    print("  - within_model_two_layer_comparison_3d.png/pdf")


if __name__ == "__main__":
    main()

