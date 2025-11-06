#!/usr/bin/env python3
"""
Cross-Species Performance Comparison: Human vs Rat vs Pig vs Chicken
Generate publication-ready comparison figures for 4.2.3
Excludes NT, SpliceBERT, DNABERT2-Logistic (poor human performance)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_auc_score, average_precision_score
import warnings
import os
warnings.filterwarnings('ignore')

# Set style
plt.style.use('seaborn-v0_8-whitegrid')

# CONSISTENT COLORS across all plots
MODEL_COLORS = {
    'Pangolin': '#d62728',              # Red
    'SpliceAI': '#9467bd',              # Purple
    'MMSplice_pathogenicity': '#bcbd22', # Olive
    'SpliceTransformer': '#e377c2',     # Pink
    'AlphaGenome': '#1f77b4',           # Blue
    'Evo2_zeroshot': '#ff7f0e',         # Orange
}

MODEL_LABELS = {
    'Pangolin': 'Pangolin',
    'SpliceAI': 'SpliceAI',
    'MMSplice_pathogenicity': 'MMSplice',
    'SpliceTransformer': 'SpliceTransformer',
    'AlphaGenome': 'AlphaGenome',
    'Evo2_zeroshot': 'Evo2',
}

# Model order: Task-specific first, then GFMs
MODEL_ORDER = ['Pangolin', 'SpliceAI', 'MMSplice_pathogenicity', 'SpliceTransformer',
               'AlphaGenome', 'Evo2_zeroshot']


def calculate_metrics(df, is_human=True):
    """Calculate AUROC and AUPRC for all models (excludes NT, SpliceBERT, DNABERT2-Logistic)"""
    # Use MODEL_ORDER for consistent ordering
    models = MODEL_ORDER
    
    if is_human:
        y_true = df['ground_truth_binary'].values
    else:
        y_true = df['ground_truth'].values
    
    results = {}
    
    for model in models:
        if model in df.columns:
            y_scores = df[model].values
            mask = ~np.isnan(y_scores)
            y_true_clean = y_true[mask]
            y_scores_clean = y_scores[mask]
            
            if len(y_true_clean) > 0 and len(np.unique(y_true_clean)) > 1:
                try:
                    auroc = roc_auc_score(y_true_clean, y_scores_clean)
                    auprc = average_precision_score(y_true_clean, y_scores_clean)
                    coverage = len(y_true_clean) / len(y_true) * 100
                    
                    results[model] = {
                        'auroc': auroc,
                        'auprc': auprc,
                        'coverage': coverage
                    }
                except Exception as e:
                    print(f"Warning: Could not calculate metrics for {model}: {e}")
    
    return results


def plot_side_by_side_comparison(human_results, rat_results, output_path):
    """Plot side-by-side AUROC comparison"""
    # Get common models in MODEL_ORDER
    common_models = [m for m in MODEL_ORDER if m in human_results and m in rat_results]
    
    model_names = [MODEL_LABELS[m].replace('\n', ' ') for m in common_models]
    human_auroc = [human_results[m]['auroc'] for m in common_models]
    rat_auroc = [rat_results[m]['auroc'] for m in common_models]
    
    # Create side-by-side bar chart
    fig, ax = plt.subplots(1, 1, figsize=(14, 8))
    
    x = np.arange(len(model_names))
    width = 0.35
    
    bars1 = ax.bar(x - width/2, human_auroc, width, label='Human (SpliceVarDB)', 
                   alpha=0.85, color='#2E86AB', edgecolor='black', linewidth=1.2)
    bars2 = ax.bar(x + width/2, rat_auroc, width, label='Rat (RatGTEx)', 
                   alpha=0.85, color='#F77F00', edgecolor='black', linewidth=1.2)
    
    ax.set_xlabel('Models', fontsize=14, fontweight='bold')
    ax.set_ylabel('AUROC', fontsize=14, fontweight='bold')
    ax.set_title('Cross-Species Performance Comparison: Human vs Rat', fontsize=16, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(model_names, rotation=45, ha='right', fontsize=11)
    ax.legend(fontsize=13, loc='upper right')
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_ylim(0, 1.05)
    
    # Add value labels
    for bar in bars1:
        height = bar.get_height()
        ax.annotate(f'{height:.3f}',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    for bar in bars2:
        height = bar.get_height()
        ax.annotate(f'{height:.3f}',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(f"{output_path}/cross_species_auroc_comparison.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_path}/cross_species_auroc_comparison.pdf", bbox_inches='tight')
    plt.close()
    print(f"  ✓ Saved cross-species AUROC comparison")


def plot_performance_drop(human_results, rat_results, output_path):
    """Plot performance drop from human to rat (models in task-specific then GFM order)"""
    # Get common models in MODEL_ORDER
    common_models = [m for m in MODEL_ORDER if m in human_results and m in rat_results]
    
    model_names = [MODEL_LABELS[m].replace('\n', ' ') for m in common_models]
    colors = [MODEL_COLORS[m] for m in common_models]
    
    # Calculate percentage drop
    auroc_drops = [(human_results[m]['auroc'] - rat_results[m]['auroc']) / human_results[m]['auroc'] * 100 
                   for m in common_models]
    
    # Create horizontal bar chart
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    
    y_pos = np.arange(len(model_names))
    bars = ax.barh(y_pos, auroc_drops, alpha=0.85, color=colors, edgecolor='black', linewidth=1.2)
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels(model_names, fontsize=12)
    ax.set_xlabel('AUROC Performance Drop (%)', fontsize=14, fontweight='bold')
    ax.set_title('Cross-Species Performance Degradation (Human → Rat)', fontsize=16, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='x')
    ax.axvline(x=0, color='red', linestyle='--', linewidth=2, alpha=0.7)
    
    # Add value labels
    for i, (bar, val) in enumerate(zip(bars, auroc_drops)):
        ax.text(val + 1, bar.get_y() + bar.get_height()/2, 
                f'{val:.1f}%', 
                va='center', ha='left', fontsize=11, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(f"{output_path}/cross_species_performance_drop.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_path}/cross_species_performance_drop.pdf", bbox_inches='tight')
    plt.close()
    print(f"  ✓ Saved performance drop analysis")


def plot_coverage_comparison(human_results, rat_results, output_path):
    """Plot coverage comparison (models in task-specific then GFM order)"""
    # Get common models in MODEL_ORDER
    common_models = [m for m in MODEL_ORDER if m in human_results and m in rat_results]
    
    model_names = [MODEL_LABELS[m].replace('\n', ' ') for m in common_models]
    human_coverage = [human_results[m]['coverage'] for m in common_models]
    rat_coverage = [rat_results[m]['coverage'] for m in common_models]
    
    # Create grouped bar chart
    fig, ax = plt.subplots(1, 1, figsize=(14, 8))
    
    x = np.arange(len(model_names))
    width = 0.35
    
    bars1 = ax.bar(x - width/2, human_coverage, width, label='Human (SpliceVarDB)', 
                   alpha=0.85, color='#2E86AB', edgecolor='black', linewidth=1.2)
    bars2 = ax.bar(x + width/2, rat_coverage, width, label='Rat (RatGTEx)', 
                   alpha=0.85, color='#F77F00', edgecolor='black', linewidth=1.2)
    
    ax.set_xlabel('Models', fontsize=14, fontweight='bold')
    ax.set_ylabel('Coverage (%)', fontsize=14, fontweight='bold')
    ax.set_title('Cross-Species Coverage Comparison', fontsize=16, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(model_names, rotation=45, ha='right', fontsize=11)
    ax.legend(fontsize=13)
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_ylim(0, 105)
    
    # Add value labels
    for bar in bars1:
        height = bar.get_height()
        if height > 0:
            ax.annotate(f'{height:.1f}%',
                        xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 3),
                        textcoords="offset points",
                        ha='center', va='bottom', fontsize=9)
    
    for bar in bars2:
        height = bar.get_height()
        if height > 0:
            ax.annotate(f'{height:.1f}%',
                        xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 3),
                        textcoords="offset points",
                        ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    plt.savefig(f"{output_path}/cross_species_coverage_comparison.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_path}/cross_species_coverage_comparison.pdf", bbox_inches='tight')
    plt.close()
    print(f"  ✓ Saved coverage comparison")


def plot_three_way_comparison(human_results, rat_results, pig_results, chicken_results, output_path):
    """Plot four-way AUROC comparison across human, rat, pig, and chicken (models in task-specific then GFM order)"""
    # Get common models in MODEL_ORDER
    common_models = [m for m in MODEL_ORDER if m in human_results and m in rat_results and m in pig_results and m in chicken_results]
    
    model_names = [MODEL_LABELS[m].replace('\n', ' ') for m in common_models]
    human_auroc = [human_results[m]['auroc'] for m in common_models]
    rat_auroc = [rat_results[m]['auroc'] for m in common_models]
    pig_auroc = [pig_results[m]['auroc'] for m in common_models]
    chicken_auroc = [chicken_results[m]['auroc'] for m in common_models]
    
    # Create grouped bar chart
    fig, ax = plt.subplots(1, 1, figsize=(18, 8))
    
    x = np.arange(len(model_names))
    width = 0.2
    
    bars1 = ax.bar(x - 1.5*width, human_auroc, width, label='Human (SpliceVarDB)', 
                   alpha=0.85, color='#7B9FC2', edgecolor='black', linewidth=1.2)
    bars2 = ax.bar(x - 0.5*width, rat_auroc, width, label='Rat (RatGTEx)', 
                   alpha=0.85, color='#8B9DAF', edgecolor='black', linewidth=1.2)
    bars3 = ax.bar(x + 0.5*width, pig_auroc, width, label='Pig (PigGTEx)', 
                   alpha=0.85, color='#E8A5B8', edgecolor='black', linewidth=1.2)
    bars4 = ax.bar(x + 1.5*width, chicken_auroc, width, label='Chicken (ChickenGTEx)', 
                   alpha=0.85, color='#F4D35E', edgecolor='black', linewidth=1.2)
    
    ax.set_xlabel('Models', fontsize=14, fontweight='bold')
    ax.set_ylabel('AUROC', fontsize=14, fontweight='bold')
    ax.set_title('Cross-Species Performance Comparison: Human vs Rat vs Pig vs Chicken', fontsize=16, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(model_names, rotation=45, ha='right', fontsize=11)
    ax.legend(fontsize=13, loc='upper right')
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_ylim(0, 1.05)
    
    # Add value labels
    for bar in bars1:
        height = bar.get_height()
        ax.annotate(f'{height:.3f}',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=7, fontweight='bold')
    
    for bar in bars2:
        height = bar.get_height()
        ax.annotate(f'{height:.3f}',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=7, fontweight='bold')
    
    for bar in bars3:
        height = bar.get_height()
        ax.annotate(f'{height:.3f}',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=7, fontweight='bold')
    
    for bar in bars4:
        height = bar.get_height()
        ax.annotate(f'{height:.3f}',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=7, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(f"{output_path}/cross_species_auroc_comparison_four_way.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_path}/cross_species_auroc_comparison_four_way.pdf", bbox_inches='tight')
    plt.close()
    print(f"  ✓ Saved four-way AUROC comparison")


def plot_performance_drop_three_way(human_results, rat_results, pig_results, chicken_results, output_path):
    """Plot performance drop from human to rat, pig, and chicken (models in task-specific then GFM order)"""
    # Get common models in MODEL_ORDER
    common_models = [m for m in MODEL_ORDER if m in human_results and m in rat_results and m in pig_results and m in chicken_results]
    
    model_names = [MODEL_LABELS[m].replace('\n', ' ') for m in common_models]
    colors = [MODEL_COLORS[m] for m in common_models]
    
    # Calculate percentage drops
    auroc_drops_rat = [(human_results[m]['auroc'] - rat_results[m]['auroc']) / human_results[m]['auroc'] * 100 
                       for m in common_models]
    auroc_drops_pig = [(human_results[m]['auroc'] - pig_results[m]['auroc']) / human_results[m]['auroc'] * 100 
                       for m in common_models]
    auroc_drops_chicken = [(human_results[m]['auroc'] - chicken_results[m]['auroc']) / human_results[m]['auroc'] * 100 
                           for m in common_models]
    
    # Create grouped horizontal bar chart
    fig, ax = plt.subplots(1, 1, figsize=(14, 8))
    
    y_pos = np.arange(len(model_names))
    height = 0.25
    
    bars1 = ax.barh(y_pos - height, auroc_drops_rat, height, 
                    label='Human → Rat', alpha=0.85, color='#8B9DAF', edgecolor='black', linewidth=1.2)
    bars2 = ax.barh(y_pos, auroc_drops_pig, height,
                    label='Human → Pig', alpha=0.85, color='#E8A5B8', edgecolor='black', linewidth=1.2)
    bars3 = ax.barh(y_pos + height, auroc_drops_chicken, height,
                    label='Human → Chicken', alpha=0.85, color='#F4D35E', edgecolor='black', linewidth=1.2)
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels(model_names, fontsize=12)
    ax.set_xlabel('AUROC Performance Drop (%)', fontsize=14, fontweight='bold')
    ax.set_title('Cross-Species Performance Degradation', fontsize=16, fontweight='bold')
    ax.legend(fontsize=13, loc='lower right')
    ax.grid(True, alpha=0.3, axis='x')
    ax.axvline(x=0, color='red', linestyle='--', linewidth=2, alpha=0.7)
    
    # Add value labels
    for i, (bar, val) in enumerate(zip(bars1, auroc_drops_rat)):
        ax.text(val + 1, bar.get_y() + bar.get_height()/2, 
                f'{val:.1f}%', 
                va='center', ha='left', fontsize=9, fontweight='bold')
    
    for i, (bar, val) in enumerate(zip(bars2, auroc_drops_pig)):
        ax.text(val + 1, bar.get_y() + bar.get_height()/2, 
                f'{val:.1f}%', 
                va='center', ha='left', fontsize=9, fontweight='bold')
    
    for i, (bar, val) in enumerate(zip(bars3, auroc_drops_chicken)):
        ax.text(val + 1, bar.get_y() + bar.get_height()/2, 
                f'{val:.1f}%', 
                va='center', ha='left', fontsize=9, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(f"{output_path}/cross_species_performance_drop_four_way.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_path}/cross_species_performance_drop_four_way.pdf", bbox_inches='tight')
    plt.close()
    print(f"  ✓ Saved four-way performance drop analysis")


def plot_coverage_comparison_three_way(human_results, rat_results, pig_results, chicken_results, output_path):
    """Plot coverage comparison across all four species (models in task-specific then GFM order)"""
    # Get common models in MODEL_ORDER
    common_models = [m for m in MODEL_ORDER if m in human_results and m in rat_results and m in pig_results and m in chicken_results]
    
    model_names = [MODEL_LABELS[m].replace('\n', ' ') for m in common_models]
    human_coverage = [human_results[m]['coverage'] for m in common_models]
    rat_coverage = [rat_results[m]['coverage'] for m in common_models]
    pig_coverage = [pig_results[m]['coverage'] for m in common_models]
    chicken_coverage = [chicken_results[m]['coverage'] for m in common_models]
    
    # Create grouped bar chart
    fig, ax = plt.subplots(1, 1, figsize=(18, 8))
    
    x = np.arange(len(model_names))
    width = 0.2
    
    bars1 = ax.bar(x - 1.5*width, human_coverage, width, label='Human (SpliceVarDB)', 
                   alpha=0.85, color='#7B9FC2', edgecolor='black', linewidth=1.2)
    bars2 = ax.bar(x - 0.5*width, rat_coverage, width, label='Rat (RatGTEx)', 
                   alpha=0.85, color='#8B9DAF', edgecolor='black', linewidth=1.2)
    bars3 = ax.bar(x + 0.5*width, pig_coverage, width, label='Pig (PigGTEx)', 
                   alpha=0.85, color='#E8A5B8', edgecolor='black', linewidth=1.2)
    bars4 = ax.bar(x + 1.5*width, chicken_coverage, width, label='Chicken (ChickenGTEx)', 
                   alpha=0.85, color='#F4D35E', edgecolor='black', linewidth=1.2)
    
    ax.set_xlabel('Models', fontsize=14, fontweight='bold')
    ax.set_ylabel('Coverage (%)', fontsize=14, fontweight='bold')
    ax.set_title('Cross-Species Coverage Comparison', fontsize=16, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(model_names, rotation=45, ha='right', fontsize=11)
    ax.legend(fontsize=13)
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_ylim(0, 105)
    
    # Add value labels
    for bar in bars1:
        height = bar.get_height()
        if height > 0:
            ax.annotate(f'{height:.1f}%',
                        xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 3),
                        textcoords="offset points",
                        ha='center', va='bottom', fontsize=7)
    
    for bar in bars2:
        height = bar.get_height()
        if height > 0:
            ax.annotate(f'{height:.1f}%',
                        xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 3),
                        textcoords="offset points",
                        ha='center', va='bottom', fontsize=7)
    
    for bar in bars3:
        height = bar.get_height()
        if height > 0:
            ax.annotate(f'{height:.1f}%',
                        xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 3),
                        textcoords="offset points",
                        ha='center', va='bottom', fontsize=7)
    
    for bar in bars4:
        height = bar.get_height()
        if height > 0:
            ax.annotate(f'{height:.1f}%',
                        xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 3),
                        textcoords="offset points",
                        ha='center', va='bottom', fontsize=7)
    
    plt.tight_layout()
    plt.savefig(f"{output_path}/cross_species_coverage_comparison_four_way.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_path}/cross_species_coverage_comparison_four_way.pdf", bbox_inches='tight')
    plt.close()
    print(f"  ✓ Saved four-way coverage comparison")


def main():
    # Paths
    human_file = "/Users/byronsun/Desktop/AS_复现模型/BIB_review/code/splicevardb/baselines/splicevardb_all_model_predictions.csv"
    rat_file = "/Users/byronsun/Desktop/AS_复现模型/BIB_review/code/ratGTEx/baselines/ratgtex_all_model_predictions.csv"
    pig_file = "/Users/byronsun/Desktop/AS_复现模型/BIB_review/code/pigGTEx/baselines/piggtex_all_model_predictions.csv"
    chicken_file = "/Users/byronsun/Desktop/AS_复现模型/BIB_review/code/chickenGTEx/baselines/chickengtex_all_model_predictions.csv"
    output_path = "/Users/byronsun/Desktop/AS_复现模型/BIB_review/code/ratGTEx/baselines/visualisation"
    
    os.makedirs(output_path, exist_ok=True)
    
    print("📊 Loading data...")
    df_human = pd.read_csv(human_file)
    df_human['ground_truth_binary'] = df_human['ground_truth'].apply(lambda x: 1 if x == 'Splice-altering' else 0)
    df_rat = pd.read_csv(rat_file)
    df_pig = pd.read_csv(pig_file)
    df_chicken = pd.read_csv(chicken_file)
    
    print(f"Human: {len(df_human)} variants")
    print(f"Rat: {len(df_rat)} variants")
    print(f"Pig: {len(df_pig)} variants")
    print(f"Chicken: {len(df_chicken)} variants")
    
    print("\n📈 Calculating metrics...")
    human_results = calculate_metrics(df_human, is_human=True)
    rat_results = calculate_metrics(df_rat, is_human=False)
    pig_results = calculate_metrics(df_pig, is_human=False)
    chicken_results = calculate_metrics(df_chicken, is_human=False)
    
    print(f"Human models: {len(human_results)}")
    print(f"Rat models: {len(rat_results)}")
    print(f"Pig models: {len(pig_results)}")
    print(f"Chicken models: {len(chicken_results)}")
    
    print("\n🎨 Generating cross-species comparison plots...")
    
    print("1. Four-way species AUROC comparison...")
    plot_three_way_comparison(human_results, rat_results, pig_results, chicken_results, output_path)
    
    print("2. Performance drop analysis (Human → Rat, Pig, Chicken)...")
    plot_performance_drop_three_way(human_results, rat_results, pig_results, chicken_results, output_path)
    
    print("3. Coverage comparison (all four species)...")
    plot_coverage_comparison_three_way(human_results, rat_results, pig_results, chicken_results, output_path)
    
    print("\n✅ All cross-species comparison plots generated!")
    print(f"📁 Plots saved to: {output_path}")
    print("Files generated:")
    print("  - cross_species_auroc_comparison_four_way.png/pdf")
    print("  - cross_species_performance_drop_four_way.png/pdf")
    print("  - cross_species_coverage_comparison_four_way.png/pdf")


if __name__ == "__main__":
    main()
