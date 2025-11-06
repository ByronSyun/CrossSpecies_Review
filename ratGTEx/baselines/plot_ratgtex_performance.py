#!/usr/bin/env python3
"""
Plot RatGTEx performance results for all models
Generate publication-ready figures for 4.2.2 Non-human Benchmark Results
Colors are matched with SpliceVarDB plots for consistency
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_curve, precision_recall_curve, auc, roc_auc_score, average_precision_score, accuracy_score, precision_score, recall_score, f1_score
import warnings
import os
warnings.filterwarnings('ignore')

# Set style for publication
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

# MODEL COLORS - MUST MATCH SPLICEVARDB for consistency across datasets
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


def plot_roc_curves(df, output_path):
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
            
            # Remove NaN values
            mask = ~np.isnan(y_scores)
            y_true_clean = y_true[mask]
            y_scores_clean = y_scores[mask]
            
            if len(y_true_clean) > 0 and len(np.unique(y_true_clean)) > 1:
                try:
                    # Calculate ROC curve
                    fpr, tpr, _ = roc_curve(y_true_clean, y_scores_clean)
                    auroc = roc_auc_score(y_true_clean, y_scores_clean)
                    
                    # Plot
                    ax.plot(fpr, tpr, 
                           color=MODEL_COLORS[model], 
                           linewidth=4.0,
                           label=f'{MODEL_LABELS[model]} (AUROC = {auroc:.3f})')
                except Exception as e:
                    print(f"Warning: Could not plot {model}: {e}")
    
    # Plot diagonal line
    ax.plot([0, 1], [0, 1], 'k--', alpha=0.5, linewidth=2.5)
    
    ax.set_xlabel('False Positive Rate', fontsize=22, fontweight='bold')
    ax.set_ylabel('True Positive Rate', fontsize=22, fontweight='bold')
    ax.set_title('ROC Curves - RatGTEx Cross-Species Benchmark', fontsize=24, fontweight='bold')
    ax.legend(loc='lower right', fontsize=20, frameon=True, fancybox=True, shadow=True, framealpha=0.95)
    ax.tick_params(axis='both', which='major', labelsize=18)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{output_path}/rat_roc_curves.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_path}/rat_roc_curves.pdf", bbox_inches='tight')
    plt.close()
    print(f"  ✓ Saved ROC curves")


def plot_precision_recall_curves(df, output_path):
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
            
            # Remove NaN values
            mask = ~np.isnan(y_scores)
            y_true_clean = y_true[mask]
            y_scores_clean = y_scores[mask]
            
            if len(y_true_clean) > 0 and len(np.unique(y_true_clean)) > 1:
                try:
                    # Calculate Precision-Recall curve
                    precision, recall, _ = precision_recall_curve(y_true_clean, y_scores_clean)
                    auprc = average_precision_score(y_true_clean, y_scores_clean)
                    
                    # Plot
                    ax.plot(recall, precision, 
                           color=MODEL_COLORS[model], 
                           linewidth=4.0,
                           label=f'{MODEL_LABELS[model]} (AUPRC = {auprc:.3f})')
                except Exception as e:
                    print(f"Warning: Could not plot {model}: {e}")
    
    # Plot baseline (random classifier) - no label in legend
    ax.axhline(y=baseline_precision, color='k', linestyle='--', alpha=0.5, linewidth=2.5)
    
    ax.set_xlabel('Recall', fontsize=26, fontweight='bold')
    ax.set_ylabel('Precision', fontsize=26, fontweight='bold')
    ax.set_title('Precision-Recall Curves - RatGTEx Cross-Species Benchmark', fontsize=28, fontweight='bold')
    ax.legend(loc='lower right', fontsize=28, frameon=True, fancybox=True, shadow=True, framealpha=0.95)
    ax.tick_params(axis='both', which='major', labelsize=22, width=2)
    ax.grid(True, alpha=0.3)
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1.05])
    
    plt.tight_layout()
    plt.savefig(f"{output_path}/rat_precision_recall_curves.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_path}/rat_precision_recall_curves.pdf", bbox_inches='tight')
    plt.close()
    print(f"  ✓ Saved Precision-Recall curves")


def plot_performance_bar_chart(df, output_path):
    """Plot AUROC and AUPRC comparison as grouped bar chart"""
    models = ['AlphaGenome', 'Pangolin', 'SpliceAI', 'Nucleotide_Transformer',
              'SpliceBERT', 'SpliceTransformer', 'DNABERT2_Logistic', 
              'Evo2', 'MMSplice_pathogenicity']
    
    y_true = df['ground_truth'].values
    
    auroc_scores = []
    auprc_scores = []
    model_names = []
    
    for model in models:
        if model in df.columns:
            y_scores = df[model].values
            
            mask = ~np.isnan(y_scores)
            y_true_clean = y_true[mask]
            y_scores_clean = y_scores[mask]
            
            if len(y_true_clean) > 0 and len(np.unique(y_true_clean)) > 1:
                if model == 'SpliceBERT':
                    y_scores_clean = -y_scores_clean
                
                try:
                    auroc = roc_auc_score(y_true_clean, y_scores_clean)
                    auprc = average_precision_score(y_true_clean, y_scores_clean)
                    
                    auroc_scores.append(auroc)
                    auprc_scores.append(auprc)
                    model_names.append(MODEL_LABELS[model])
                except Exception as e:
                    print(f"Warning: Could not calculate metrics for {model}: {e}")
    
    # Create bar chart
    fig, ax = plt.subplots(1, 1, figsize=(14, 8))
    
    x = np.arange(len(model_names))
    width = 0.35
    
    bars1 = ax.bar(x - width/2, auroc_scores, width, label='AUROC', alpha=0.8, color='#2E86AB')
    bars2 = ax.bar(x + width/2, auprc_scores, width, label='AUPRC', alpha=0.8, color='#A23B72')
    
    ax.set_xlabel('Models', fontsize=14, fontweight='bold')
    ax.set_ylabel('Performance Score', fontsize=14, fontweight='bold')
    ax.set_title('Model Performance Comparison - RatGTEx Cross-Species Benchmark', fontsize=16, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(model_names, rotation=45, ha='right', fontsize=12)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_ylim(0, 1.05)
    
    # Add value labels on bars
    for bar in bars1:
        height = bar.get_height()
        ax.annotate(f'{height:.3f}',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=10)
    
    for bar in bars2:
        height = bar.get_height()
        ax.annotate(f'{height:.3f}',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=10)
    
    plt.tight_layout()
    plt.savefig(f"{output_path}/rat_performance_comparison.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_path}/rat_performance_comparison.pdf", bbox_inches='tight')
    plt.close()
    print(f"  ✓ Saved performance comparison bar chart")


def plot_coverage_analysis(df, output_path):
    """Plot coverage analysis for all models"""
    models = ['AlphaGenome', 'Pangolin', 'SpliceAI', 'SpliceBERT', 'Evo2',
              'Nucleotide_Transformer', 'MMSplice_pathogenicity', 
              'DNABERT2_Logistic', 'SpliceTransformer']
    
    coverage_data = []
    model_names = []
    colors = []
    
    total_variants = len(df)
    
    for model in models:
        if model in df.columns:
            predictions = df[model].notna().sum()
            coverage_pct = (predictions / total_variants) * 100
            
            coverage_data.append(coverage_pct)
            model_names.append(MODEL_LABELS[model].replace('\n', ' '))
            colors.append(MODEL_COLORS[model])
    
    # Create coverage plot
    fig, ax = plt.subplots(1, 1, figsize=(14, 6))
    
    bars = ax.barh(model_names, coverage_data, alpha=0.8, color=colors, edgecolor='black', linewidth=1.2)
    
    ax.set_xlabel('Coverage (%)', fontsize=14, fontweight='bold')
    ax.set_title('Model Coverage - RatGTEx Cross-Species Benchmark', fontsize=16, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='x')
    ax.set_xlim(0, 105)
    
    # Add value labels
    for i, (bar, val) in enumerate(zip(bars, coverage_data)):
        ax.text(val + 1, bar.get_y() + bar.get_height()/2, 
                f'{val:.1f}%', 
                va='center', ha='left', fontsize=11, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(f"{output_path}/rat_coverage_analysis.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_path}/rat_coverage_analysis.pdf", bbox_inches='tight')
    plt.close()
    print(f"  ✓ Saved coverage analysis")


def main():
    # Paths
    data_file = "/Users/byronsun/Desktop/AS_复现模型/BIB_review/code/ratGTEx/baselines/ratgtex_all_model_predictions.csv"
    output_path = "/Users/byronsun/Desktop/AS_复现模型/BIB_review/code/ratGTEx/baselines/visualisation"
    
    # Create output directory
    os.makedirs(output_path, exist_ok=True)
    
    print("📊 Loading RatGTEx results...")
    df = pd.read_csv(data_file)
    
    print(f"Loaded {len(df)} variants")
    print(f"Splice-altering: {df['ground_truth'].sum()} ({df['ground_truth'].mean()*100:.1f}%)")
    print(f"Normal: {(1 - df['ground_truth']).sum()} ({(1 - df['ground_truth'].mean())*100:.1f}%)")
    
    print("\n🎨 Generating performance plots...")
    
    # Generate all plots
    print("1. ROC Curves...")
    plot_roc_curves(df, output_path)
    
    print("2. Precision-Recall Curves...")
    plot_precision_recall_curves(df, output_path)
    
    print("3. Performance Comparison Bar Chart (AUROC & AUPRC)...")
    plot_performance_bar_chart(df, output_path)
    
    print("4. Coverage Analysis...")
    plot_coverage_analysis(df, output_path)
    
    print("\n✅ All plots generated successfully!")
    print(f"📁 Plots saved to: {output_path}")
    print("Files generated:")
    print("  - rat_roc_curves.png/pdf")
    print("  - rat_precision_recall_curves.png/pdf")
    print("  - rat_performance_comparison.png/pdf")
    print("  - rat_coverage_analysis.png/pdf")


if __name__ == "__main__":
    main()

