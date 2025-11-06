#!/usr/bin/env python3
"""
Plot SpliceVarDB performance results for all models
Generate publication-ready figures for 4.2.1 Human Benchmark Results
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


def plot_roc_curves(df, output_path):
    """Plot ROC curves for all models"""
    fig, ax = plt.subplots(1, 1, figsize=(14, 11))
    
    models = ['AlphaGenome', 'Evo2_zeroshot', 'Nucleotide_Transformer', 'Pangolin', 
              'SpliceAI', 'SpliceBERT', 'SpliceTransformer',
              'DNABERT2_Logistic', 'MMSplice_pathogenicity']
    
    model_colors = {
        'AlphaGenome': '#1f77b4',
        'Evo2_zeroshot': '#ff7f0e', 
        'Nucleotide_Transformer': '#2ca02c',
        'Pangolin': '#d62728',
        'SpliceAI': '#9467bd',
        'SpliceBERT': '#8c564b',
        'SpliceTransformer': '#e377c2',
        'DNABERT2_Logistic': '#17becf',
        'MMSplice_pathogenicity': '#bcbd22'
    }
    
    model_labels = {
        'AlphaGenome': 'AlphaGenome',
        'Evo2_zeroshot': 'Evo2',
        'Nucleotide_Transformer': 'Nucleotide Transformer',
        'Pangolin': 'Pangolin',
        'SpliceAI': 'SpliceAI',
        'SpliceBERT': 'SpliceBERT',
        'SpliceTransformer': 'SpliceTransformer',
        'DNABERT2_Logistic': 'DNABERT-2 (Logistic)',
        'MMSplice_pathogenicity': 'MMSplice (pathogenicity)'
    }
    
    y_true = df['ground_truth_binary'].values
    
    for model in models:
        if model in df.columns:
            y_scores = df[model].values
            
            # Remove NaN values
            mask = ~np.isnan(y_scores)
            y_true_clean = y_true[mask]
            y_scores_clean = y_scores[mask]
            
            if len(y_true_clean) > 0:
                # For NT and SpliceBERT, we need to handle the scoring appropriately
                if model == 'Nucleotide_Transformer':
                    # NT uses cosine distance, higher distance = more different = more splice-altering
                    pass  # Keep as is since higher should mean more splice-altering
                elif model == 'SpliceBERT':
                    # SpliceBERT KL-context score is negative, invert it
                    y_scores_clean = -y_scores_clean
                
                # Calculate ROC curve
                fpr, tpr, _ = roc_curve(y_true_clean, y_scores_clean)
                auroc = roc_auc_score(y_true_clean, y_scores_clean)
                
                # Plot
                ax.plot(fpr, tpr, 
                       color=model_colors[model], 
                       linewidth=4.0,
                       label=f'{model_labels[model]} (AUROC = {auroc:.3f})')
    
    # Plot diagonal line
    ax.plot([0, 1], [0, 1], 'k--', alpha=0.5, linewidth=2.5)
    
    ax.set_xlabel('False Positive Rate', fontsize=22, fontweight='bold')
    ax.set_ylabel('True Positive Rate', fontsize=22, fontweight='bold')
    ax.set_title('ROC Curves - SpliceVarDB Human Benchmark', fontsize=24, fontweight='bold')
    ax.legend(loc='lower right', fontsize=20, frameon=True, fancybox=True, shadow=True, framealpha=0.95)
    ax.tick_params(axis='both', which='major', labelsize=18)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{output_path}/roc_curves.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_path}/roc_curves.pdf", bbox_inches='tight')
    plt.show()


def plot_precision_recall_curves(df, output_path):
    """Plot Precision-Recall curves for all models"""
    fig, ax = plt.subplots(1, 1, figsize=(14, 11))
    
    models = ['AlphaGenome', 'Evo2_zeroshot', 'Nucleotide_Transformer', 'Pangolin', 
              'SpliceAI', 'SpliceBERT', 'SpliceTransformer',
              'DNABERT2_Logistic', 'MMSplice_pathogenicity']
    
    model_colors = {
        'AlphaGenome': '#1f77b4',
        'Evo2_zeroshot': '#ff7f0e', 
        'Nucleotide_Transformer': '#2ca02c',
        'Pangolin': '#d62728',
        'SpliceAI': '#9467bd',
        'SpliceBERT': '#8c564b',
        'SpliceTransformer': '#e377c2',
        'DNABERT2_Logistic': '#17becf',
        'MMSplice_pathogenicity': '#bcbd22'
    }
    
    model_labels = {
        'AlphaGenome': 'AlphaGenome',
        'Evo2_zeroshot': 'Evo2',
        'Nucleotide_Transformer': 'Nucleotide Transformer',
        'Pangolin': 'Pangolin',
        'SpliceAI': 'SpliceAI',
        'SpliceBERT': 'SpliceBERT',
        'SpliceTransformer': 'SpliceTransformer',
        'DNABERT2_Logistic': 'DNABERT-2 (Logistic)',
        'MMSplice_pathogenicity': 'MMSplice (pathogenicity)'
    }
    
    y_true = df['ground_truth_binary'].values
    
    for model in models:
        if model in df.columns:
            y_scores = df[model].values
            
            # Remove NaN values
            mask = ~np.isnan(y_scores)
            y_true_clean = y_true[mask]
            y_scores_clean = y_scores[mask]
            
            if len(y_true_clean) > 0:
                # Handle scoring for different models
                if model == 'Nucleotide_Transformer':
                    # NT uses cosine distance, higher distance = more different = more splice-altering
                    pass  # Keep as is
                elif model == 'SpliceBERT':
                    # SpliceBERT KL-context score is negative, invert it
                    y_scores_clean = -y_scores_clean
                
                # Calculate Precision-Recall curve
                precision, recall, _ = precision_recall_curve(y_true_clean, y_scores_clean)
                auprc = average_precision_score(y_true_clean, y_scores_clean)
                
                # Plot
                ax.plot(recall, precision, 
                       color=model_colors[model], 
                       linewidth=4.0,
                       label=f'{model_labels[model]} (AUPRC = {auprc:.3f})')
    
    # Plot baseline (random classifier)
    baseline = sum(y_true) / len(y_true)
    ax.axhline(y=baseline, color='k', linestyle='--', alpha=0.5, linewidth=2.5, 
               label=f'Random (AUPRC = {baseline:.3f})')
    
    ax.set_xlabel('Recall', fontsize=22, fontweight='bold')
    ax.set_ylabel('Precision', fontsize=22, fontweight='bold')
    ax.set_title('Precision-Recall Curves - SpliceVarDB Human Benchmark', fontsize=24, fontweight='bold')
    ax.legend(loc='lower right', fontsize=20, frameon=True, fancybox=True, shadow=True, framealpha=0.95)
    ax.tick_params(axis='both', which='major', labelsize=18)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{output_path}/precision_recall_curves.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_path}/precision_recall_curves.pdf", bbox_inches='tight')
    plt.show()


def plot_performance_bar_chart(df, output_path):
    """Plot bar chart comparing AUROC and AUPRC"""
    models = ['AlphaGenome', 'Evo2_zeroshot', 'Nucleotide_Transformer', 'Pangolin', 
              'SpliceAI', 'SpliceBERT', 'SpliceTransformer',
              'DNABERT2_Logistic', 'MMSplice_pathogenicity']
    
    model_labels = {
        'AlphaGenome': 'AlphaGenome',
        'Evo2_zeroshot': 'Evo2',
        'Nucleotide_Transformer': 'Nucleotide\nTransformer',
        'Pangolin': 'Pangolin',
        'SpliceAI': 'SpliceAI',
        'SpliceBERT': 'SpliceBERT',
        'SpliceTransformer': 'SpliceTransformer',
        'DNABERT2_Logistic': 'DNABERT-2\n(Logistic)',
        'MMSplice_pathogenicity': 'MMSplice\n(pathogenicity)'
    }
    
    y_true = df['ground_truth_binary'].values
    
    auroc_scores = []
    auprc_scores = []
    model_names = []
    
    for model in models:
        if model in df.columns:
            y_scores = df[model].values
            
            # Remove NaN values
            mask = ~np.isnan(y_scores)
            y_true_clean = y_true[mask]
            y_scores_clean = y_scores[mask]
            
            if len(y_true_clean) > 0:
                # Handle scoring for different models
                if model == 'Nucleotide_Transformer':
                    pass  # Keep as is
                elif model == 'SpliceBERT':
                    y_scores_clean = -y_scores_clean
                
                auroc = roc_auc_score(y_true_clean, y_scores_clean)
                auprc = average_precision_score(y_true_clean, y_scores_clean)
                
                auroc_scores.append(auroc)
                auprc_scores.append(auprc)
                model_names.append(model_labels[model])
    
    # Create bar chart
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    
    x = np.arange(len(model_names))
    width = 0.35
    
    bars1 = ax.bar(x - width/2, auroc_scores, width, label='AUROC', alpha=0.8, color='#2E86AB')
    bars2 = ax.bar(x + width/2, auprc_scores, width, label='AUPRC', alpha=0.8, color='#A23B72')
    
    ax.set_xlabel('Models', fontsize=14)
    ax.set_ylabel('Performance Score', fontsize=14)
    ax.set_title('Model Performance Comparison - SpliceVarDB Human Benchmark', fontsize=16, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(model_names, rotation=45, ha='right')
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
    plt.savefig(f"{output_path}/performance_comparison.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_path}/performance_comparison.pdf", bbox_inches='tight')
    plt.show()


def plot_coverage_analysis(df, output_path):
    """Plot coverage analysis for all models"""
    models = ['AlphaGenome', 'Evo2_zeroshot', 'Nucleotide_Transformer', 'Pangolin', 
              'SpliceAI', 'SpliceBERT', 'SpliceTransformer',
              'DNABERT2_Logistic', 'MMSplice_pathogenicity']
    
    model_labels = {
        'AlphaGenome': 'AlphaGenome',
        'Evo2_zeroshot': 'Evo2',
        'Nucleotide_Transformer': 'Nucleotide\nTransformer',
        'Pangolin': 'Pangolin',
        'SpliceAI': 'SpliceAI',
        'SpliceBERT': 'SpliceBERT',
        'SpliceTransformer': 'SpliceTransformer',
        'DNABERT2_Logistic': 'DNABERT-2\n(Logistic)',
        'MMSplice_pathogenicity': 'MMSplice\n(pathogenicity)'
    }
    
    coverage_data = []
    model_names = []
    
    total_variants = len(df)
    
    for model in models:
        if model in df.columns:
            predictions = df[model].notna().sum()
            coverage_pct = (predictions / total_variants) * 100
            
            coverage_data.append(coverage_pct)
            model_names.append(model_labels[model])
    
    # Create coverage plot
    fig, ax = plt.subplots(1, 1, figsize=(12, 6))
    
    bars = ax.bar(model_names, coverage_data, alpha=0.8, color='#F18F01')
    
    ax.set_xlabel('Models', fontsize=14)
    ax.set_ylabel('Coverage (%)', fontsize=14)
    ax.set_title('Model Coverage - SpliceVarDB Human Benchmark', fontsize=16, fontweight='bold')
    ax.set_ylim(95, 101)
    ax.grid(True, alpha=0.3, axis='y')
    
    # Add value labels on bars
    for bar in bars:
        height = bar.get_height()
        ax.annotate(f'{height:.2f}%',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=11)
    
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(f"{output_path}/coverage_analysis.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_path}/coverage_analysis.pdf", bbox_inches='tight')
    plt.show()

def plot_threshold_metrics(df, output_path):
    """Plot Accuracy, Precision, Recall, F1 for all models using optimal thresholds"""
    models = ['AlphaGenome', 'Evo2_zeroshot', 'Nucleotide_Transformer', 'Pangolin', 
              'SpliceAI', 'SpliceBERT', 'SpliceTransformer']
    
    model_labels = {
        'AlphaGenome': 'AlphaGenome',
        'Evo2_zeroshot': 'Evo2',
        'Nucleotide_Transformer': 'Nucleotide\nTransformer',
        'Pangolin': 'Pangolin',
        'SpliceAI': 'SpliceAI',
        'SpliceBERT': 'SpliceBERT',
        'SpliceTransformer': 'SpliceTransformer',
        'DNABERT2_Logistic': 'DNABERT-2\n(Logistic)',
        'MMSplice_pathogenicity': 'MMSplice\n(pathogenicity)'
    }
    
    y_true = df['ground_truth_binary'].values
    
    # Calculate metrics for each model
    metrics_data = {
        'Accuracy': [],
        'Precision': [],
        'Recall': [],
        'F1': []
    }
    model_names = []
    
    for model in models:
        if model in df.columns:
            y_scores = df[model].values
            
            # Remove NaN values
            mask = ~np.isnan(y_scores)
            y_true_clean = y_true[mask]
            y_scores_clean = y_scores[mask]
            
            if len(y_true_clean) > 0:
                # Handle scoring for different models
                if model == 'Nucleotide_Transformer':
                    pass  # Keep as is
                elif model == 'SpliceBERT':
                    y_scores_clean = -y_scores_clean
                
                # Find optimal threshold using F1-score, with special handling for poor-performing models
                precision_curve, recall_curve, thresholds = precision_recall_curve(y_true_clean, y_scores_clean)
                
                if len(thresholds) > 0:
                    f1_scores = 2 * precision_curve[:-1] * recall_curve[:-1] / (precision_curve[:-1] + recall_curve[:-1] + 1e-8)
                    optimal_idx = np.argmax(f1_scores)
                    optimal_threshold = thresholds[optimal_idx]
                    
                    # Check for problematic models and adjust threshold
                    if model in ['Evo2_zeroshot', 'Nucleotide_Transformer', 'SpliceBERT']:
                        # For models with poor discrimination, use median or percentile-based threshold
                        if model == 'Evo2_zeroshot':
                            # Use 75th percentile for Evo2 to avoid too many zero predictions
                            optimal_threshold = np.percentile(y_scores_clean[y_scores_clean > 0], 75)
                        elif model == 'Nucleotide_Transformer':
                            # Use median for NT since scores are all near zero
                            optimal_threshold = np.median(y_scores_clean)
                        elif model == 'SpliceBERT':
                            # Use 25th percentile for SpliceBERT (inverted scores)
                            optimal_threshold = np.percentile(y_scores_clean, 25)
                else:
                    optimal_threshold = 0.5
                
                # Make binary predictions using optimal threshold
                y_pred = (y_scores_clean >= optimal_threshold).astype(int)
                
                # Calculate metrics
                accuracy = accuracy_score(y_true_clean, y_pred)
                precision = precision_score(y_true_clean, y_pred, zero_division=0)
                recall = recall_score(y_true_clean, y_pred, zero_division=0)
                f1 = f1_score(y_true_clean, y_pred, zero_division=0)
                
                metrics_data['Accuracy'].append(accuracy)
                metrics_data['Precision'].append(precision)
                metrics_data['Recall'].append(recall)
                metrics_data['F1'].append(f1)
                model_names.append(model_labels[model])
    
    # Create subplots for each metric
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    axes = axes.flatten()
    
    metrics = ['Accuracy', 'Precision', 'Recall', 'F1']
    colors = ['#2E86AB', '#A23B72', '#F18F01', '#54A24B']
    
    for i, metric in enumerate(metrics):
        ax = axes[i]
        
        bars = ax.bar(model_names, metrics_data[metric], alpha=0.8, color=colors[i])
        
        ax.set_xlabel('Models', fontsize=12)
        ax.set_ylabel(metric, fontsize=12)
        ax.set_title(f'{metric} - SpliceVarDB Human Benchmark', fontsize=14, fontweight='bold')
        ax.set_ylim(0, 1.05)
        ax.grid(True, alpha=0.3, axis='y')
        
        # Add value labels on bars
        for bar in bars:
            height = bar.get_height()
            ax.annotate(f'{height:.3f}',
                        xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 3),
                        textcoords="offset points",
                        ha='center', va='bottom', fontsize=10)
        
        ax.tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig(f"{output_path}/threshold_metrics.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_path}/threshold_metrics.pdf", bbox_inches='tight')
    plt.show()


def plot_grouped_threshold_metrics(df, output_path):
    """Plot grouped bar chart with metrics as groups and models as bars within each group"""
    models = ['AlphaGenome', 'Evo2_zeroshot', 'Nucleotide_Transformer', 'Pangolin', 
              'SpliceAI', 'SpliceBERT', 'SpliceTransformer',
              'DNABERT2_Logistic', 'MMSplice_pathogenicity']
    
    model_labels = ['AlphaGenome', 'Evo2', 'Nucleotide\nTransformer', 'Pangolin', 
                   'SpliceAI', 'SpliceBERT', 'SpliceTransformer',
                   'DNABERT-2\n(Logistic)', 'MMSplice\n(pathogenicity)']
    
    y_true = df['ground_truth_binary'].values
    
    # Calculate metrics for each model
    metrics_data = {
        'Accuracy': [],
        'Precision': [],
        'Recall': [],
        'F1': []
    }
    
    for model in models:
        if model in df.columns:
            y_scores = df[model].values
            
            # Remove NaN values
            mask = ~np.isnan(y_scores)
            y_true_clean = y_true[mask]
            y_scores_clean = y_scores[mask]
            
            if len(y_true_clean) > 0:
                # Handle scoring for different models
                if model == 'Nucleotide_Transformer':
                    pass  # Keep as is
                elif model == 'SpliceBERT':
                    y_scores_clean = -y_scores_clean
                
                # Find optimal threshold using F1-score, with special handling for poor-performing models
                precision_curve, recall_curve, thresholds = precision_recall_curve(y_true_clean, y_scores_clean)
                
                if len(thresholds) > 0:
                    f1_curve = 2 * precision_curve[:-1] * recall_curve[:-1] / (precision_curve[:-1] + recall_curve[:-1] + 1e-8)
                    optimal_idx = np.argmax(f1_curve)
                    optimal_threshold = thresholds[optimal_idx]
                    
                    # Check for problematic models and adjust threshold
                    if model in ['Evo2_zeroshot', 'Nucleotide_Transformer', 'SpliceBERT']:
                        # For models with poor discrimination, use median or percentile-based threshold
                        if model == 'Evo2_zeroshot':
                            # Use 75th percentile for Evo2 to avoid too many zero predictions
                            optimal_threshold = np.percentile(y_scores_clean[y_scores_clean > 0], 75)
                        elif model == 'Nucleotide_Transformer':
                            # Use median for NT since scores are all near zero
                            optimal_threshold = np.median(y_scores_clean)
                        elif model == 'SpliceBERT':
                            # Use 25th percentile for SpliceBERT (inverted scores)
                            optimal_threshold = np.percentile(y_scores_clean, 25)
                else:
                    optimal_threshold = 0.5
                
                # Make binary predictions using optimal threshold
                y_pred = (y_scores_clean >= optimal_threshold).astype(int)
                
                # Calculate metrics
                accuracy = accuracy_score(y_true_clean, y_pred)
                precision = precision_score(y_true_clean, y_pred, zero_division=0)
                recall = recall_score(y_true_clean, y_pred, zero_division=0)
                f1 = f1_score(y_true_clean, y_pred, zero_division=0)
                
                metrics_data['Accuracy'].append(accuracy)
                metrics_data['Precision'].append(precision)
                metrics_data['Recall'].append(recall)
                metrics_data['F1'].append(f1)
    
    # Create grouped bar chart with metrics as groups
    fig, ax = plt.subplots(1, 1, figsize=(16, 8))
    
    metrics = ['Accuracy', 'Precision', 'Recall', 'F1']
    n_models = len(model_labels)
    n_metrics = len(metrics)
    
    # Set the width of bars and positions
    bar_width = 0.1
    spacing = 0.02  # Space between bars within a group
    group_width = n_models * bar_width + (n_models - 1) * spacing
    group_spacing = 0.3  # Space between groups
    
    # Calculate x positions for each metric group
    metric_positions = []
    for i in range(n_metrics):
        metric_positions.append(i * (group_width + group_spacing))
    
    # Colors for each model
    model_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#17becf', '#bcbd22']
    
    # Plot bars for each model in each metric group
    for model_idx, model_label in enumerate(model_labels):
        model_x_positions = []
        model_values = []
        
        for metric_idx, metric in enumerate(metrics):
            x_pos = metric_positions[metric_idx] + model_idx * (bar_width + spacing)
            model_x_positions.append(x_pos)
            model_values.append(metrics_data[metric][model_idx])
        
        bars = ax.bar(model_x_positions, model_values, bar_width, 
                     label=model_label, alpha=0.8, color=model_colors[model_idx])
        
        # Add value labels on bars
        for bar, value in zip(bars, model_values):
            height = bar.get_height()
            ax.annotate(f'{value:.3f}',
                       xy=(bar.get_x() + bar.get_width() / 2, height),
                       xytext=(0, 3),
                       textcoords="offset points",
                       ha='center', va='bottom', fontsize=8, rotation=90)
    
    # Set x-axis labels and ticks
    group_centers = [pos + group_width/2 - (n_models-1)*(bar_width+spacing)/2 for pos in metric_positions]
    ax.set_xticks(group_centers)
    ax.set_xticklabels(metrics, fontsize=14)
    
    ax.set_xlabel('Metrics', fontsize=20, fontweight='bold')
    ax.set_ylabel('Score', fontsize=20, fontweight='bold')
    ax.set_title('Threshold-Based Metrics Comparison - SpliceVarDB Human Benchmark', fontsize=22, fontweight='bold')
    ax.legend(bbox_to_anchor=(1.02, 0.5), loc='center left', fontsize=18, frameon=True, fancybox=True, shadow=True, framealpha=0.95)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_ylim(0, 1.05)
    
    plt.tight_layout()
    plt.savefig(f"{output_path}/grouped_threshold_metrics.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_path}/grouped_threshold_metrics.pdf", bbox_inches='tight')
    plt.show()


def main():
    # Load data
    data_file = "/Users/byronsun/Desktop/AS_复现模型/BIB_review/code/splicevardb/baselines/splicevardb_all_model_predictions.csv"
    output_path = "/Users/byronsun/Desktop/AS_复现模型/BIB_review/code/splicevardb/baselines/visualisation"
    
    # Create output directory if it doesn't exist
    os.makedirs(output_path, exist_ok=True)
    
    print("📊 Loading SpliceVarDB results...")
    df = pd.read_csv(data_file)
    
    # Convert ground_truth to binary
    df['ground_truth_binary'] = df['ground_truth'].apply(lambda x: 1 if x == 'Splice-altering' else 0)
    
    print(f"Loaded {len(df)} variants")
    print(f"Splice-altering: {df['ground_truth_binary'].sum()} ({df['ground_truth_binary'].mean()*100:.1f}%)")
    print(f"Normal: {(1 - df['ground_truth_binary']).sum()} ({(1 - df['ground_truth_binary'].mean())*100:.1f}%)")
    
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
    
    print("5. Threshold-Based Metrics (Individual)...")
    plot_threshold_metrics(df, output_path)
    
    print("6. Threshold-Based Metrics (Grouped)...")
    plot_grouped_threshold_metrics(df, output_path)
    
    print("\n✅ All plots generated successfully!")
    print(f"📁 Plots saved to: {output_path}")
    print("Files generated:")
    print("  - roc_curves.png/pdf")
    print("  - precision_recall_curves.png/pdf") 
    print("  - performance_comparison.png/pdf")
    print("  - coverage_analysis.png/pdf")
    print("  - threshold_metrics.png/pdf")
    print("  - grouped_threshold_metrics.png/pdf")


if __name__ == "__main__":
    main()


