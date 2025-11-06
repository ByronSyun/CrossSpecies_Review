#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot AlphaGenome cross-species performance comparison
Shows AUROC on human (SpliceVarDB) vs rat (ratGTEx) datasets
"""

import os
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve, confusion_matrix, precision_recall_curve

# Configure plotting style
plt.style.use('default')
sns.set_palette("husl")

def load_splicevardb_results(json_path):
    """Load SpliceVarDB AlphaGenome results."""
    with open(json_path, 'r') as f:
        data = json.load(f)
    
    # Extract successful predictions
    successful = [r for r in data['results'] if r['status'] == 'success']
    scores = np.array([r['analysis']['differential']['max_abs_change'] for r in successful])
    labels = np.array([1 if r['classification'].lower() == 'splice-altering' else 0 for r in successful])
    
    return scores, labels

def load_ratgtex_results(tsv_path):
    """Load ratGTEx AlphaGenome results from the final TSV output."""
    df = pd.read_csv(tsv_path, sep='\t')
    
    # Filter out any rows with processing errors (score == -1.0)
    valid_df = df[df['alphagenome_score'] != -1.0].copy()
    
    scores = valid_df['alphagenome_score'].values
    labels = valid_df['label'].values
    
    return scores, labels

def main():
    # File paths
    splicevardb_path = "/Users/byronsun/Desktop/AS_复现模型/BIB_review/code/splicevardb/baselines/alphagenome/alphagenome_COMPLETE_results_20250904_094904.json"
    ratgtex_path = "/Users/byronsun/Desktop/AS_复现模型/BIB_review/code/ratGTEx/baselines/alphagenome/alphagenome_scores_20250915_125255.tsv"
    output_dir = "/Users/byronsun/Desktop/AS_复现模型/BIB_review/code/ratGTEx/baselines/alphagenome/"
    
    # Load results
    print("Loading SpliceVarDB (human) results...")
    human_scores, human_labels = load_splicevardb_results(splicevardb_path)
    
    print("Loading ratGTEx (rat) results...")
    rat_scores, rat_labels = load_ratgtex_results(ratgtex_path)
    
    # Calculate metrics
    human_auroc = roc_auc_score(human_labels, human_scores)
    human_auprc = average_precision_score(human_labels, human_scores)
    
    rat_auroc = roc_auc_score(rat_labels, rat_scores)
    rat_auprc = average_precision_score(rat_labels, rat_scores)
    
    print(f"\n=== Performance Summary ===")
    print(f"Human (SpliceVarDB): AUROC={human_auroc:.4f}, AUPRC={human_auprc:.4f}")
    print(f"Rat (ratGTEx):       AUROC={rat_auroc:.4f}, AUPRC={rat_auprc:.4f}")
    print(f"Performance gap:     ΔAUROC={human_auroc - rat_auroc:.4f}")
    
    # Create comparison plots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
    
    # 1. AUROC Bar Chart
    species = ['Human\n(SpliceVarDB)', 'Rat\n(ratGTEx)']
    aurocs = [human_auroc, rat_auroc]
    colors = ['#4C78A8', '#E45756']
    
    bars = ax1.bar(species, aurocs, color=colors, alpha=0.8)
    ax1.set_ylabel('AUROC', fontsize=12)
    ax1.set_title('AlphaGenome Cross-Species Performance', fontsize=14, fontweight='bold')
    ax1.set_ylim(0, 1.0)
    ax1.grid(axis='y', alpha=0.3)
    
    # Add value labels on bars
    for bar, val in zip(bars, aurocs):
        ax1.text(bar.get_x() + bar.get_width()/2, val + 0.02, 
                f'{val:.3f}', ha='center', va='bottom', fontweight='bold')
    
    # 2. Rat Confusion Matrix
    # Find F1-optimal threshold for rat data
    precision, recall, thresholds = precision_recall_curve(rat_labels, rat_scores)
    f1_scores = 2 * precision[:-1] * recall[:-1] / (precision[:-1] + recall[:-1] + 1e-12)
    optimal_idx = np.argmax(f1_scores)
    optimal_threshold = thresholds[optimal_idx]
    
    rat_pred = (rat_scores >= optimal_threshold).astype(int)
    cm = confusion_matrix(rat_labels, rat_pred)
    
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=ax2,
                xticklabels=['Negative', 'Positive'],
                yticklabels=['Negative', 'Positive'],
                annot_kws={"size": 12})
    ax2.set_title(f'Rat (ratGTEx) Confusion Matrix\nThreshold: {optimal_threshold:.6f}', fontsize=14, fontweight='bold')
    ax2.set_xlabel('Predicted Label', fontsize=12)
    ax2.set_ylabel('True Label', fontsize=12)
    
    # 3. ROC Curves
    human_fpr, human_tpr, _ = roc_curve(human_labels, human_scores)
    rat_fpr, rat_tpr, _ = roc_curve(rat_labels, rat_scores)
    
    ax3.plot(human_fpr, human_tpr, color=colors[0], linewidth=2, 
             label=f'Human (AUROC={human_auroc:.3f})')
    ax3.plot(rat_fpr, rat_tpr, color=colors[1], linewidth=2, 
             label=f'Rat (AUROC={rat_auroc:.3f})')
    ax3.plot([0, 1], [0, 1], 'k--', alpha=0.5, linewidth=1)
    ax3.set_xlabel('False Positive Rate', fontsize=12)
    ax3.set_ylabel('True Positive Rate', fontsize=12)
    ax3.set_title('ROC Curves Comparison', fontsize=14, fontweight='bold')
    ax3.legend(fontsize=10)
    ax3.grid(alpha=0.3)
    
    # 4. Score Distributions
    bins = np.linspace(0, max(human_scores.max(), rat_scores.max()), 50)
    
    ax4.hist(human_scores, bins=bins, alpha=0.6, color=colors[0], 
             label=f'Human (n={len(human_scores):,})', density=True)
    ax4.hist(rat_scores, bins=bins, alpha=0.6, color=colors[1], 
             label=f'Rat (n={len(rat_scores):,})', density=True)
    ax4.set_xlabel('AlphaGenome Score', fontsize=12)
    ax4.set_ylabel('Density', fontsize=12)
    ax4.set_title('Score Distribution Comparison', fontsize=14, fontweight='bold')
    ax4.legend(fontsize=10)
    ax4.grid(alpha=0.3)
    
    plt.tight_layout()
    
    # Save plots
    output_path = os.path.join(output_dir, 'alphagenome_cross_species_performance.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.replace('.png', '.pdf'), bbox_inches='tight')
    
    print(f"\nPlots saved to: {output_path}")
    
    # Create summary table
    tn, fp, fn, tp = cm.ravel()
    accuracy = (tp + tn) / (tp + tn + fp + fn)
    precision_score = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall_score = tp / (tp + fn) if (tp + fn) > 0 else 0
    
    summary_data = {
        'Dataset': ['Human (SpliceVarDB)', 'Rat (ratGTEx)', 'Performance Gap'],
        'Species': ['Homo sapiens', 'Rattus norvegicus', '-'],
        'Sample_Size': [len(human_scores), len(rat_scores), '-'],
        'AUROC': [f'{human_auroc:.4f}', f'{rat_auroc:.4f}', f'{human_auroc - rat_auroc:.4f}'],
        'AUPRC': [f'{human_auprc:.4f}', f'{rat_auprc:.4f}', f'{human_auprc - rat_auprc:.4f}'],
        'Rat_Accuracy': ['-', f'{accuracy:.4f}', '-'],
        'Rat_Precision': ['-', f'{precision_score:.4f}', '-'],
        'Rat_Recall': ['-', f'{recall_score:.4f}', '-'],
        'Score_Mean': [f'{human_scores.mean():.6f}', f'{rat_scores.mean():.6f}', '-'],
        'Score_Std': [f'{human_scores.std():.6f}', f'{rat_scores.std():.6f}', '-']
    }
    
    summary_df = pd.DataFrame(summary_data)
    summary_path = os.path.join(output_dir, 'alphagenome_cross_species_summary.tsv')
    summary_df.to_csv(summary_path, sep='\t', index=False)
    
    print(f"Summary table saved to: {summary_path}")
    print("\n=== Summary Table ===")
    print(summary_df.to_string(index=False))
    
    plt.show()

if __name__ == '__main__':
    main()
