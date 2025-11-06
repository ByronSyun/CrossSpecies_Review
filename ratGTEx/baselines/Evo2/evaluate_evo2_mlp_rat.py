#!/usr/bin/env python3
import numpy as np
from sklearn.metrics import (accuracy_score, classification_report, 
                              roc_auc_score, average_precision_score, 
                              confusion_matrix, roc_curve, precision_recall_curve)
import joblib
import argparse
import json
import matplotlib.pyplot as plt
import seaborn as sns


def plot_curves(y_true, y_pred_proba, output_prefix):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # ROC curve
    fpr, tpr, _ = roc_curve(y_true, y_pred_proba)
    auroc = roc_auc_score(y_true, y_pred_proba)
    ax1.plot(fpr, tpr, label=f'AUROC = {auroc:.4f}', linewidth=2)
    ax1.plot([0, 1], [0, 1], 'k--', label='Random', linewidth=1)
    ax1.set_xlabel('False Positive Rate', fontsize=12)
    ax1.set_ylabel('True Positive Rate', fontsize=12)
    ax1.set_title('ROC Curve - Evo2 MLP (Rat)', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=11)
    ax1.grid(alpha=0.3)
    
    # PR curve
    precision, recall, _ = precision_recall_curve(y_true, y_pred_proba)
    auprc = average_precision_score(y_true, y_pred_proba)
    baseline = np.sum(y_true) / len(y_true)
    ax2.plot(recall, precision, label=f'AUPRC = {auprc:.4f}', linewidth=2)
    ax2.axhline(y=baseline, color='k', linestyle='--', 
                label=f'Baseline = {baseline:.4f}', linewidth=1)
    ax2.set_xlabel('Recall', fontsize=12)
    ax2.set_ylabel('Precision', fontsize=12)
    ax2.set_title('Precision-Recall Curve - Evo2 MLP (Rat)', 
                  fontsize=14, fontweight='bold')
    ax2.legend(fontsize=11)
    ax2.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_curves.png", dpi=300, bbox_inches='tight')
    plt.close()


def plot_confusion_matrix(cm, output_prefix):
    plt.figure(figsize=(8, 6))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', cbar=True,
                xticklabels=['Normal', 'Splice-Altering'],
                yticklabels=['Normal', 'Splice-Altering'])
    plt.xlabel('Predicted Label', fontsize=12)
    plt.ylabel('True Label', fontsize=12)
    plt.title('Confusion Matrix - Evo2 MLP (Rat)', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_confusion_matrix.png", dpi=300, bbox_inches='tight')
    plt.close()


def main(args):
    data = np.load(args.input_file, allow_pickle=True)
    embeddings = data['embeddings']
    labels = data['labels']
    variant_ids = data['variant_ids']
    
    print(f"Loaded {len(embeddings):,} rat variants, embedding dim: {embeddings.shape[1]:,}D")
    
    model = joblib.load(args.model_file)
    scaler = joblib.load(args.scaler_file)
    
    embeddings_scaled = scaler.transform(embeddings)
    
    y_pred = model.predict(embeddings_scaled)
    y_pred_proba = model.predict_proba(embeddings_scaled)[:, 1]
    
    auroc = roc_auc_score(labels, y_pred_proba)
    auprc = average_precision_score(labels, y_pred_proba)
    accuracy = accuracy_score(labels, y_pred)
    cm = confusion_matrix(labels, y_pred)
    
    print(f"\nAUROC: {auroc:.4f}")
    print(f"AUPRC: {auprc:.4f}")
    print(f"Accuracy: {accuracy:.4f}")
    print(f"\nConfusion Matrix:")
    print(f"  [[TN={cm[0,0]:,}  FP={cm[0,1]:,}]")
    print(f"   [FN={cm[1,0]:,}  TP={cm[1,1]:,}]]")
    print(f"\nClassification Report:")
    print(classification_report(labels, y_pred,
                                target_names=['Normal', 'Splice-Altering']))
    
    metrics = {
        'auroc': float(auroc),
        'auprc': float(auprc),
        'accuracy': float(accuracy),
        'confusion_matrix': cm.tolist(),
        'n_variants': int(len(embeddings))
    }
    
    with open(f"{args.output_dir}/evo2_mlp_rat_metrics.json", 'w') as f:
        json.dump(metrics, f, indent=2)
    
    import pandas as pd
    predictions_df = pd.DataFrame({
        'variant_id': variant_ids,
        'true_label': labels,
        'predicted_label': y_pred,
        'prediction_probability': y_pred_proba
    })
    output_tsv = f"{args.output_dir}/evo2_mlp_rat_predictions.tsv"
    predictions_df.to_csv(output_tsv, sep='\t', index=False)
    
    plot_curves(labels, y_pred_proba, f"{args.output_dir}/evo2_mlp_rat")
    plot_confusion_matrix(cm, f"{args.output_dir}/evo2_mlp_rat")
    
    print(f"\nSaved: {output_tsv}, metrics.json, curves.png, confusion_matrix.png")
    print(f"\nCompare: DNABERT-2 MLP 0.7112 vs Evo2 delta_logp 0.4848")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Evaluate Evo2 MLP on rat data (zero-shot)"
    )
    parser.add_argument("--input_file", type=str, required=True,
                        help="Rat embeddings NPZ file")
    parser.add_argument("--model_file", type=str, required=True,
                        help="Trained MLP model (from human data)")
    parser.add_argument("--scaler_file", type=str, required=True,
                        help="Trained StandardScaler (from human data)")
    parser.add_argument("--output_dir", type=str, required=True,
                        help="Output directory for results")
    
    args = parser.parse_args()
    main(args)

