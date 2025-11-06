#!/usr/bin/env python3
"""Train Logistic Regression classifier on DNABERT-2 differential embeddings for chicken data."""
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, classification_report, roc_auc_score, average_precision_score, confusion_matrix
import joblib
import argparse
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def train_and_evaluate(args):
    print(f"Loading data from {args.input_file}...")
    data = np.load(args.input_file, allow_pickle=True)
    embeddings = data['embeddings']
    labels = data['labels']
    variant_ids = data['variant_ids']
    
    print(f"Loaded {len(embeddings)} samples, embedding dim: {embeddings.shape[1]}")
    
    X_train, X_test, y_train, y_test, ids_train, ids_test = train_test_split(
        embeddings, labels, variant_ids, test_size=0.2, random_state=42, stratify=labels
    )
    
    print(f"Train: {len(X_train)}, Test: {len(X_test)}")
    
    print("Training Logistic Regression...")
    model = LogisticRegression(random_state=42, max_iter=1000, class_weight='balanced')
    model.fit(X_train, y_train)
    
    print("Evaluating...")
    y_pred = model.predict(X_test)
    y_pred_proba = model.predict_proba(X_test)[:, 1]
    
    accuracy = accuracy_score(y_test, y_pred)
    roc_auc = roc_auc_score(y_test, y_pred_proba)
    auprc = average_precision_score(y_test, y_pred_proba)
    
    print(f"\nResults:")
    print(f"Accuracy: {accuracy:.4f}, AUROC: {roc_auc:.4f}, AUPRC: {auprc:.4f}")
    print(f"\n{classification_report(y_test, y_pred)}")
    
    cm = confusion_matrix(y_test, y_pred)
    print(f"Confusion Matrix:\n{cm}")
    
    plt.figure(figsize=(8, 6))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues')
    plt.title(f'DNABERT-2 Logistic Chicken\nAUROC={roc_auc:.4f}')
    plt.ylabel('True')
    plt.xlabel('Predicted')
    plt.savefig(args.confusion_matrix_file, dpi=300, bbox_inches='tight')
    print(f"Saved confusion matrix: {args.confusion_matrix_file}")
    
    joblib.dump(model, args.model_output_file)
    print(f"Saved model: {args.model_output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train Logistic Regression on DNABERT-2 embeddings for chicken")
    parser.add_argument("--input_file", type=str, required=True, help="Input NPZ file")
    parser.add_argument("--model_output_file", type=str, default="dnabert2_chicken_logistic_classifier.joblib", help="Output model file")
    parser.add_argument("--confusion_matrix_file", type=str, default="dnabert2_chicken_logistic_confusion_matrix.png", help="Output confusion matrix")
    args = parser.parse_args()
    train_and_evaluate(args)

