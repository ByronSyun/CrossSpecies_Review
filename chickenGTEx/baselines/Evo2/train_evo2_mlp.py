#!/usr/bin/env python3
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.metrics import (accuracy_score, classification_report, 
                              roc_auc_score, average_precision_score, confusion_matrix)
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
import joblib
import argparse
import json


def main(args):
    data = np.load(args.input_file, allow_pickle=True)
    embeddings = data['embeddings']
    labels = data['labels']
    variant_ids = data['variant_ids']
    
    print(f"Loaded {len(embeddings):,} samples, embedding dim: {embeddings.shape[1]:,}D")
    
    X_train, X_test, y_train, y_test, ids_train, ids_test = train_test_split(
        embeddings, labels, variant_ids, 
        test_size=0.2, random_state=42, stratify=labels
    )
    
    print(f"Train: {len(X_train):,}, Test: {len(X_test):,}")
    
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    print("Training MLP [256, 128]...")
    model = MLPClassifier(
        hidden_layer_sizes=(256, 128),
        activation='relu',
        solver='adam',
        alpha=0.0001,
        batch_size=64,
        learning_rate='adaptive',
        learning_rate_init=0.001,
        max_iter=1000,
        random_state=42,
        early_stopping=True,
        validation_fraction=0.1,
        n_iter_no_change=20,
        verbose=True
    )
    
    model.fit(X_train_scaled, y_train)
    
    y_pred = model.predict(X_test_scaled)
    y_pred_proba = model.predict_proba(X_test_scaled)[:, 1]
    
    auroc = roc_auc_score(y_test, y_pred_proba)
    auprc = average_precision_score(y_test, y_pred_proba)
    accuracy = accuracy_score(y_test, y_pred)
    cm = confusion_matrix(y_test, y_pred)
    
    print(f"\nAUROC: {auroc:.4f}")
    print(f"AUPRC: {auprc:.4f}")
    print(f"Accuracy: {accuracy:.4f}")
    print(f"\nConfusion Matrix:")
    print(f"  [[TN={cm[0,0]:,}  FP={cm[0,1]:,}]")
    print(f"   [FN={cm[1,0]:,}  TP={cm[1,1]:,}]]")
    
    print(f"\nClassification Report:")
    print(classification_report(y_test, y_pred, 
                                target_names=['Normal', 'Splice-Altering']))
    
    metrics = {
        'auroc': float(auroc),
        'auprc': float(auprc),
        'accuracy': float(accuracy),
        'confusion_matrix': cm.tolist(),
        'test_size': int(len(X_test)),
        'train_size': int(len(X_train)),
        'embedding_dim': int(embeddings.shape[1])
    }
    
    joblib.dump(model, args.model_file)
    joblib.dump(scaler, args.scaler_file)
    with open(args.metrics_file, 'w') as f:
        json.dump(metrics, f, indent=2)
    
    print(f"\nSaved: {args.model_file}, {args.scaler_file}, {args.metrics_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Train MLP classifier on Evo2 embeddings"
    )
    parser.add_argument("--input_file", type=str, required=True,
                        help="Input NPZ file with embeddings")
    parser.add_argument("--model_file", type=str, required=True,
                        help="Output file for trained MLP model")
    parser.add_argument("--scaler_file", type=str, required=True,
                        help="Output file for StandardScaler")
    parser.add_argument("--metrics_file", type=str, required=True,
                        help="Output JSON file for metrics")
    
    args = parser.parse_args()
    main(args)

