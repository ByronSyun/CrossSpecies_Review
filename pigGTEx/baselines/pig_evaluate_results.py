#!/usr/bin/env python3
"""
Evaluate pigGTEx consolidated predictions.
Computes AUROC, AUPRC, and threshold-based metrics with robust thresholds for poorly separable scores.
"""

from pathlib import Path
import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score, precision_recall_curve, accuracy_score, precision_score, recall_score, f1_score


def evaluate_series(y_true: np.ndarray, y_scores: np.ndarray, model_name: str) -> dict:
    mask = ~np.isnan(y_scores)
    y = y_true[mask]
    s = y_scores[mask]
    if len(y) == 0 or len(np.unique(y)) < 2:
        return dict(model=model_name, n=len(y), coverage=0.0, auroc=np.nan, auprc=np.nan,
                    accuracy=np.nan, precision=np.nan, recall=np.nan, f1=np.nan, threshold=np.nan)

    auroc = roc_auc_score(y, s)
    auprc = average_precision_score(y, s)

    # F1-optimal threshold; adjust for poorly separable distributions
    prec, rec, thr = precision_recall_curve(y, s)
    if len(thr) > 0:
        f1_curve = 2 * prec[:-1] * rec[:-1] / (prec[:-1] + rec[:-1] + 1e-12)
        best_idx = int(np.argmax(f1_curve))
        t = float(thr[best_idx])
        # Robust adjustment: if best_idx indicates degenerate prediction, switch to median/percentile
        pos_rate = y.mean()
        frac_pos = float((s >= t).mean())
        if abs(frac_pos - 1.0) < 1e-6 or abs(frac_pos - 0.0) < 1e-6:
            # median fallback
            t = float(np.median(s))
    else:
        t = float(np.median(s))

    y_pred = (s >= t).astype(int)
    acc = accuracy_score(y, y_pred)
    pre = precision_score(y, y_pred, zero_division=0)
    rec_ = recall_score(y, y_pred, zero_division=0)
    f1 = f1_score(y, y_pred, zero_division=0)

    return dict(model=model_name, n=len(y), coverage=len(y)/len(y_true)*100,
                auroc=auroc, auprc=auprc, accuracy=acc, precision=pre, recall=rec_, f1=f1, threshold=t)


def main():
    base = Path('/Users/byronsun/Desktop/AS_复现模型/BIB_review/code/pigGTEx/baselines')
    summary_file = base / 'piggtex_all_model_predictions.csv'
    out_file = base / 'piggtex_model_performance.csv'

    df = pd.read_csv(summary_file)
    y_true = df['ground_truth'].values.astype(int)

    # Evaluate 9 zero-shot models: 4 task-specific + 5 GFMs (MLP models evaluated separately on test sets)
    # Order: Task-specific first, then GFMs
    available_models = ['Pangolin', 'SpliceAI', 'MMSplice_pathogenicity', 'SpliceTransformer',
                        'AlphaGenome', 'Evo2_zeroshot', 'Nucleotide_Transformer', 'SpliceBERT', 'DNABERT2_Logistic']
    models = [c for c in available_models if c in df.columns]
    
    results = []
    for m in models:
        print(f'Evaluating {m} ...')
        res = evaluate_series(y_true, df[m].values.astype(float), m)
        results.append(res)
        print({k: (round(v,4) if isinstance(v, float) else v) for k, v in res.items()})

    perf = pd.DataFrame(results)
    perf_rounded = perf.copy()
    for col in ['coverage','auroc','auprc','accuracy','precision','recall','f1','threshold']:
        if col in perf_rounded.columns:
            perf_rounded[col] = perf_rounded[col].round(4)

    print('\nPigGTEx Performance Summary:')
    print(perf_rounded.to_string(index=False))
    perf_rounded.to_csv(out_file, index=False)
    print(f'\nSaved: {out_file}')


if __name__ == '__main__':
    main()


