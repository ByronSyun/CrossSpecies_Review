#!/bin/bash
set -euo pipefail

# Host paths
HOST_CURRENT_PROJECT_DIR="/mnt/userdata4/splicing/Evo2_pig"
HOST_BASE_PROJECT_DIR="/mnt/userdata4/splicing/Evo2Splicevardb"

# Files
EMB_FILE="evo2_pig_embeddings.npz"
RESULTS_DIR="results"
MODEL_FILE="$RESULTS_DIR/evo2_mlp_pig_model.pkl"
SCALER_FILE="$RESULTS_DIR/evo2_mlp_pig_scaler.pkl"
TRAIN_METRICS="$RESULTS_DIR/evo2_mlp_pig_train_test_metrics.json"
EVAL_METRICS="$RESULTS_DIR/evo2_mlp_pig_metrics.json"
PRED_TSV="$RESULTS_DIR/evo2_mlp_pig_predictions.tsv"
PRED_METRICS_NO_PLOT="$RESULTS_DIR/evo2_mlp_pig_metrics_noplot.json"

# Checks
if [ ! -d "$HOST_CURRENT_PROJECT_DIR" ]; then
  echo "ERROR: Project dir not found: $HOST_CURRENT_PROJECT_DIR" >&2
  exit 1
fi

if [ ! -f "$HOST_CURRENT_PROJECT_DIR/$EMB_FILE" ]; then
  echo "ERROR: Embedding file not found: $HOST_CURRENT_PROJECT_DIR/$EMB_FILE" >&2
  exit 1
fi

mkdir -p "$HOST_CURRENT_PROJECT_DIR/$RESULTS_DIR"
# Commands
# Train only (no evaluation; container lacks matplotlib)
TRAIN_CMD="python /app/current_project/train_evo2_mlp.py \
  --input_file /app/current_project/$EMB_FILE \
  --model_file /app/current_project/$MODEL_FILE \
  --scaler_file /app/current_project/$SCALER_FILE \
  --metrics_file /app/current_project/$TRAIN_METRICS"

# Run in container (train only)
apptainer exec \
  --bind "$HOST_CURRENT_PROJECT_DIR:/app/current_project" \
  --bind "$HOST_BASE_PROJECT_DIR:/app/base_project" \
  "$HOST_BASE_PROJECT_DIR/vortex.sif" \
  bash -c "$TRAIN_CMD"

echo "Done. Results in $HOST_CURRENT_PROJECT_DIR/$RESULTS_DIR"
echo "Model:    $HOST_CURRENT_PROJECT_DIR/$MODEL_FILE"
echo "Scaler:   $HOST_CURRENT_PROJECT_DIR/$SCALER_FILE"
echo "Train/Test Metrics:  $HOST_CURRENT_PROJECT_DIR/$TRAIN_METRICS"

# Optional: generate prediction TSV and minimal metrics without matplotlib
apptainer exec \
  --bind "$HOST_CURRENT_PROJECT_DIR:/app/current_project" \
  --bind "$HOST_BASE_PROJECT_DIR:/app/base_project" \
  "$HOST_BASE_PROJECT_DIR/vortex.sif" \
  bash -lc "python - << \"PY\"
import numpy as np, pandas as pd, json, joblib
from sklearn.metrics import roc_auc_score, average_precision_score, accuracy_score, confusion_matrix
emb_path = '/app/current_project/$EMB_FILE'
model_path = '/app/current_project/$MODEL_FILE'
scaler_path = '/app/current_project/$SCALER_FILE'
out_tsv = '/app/current_project/$PRED_TSV'
out_metrics = '/app/current_project/$PRED_METRICS_NO_PLOT'
d = np.load(emb_path, allow_pickle=True)
variant_ids = d['variant_ids']; X = d['embeddings']; y = d['labels']
scaler = joblib.load(scaler_path); model = joblib.load(model_path)
X_scaled = scaler.transform(X)
proba = model.predict_proba(X_scaled)[:, 1]
pred = (proba >= 0.5).astype(int)
df = pd.DataFrame({'variant_id': variant_ids, 'true_label': y, 'prediction_probability': proba, 'predicted_label': pred})
df.to_csv(out_tsv, sep='\t', index=False)
auroc = float(roc_auc_score(y, proba))
auprc = float(average_precision_score(y, proba))
acc = float(accuracy_score(y, pred))
cm = confusion_matrix(y, pred).tolist()
json.dump({'auroc': auroc, 'auprc': auprc, 'accuracy': acc, 'confusion_matrix': cm, 'n_variants': int(len(y))}, open(out_metrics, 'w'), indent=2)
print('Saved:', out_tsv, 'and', out_metrics)
PY"

echo "Prediction TSV:     $HOST_CURRENT_PROJECT_DIR/$PRED_TSV"
echo "Prediction Metrics: $HOST_CURRENT_PROJECT_DIR/$PRED_METRICS_NO_PLOT"

