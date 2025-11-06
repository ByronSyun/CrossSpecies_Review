#!/bin/bash
set -euo pipefail

CONDA_ENV_PATH="/mnt/userdata4/splicing/conda_envs/dnabert2_env"
ROOT="/mnt/userdata4/splicing"
PYTHON_EXE="$CONDA_ENV_PATH/bin/python"

PROJECT_DIR="$ROOT/DNABert_2_chicken"
OUT_DIR="$PROJECT_DIR/results"
INPUT_TSV="$ROOT/ChickenGTEx/processed_data/chickengtex_silver_benchmark_balanced.tsv"

EMB_NPZ="$OUT_DIR/dnabert2_chicken_concat_embeddings.npz"
MODEL_PATH="$OUT_DIR/dnabert2_chicken_mlp_classifier.joblib"
SCALER_PATH="$OUT_DIR/dnabert2_chicken_mlp_classifier_scaler.joblib"
CONF_IMG="$OUT_DIR/dnabert2_chicken_mlp_confusion_matrix.png"
PROBA_TSV="$OUT_DIR/dnabert2_chicken_mlp_variant_scores.tsv"
METRICS_JSON="$OUT_DIR/dnabert2_chicken_mlp_metrics.json"
TEST_PROBA_TSV="$OUT_DIR/dnabert2_chicken_mlp_variant_scores_test.tsv"

mkdir -p "$PROJECT_DIR" "$OUT_DIR"

echo "=== DNABERT-2 MLP Baseline (Chicken) ==="
echo "Input: $INPUT_TSV"

if [ ! -f "$INPUT_TSV" ]; then
    echo "ERROR: Input not found: $INPUT_TSV"
    exit 1
fi

echo "[1/3] Generating concatenated embeddings..."
if [ -f "$EMB_NPZ" ]; then
    echo "✓ Embeddings exist, skipping"
else
    "$PYTHON_EXE" "$PROJECT_DIR/predict_dnabert2_chicken.py" \
      --input_file "$INPUT_TSV" \
      --output_file "$EMB_NPZ" \
      --batch_size 16
fi

echo "[2/3] Training MLP..."
"$PYTHON_EXE" "$PROJECT_DIR/train_mlp_chicken.py" \
  --input_file "$EMB_NPZ" \
  --model_output_file "$MODEL_PATH" \
  --confusion_matrix_file "$CONF_IMG" \
  --metrics_output_file "$METRICS_JSON" \
  --proba_output_file "$TEST_PROBA_TSV"

echo "[3/3] Generating full predictions..."
"$PYTHON_EXE" - <<PY
import numpy as np
import joblib
from sklearn.preprocessing import StandardScaler

data = np.load("${EMB_NPZ}", allow_pickle=True)
X, y, vids = data['embeddings'], data['labels'], data['variant_ids']

scaler = joblib.load("${SCALER_PATH}")
model = joblib.load("${MODEL_PATH}")

X_scaled = scaler.transform(X)
proba = model.predict_proba(X_scaled)[:, 1]

with open("${PROBA_TSV}", 'w') as f:
    f.write("variant_id\tprob_splice_altering\tlabel\n")
    for v, p, l in zip(vids, proba, y):
        f.write(f"{v}\t{p:.6f}\t{int(l)}\n")

print(f"Saved full predictions: ${PROBA_TSV}")
PY

echo "=== Pipeline completed ==="
echo "Scores TSV: $PROBA_TSV"
echo "Metrics JSON: $METRICS_JSON"

