#!/bin/bash
set -euo pipefail

CONDA_ENV_PATH="/mnt/userdata4/splicing/conda_envs/dnabert2_env"
ROOT="/mnt/userdata4/splicing"
PYTHON_EXE="$CONDA_ENV_PATH/bin/python"

PROJECT_DIR="$ROOT/DNABert_2"
OUT_DIR="$PROJECT_DIR/data"
INPUT_TSV="$ROOT/Evo2Splicevardb/splicevardb_data/evo2_splicevardb_dataset_dedup.tsv"

EMB_NPZ="$OUT_DIR/dnabert2_concat_embeddings.npz"
MODEL_PATH="$OUT_DIR/dnabert2_mlp_model.joblib"
SCALER_PATH="$OUT_DIR/dnabert2_mlp_scaler.joblib"
METRICS_JSON="$OUT_DIR/dnabert2_mlp_metrics.json"
PROBA_TSV="$OUT_DIR/dnabert2_mlp_variant_scores.tsv"

mkdir -p "$OUT_DIR"

echo "[1/3] Generating concatenate embeddings [WT, MT, diff]"
"${PYTHON_EXE}" "$PROJECT_DIR/predict_dnabert2_concat.py" \
  --input_file "$INPUT_TSV" \
  --output_file "$EMB_NPZ" \
  --batch_size 16

if [ ! -s "$EMB_NPZ" ]; then
  echo "ERROR: Embeddings not created" >&2
  exit 1
fi

echo "[2/3] Training MLP classifier"
"${PYTHON_EXE}" "$PROJECT_DIR/train_mlp_concat.py" \
  --input_file "$EMB_NPZ" \
  --model_file "$MODEL_PATH" \
  --scaler_file "$SCALER_PATH" \
  --metrics_file "$METRICS_JSON"

if [ ! -s "$MODEL_PATH" ]; then
  echo "ERROR: Model not created" >&2
  exit 1
fi

echo "[3/3] Exporting variant probabilities"
"${PYTHON_EXE}" - <<PY
import numpy as np, joblib
from pathlib import Path

data = np.load("${EMB_NPZ}", allow_pickle=True)
X = data['embeddings']
y = data['labels']
vids = data['variant_ids']

scaler = joblib.load("${SCALER_PATH}")
model = joblib.load("${MODEL_PATH}")

X_scaled = scaler.transform(X)
proba = model.predict_proba(X_scaled)[:, 1]

with open("${PROBA_TSV}", 'w') as f:
    f.write("variant_id\tprob_splice_altering\tlabel\n")
    for v, p, l in zip(vids, proba, y):
        f.write(f"{v}\t{p:.6f}\t{int(l)}\n")
print(f"Saved: ${PROBA_TSV}")
PY

if [ ! -s "$PROBA_TSV" ]; then
  echo "ERROR: Probability TSV not created" >&2
  exit 1
fi

echo "Completed. Results:"
echo "  Embeddings: $EMB_NPZ"
echo "  Model:      $MODEL_PATH"
echo "  Scaler:     $SCALER_PATH"
echo "  Metrics:    $METRICS_JSON"
echo "  Scores:     $PROBA_TSV"

