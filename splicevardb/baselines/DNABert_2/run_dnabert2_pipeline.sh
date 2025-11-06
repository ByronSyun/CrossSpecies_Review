#!/bin/bash
# DNABERT-2 baseline pipeline (server):
# 1) Generate DNABERT-2 differential embeddings (ALT - REF)
# 2) Train/evaluate a simple classifier (Logistic Regression)
# 3) Export per-variant probabilities for benchmarking

set -euo pipefail

# --- Configuration (absolute paths on server, mirror SpliceAI style) ---
CONDA_ENV_PATH="/mnt/userdata4/splicing/conda_envs/dnabert2_env"
ROOT="/mnt/userdata4/splicing"
PYTHON_EXE="$CONDA_ENV_PATH/bin/python"

# Code directory on server (not local workspace path)
PROJECT_DIR="$ROOT/DNABert_2"
OUT_DIR="$PROJECT_DIR/data"

# Input benchmark (same as SpliceAI script source)
INPUT_TSV="$ROOT/Evo2Splicevardb/splicevardb_data/evo2_splicevardb_dataset_dedup.tsv"

# Outputs
EMB_NPZ="$OUT_DIR/dnabert2_predictions.npz"
MODEL_PATH="$OUT_DIR/dnabert2_classifier.joblib"
CONF_IMG="$OUT_DIR/confusion_matrix.png"
PROBA_TSV="$OUT_DIR/dnabert2_variant_scores.tsv"

mkdir -p "$OUT_DIR"

echo "[1/3] Generating DNABERT-2 differential embeddings -> $EMB_NPZ"
"${PYTHON_EXE:-python}" "$PROJECT_DIR/predict_dnabert2.py" \
  --input_file "$INPUT_TSV" \
  --output_file "$EMB_NPZ" \
  --batch_size 16

if [ ! -s "$EMB_NPZ" ]; then
  echo "ERROR: Embeddings NPZ not created: $EMB_NPZ" >&2
  exit 1
fi

echo "[2/3] Training and evaluating classifier -> $MODEL_PATH"
"${PYTHON_EXE:-python}" "$PROJECT_DIR/train_classifier.py" \
  --input_file "$EMB_NPZ" \
  --model_output_file "$MODEL_PATH" \
  --confusion_matrix_file "$CONF_IMG"

if [ ! -s "$MODEL_PATH" ]; then
  echo "ERROR: Model not created: $MODEL_PATH" >&2
  exit 1
fi

echo "[3/3] Exporting per-variant probabilities -> $PROBA_TSV"
"${PYTHON_EXE:-python}" - <<PY
import numpy as np, joblib, sys
from pathlib import Path

emb_npz = Path("${EMB_NPZ}")
model_path = Path("${MODEL_PATH}")
out_tsv = Path("${PROBA_TSV}")

data = np.load(emb_npz, allow_pickle=True)
X = data['embeddings']
y = data['labels']
vids = data['variant_ids']

model = joblib.load(model_path)
proba = model.predict_proba(X)[:, 1]

with out_tsv.open('w') as f:
    f.write("variant_id\tprob_splice_altering\tlabel\n")
    for v, p, l in zip(vids, proba, y):
        f.write(f"{v}\t{p:.6f}\t{int(l)}\n")
print(f"Saved: {out_tsv}")
PY

if [ ! -s "$PROBA_TSV" ]; then
  echo "ERROR: Probability TSV not created: $PROBA_TSV" >&2
  exit 1
fi

echo "Pipeline completed. Outputs:"
echo " - Embeddings: $EMB_NPZ"
echo " - Model:      $MODEL_PATH"
echo " - Confusion:  $CONF_IMG"
echo " - Prob TSV:   $PROBA_TSV"


