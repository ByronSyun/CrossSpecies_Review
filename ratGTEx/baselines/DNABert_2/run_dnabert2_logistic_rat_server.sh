#!/bin/bash
# DNABERT-2 Logistic Regression baseline for rat data (server script)
# 1) Generate DNABERT-2 differential embeddings (ALT - REF)
# 2) Train/evaluate Logistic Regression classifier
# 3) Export per-variant probabilities for benchmarking
#
# Run on server at: /mnt/userdata4/splicing/DNABert_2_rat/

set -euo pipefail

# --- Configuration (server paths) ---
CONDA_ENV_PATH="/mnt/userdata4/splicing/conda_envs/dnabert2_env"
ROOT="/mnt/userdata4/splicing"
PYTHON_EXE="$CONDA_ENV_PATH/bin/python"

# Project directory on server
PROJECT_DIR="$ROOT/DNABert_2_rat"
OUT_DIR="$PROJECT_DIR/results"

# Input data (rat benchmark)
INPUT_TSV="$ROOT/ratgetx/processed_data/ratgtex_silver_benchmark_balanced.tsv"

# Outputs (Logistic-specific naming to distinguish from MLP method)
EMB_NPZ="$OUT_DIR/dnabert2_rat_logistic_embeddings.npz"
MODEL_PATH="$OUT_DIR/dnabert2_rat_logistic_classifier.joblib"
CONF_IMG="$OUT_DIR/dnabert2_rat_logistic_confusion_matrix.png"
PROBA_TSV="$OUT_DIR/dnabert2_rat_logistic_scores.tsv"

mkdir -p "$OUT_DIR"

echo "=== DNABERT-2 Logistic Regression Baseline (Rat) ==="
echo "Input: $INPUT_TSV"
echo "Output dir: $OUT_DIR"

# Check input file
if [ ! -f "$INPUT_TSV" ]; then
  echo "ERROR: Input file not found: $INPUT_TSV" >&2
  exit 1
fi

echo ""
echo "[1/3] Generating DNABERT-2 differential embeddings (ALT - REF) -> $EMB_NPZ"
"${PYTHON_EXE:-python}" "$PROJECT_DIR/predict_dnabert2_logistic_rat.py" \
  --input_file "$INPUT_TSV" \
  --output_file "$EMB_NPZ" \
  --batch_size 16

if [ ! -s "$EMB_NPZ" ]; then
  echo "ERROR: Embeddings NPZ not created: $EMB_NPZ" >&2
  exit 1
fi

echo ""
echo "[2/3] Training and evaluating Logistic Regression classifier -> $MODEL_PATH"
"${PYTHON_EXE:-python}" "$PROJECT_DIR/train_classifier_rat.py" \
  --input_file "$EMB_NPZ" \
  --model_output_file "$MODEL_PATH" \
  --confusion_matrix_file "$CONF_IMG"

if [ ! -s "$MODEL_PATH" ]; then
  echo "ERROR: Model not created: $MODEL_PATH" >&2
  exit 1
fi

echo ""
echo "[3/3] Exporting per-variant probabilities -> $PROBA_TSV"
"${PYTHON_EXE:-python}" - <<PY
import numpy as np
import joblib
import sys
from pathlib import Path

emb_npz = Path("${EMB_NPZ}")
model_path = Path("${MODEL_PATH}")
out_tsv = Path("${PROBA_TSV}")

print("Loading embeddings and model...")
data = np.load(emb_npz, allow_pickle=True)
X = data['embeddings']
y = data['labels']
vids = data['variant_ids']

model = joblib.load(model_path)

print(f"Generating predictions for {len(vids)} variants...")
proba = model.predict_proba(X)[:, 1]

print(f"Writing to {out_tsv}...")
with out_tsv.open('w') as f:
    f.write("variant_id\tprob_splice_altering\tlabel\n")
    for v, p, l in zip(vids, proba, y):
        f.write(f"{v}\t{p:.6f}\t{int(l)}\n")

print(f"Saved: {out_tsv}")
print(f"Total variants: {len(vids)}")
PY

if [ ! -s "$PROBA_TSV" ]; then
  echo "ERROR: Probability TSV not created: $PROBA_TSV" >&2
  exit 1
fi

echo ""
echo "=== Pipeline completed successfully ==="
echo "Outputs:"
echo " - Embeddings:  $EMB_NPZ"
echo " - Model:       $MODEL_PATH"
echo " - Confusion:   $CONF_IMG"
echo " - Scores TSV:  $PROBA_TSV"
echo ""
echo "Use this TSV for benchmarking: $PROBA_TSV"

