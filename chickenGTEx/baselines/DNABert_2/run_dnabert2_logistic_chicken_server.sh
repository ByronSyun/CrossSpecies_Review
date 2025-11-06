#!/bin/bash
set -euo pipefail

CONDA_ENV_PATH="/mnt/userdata4/splicing/conda_envs/dnabert2_env"
ROOT="/mnt/userdata4/splicing"
PYTHON_EXE="$CONDA_ENV_PATH/bin/python"

PROJECT_DIR="$ROOT/DNABert_2_chicken"
OUT_DIR="$PROJECT_DIR/results"
INPUT_TSV="$ROOT/ChickenGTEx/processed_data/chickengtex_silver_benchmark_balanced.tsv"

EMB_NPZ="$OUT_DIR/dnabert2_chicken_logistic_embeddings.npz"
MODEL_PATH="$OUT_DIR/dnabert2_chicken_logistic_classifier.joblib"
CONF_IMG="$OUT_DIR/dnabert2_chicken_logistic_confusion_matrix.png"
PROBA_TSV="$OUT_DIR/dnabert2_chicken_logistic_scores.tsv"

mkdir -p "$OUT_DIR"

echo "=== DNABERT-2 Logistic Regression Baseline (Chicken) ==="
echo "Input: $INPUT_TSV"

if [ ! -f "$INPUT_TSV" ]; then
  echo "ERROR: Input file not found: $INPUT_TSV" >&2
  exit 1
fi

echo "[1/3] Generating embeddings..."
"$PYTHON_EXE" "$PROJECT_DIR/predict_dnabert2_logistic_chicken.py" \
  --input_file "$INPUT_TSV" \
  --output_file "$EMB_NPZ" \
  --batch_size 16

echo "[2/3] Training classifier..."
"$PYTHON_EXE" "$PROJECT_DIR/train_classifier_chicken.py" \
  --input_file "$EMB_NPZ" \
  --model_output_file "$MODEL_PATH" \
  --confusion_matrix_file "$CONF_IMG"

echo "[3/3] Exporting probabilities..."
"$PYTHON_EXE" - <<PY
import numpy as np
import joblib

data = np.load("${EMB_NPZ}", allow_pickle=True)
X, y, vids = data['embeddings'], data['labels'], data['variant_ids']

model = joblib.load("${MODEL_PATH}")
proba = model.predict_proba(X)[:, 1]

with open("${PROBA_TSV}", 'w') as f:
    f.write("variant_id\tprob_splice_altering\tlabel\n")
    for v, p, l in zip(vids, proba, y):
        f.write(f"{v}\t{p:.6f}\t{int(l)}\n")

print(f"Saved: ${PROBA_TSV}")
PY

echo "=== Pipeline completed ==="
echo "Scores TSV: $PROBA_TSV"

