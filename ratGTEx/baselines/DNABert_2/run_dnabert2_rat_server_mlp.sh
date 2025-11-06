#!/bin/bash
set -euo pipefail

echo "Start: $(date)"

CONDA_ENV_PATH="/mnt/userdata4/splicing/conda_envs/dnabert2_env"
ROOT="/mnt/userdata4/splicing"
PYTHON_EXE="$CONDA_ENV_PATH/bin/python"

PROJECT_DIR="$ROOT/DNABert_2_rat"
OUT_DIR="$PROJECT_DIR/results"
INPUT_TSV="$ROOT/ratgetx/processed_data/ratgtex_silver_benchmark_balanced.tsv"

EMB_NPZ="$OUT_DIR/dnabert2_rat_concat_embeddings.npz"
MODEL_PATH="$OUT_DIR/dnabert2_rat_mlp_classifier.joblib"
CONF_IMG="$OUT_DIR/dnabert2_rat_mlp_confusion_matrix.png"
PROBA_TSV="$OUT_DIR/dnabert2_rat_mlp_variant_scores.tsv"
METRICS_JSON="$OUT_DIR/dnabert2_rat_mlp_metrics.json"

mkdir -p "$PROJECT_DIR"
mkdir -p "$OUT_DIR"

if [ ! -f "$INPUT_TSV" ]; then
    echo "ERROR: Input not found: $INPUT_TSV"
    exit 1
fi

if [ ! -d "$CONDA_ENV_PATH" ]; then
    echo "ERROR: Conda env not found: $CONDA_ENV_PATH"
    exit 1
fi

echo "[1/2] Generating embeddings..."
if [ -f "$EMB_NPZ" ]; then
    echo "✓ Embeddings already exist, skipping generation"
else
    "${PYTHON_EXE}" "$PROJECT_DIR/predict_dnabert2_rat.py" \
      --input_file "$INPUT_TSV" \
      --output_file "$EMB_NPZ" \
      --batch_size 16
    
    if [ ! -f "$EMB_NPZ" ]; then
        echo "ERROR: Embedding generation failed"
        exit 1
    fi
fi

echo "[2/3] Training MLP..."
SCALER_FILE="${MODEL_PATH%.joblib}_scaler.joblib"
TEST_PROBA_TSV="$OUT_DIR/dnabert2_rat_mlp_variant_scores_test.tsv"
"${PYTHON_EXE}" "$PROJECT_DIR/train_mlp_rat.py" \
  --input_file "$EMB_NPZ" \
  --model_output_file "$MODEL_PATH" \
  --confusion_matrix_file "$CONF_IMG" \
  --metrics_output_file "$METRICS_JSON" \
  --proba_output_file "$TEST_PROBA_TSV"

if [ ! -f "$METRICS_JSON" ]; then
    echo "ERROR: Training failed"
    exit 1
fi

echo "[3/3] Predicting on all variants using trained model..."
"${PYTHON_EXE}" - <<PY
import numpy as np
import joblib
import sys

try:
    data = np.load("${EMB_NPZ}", allow_pickle=True)
    X = data['embeddings']
    y = data['labels']
    vids = data['variant_ids']
    
    scaler = joblib.load("${SCALER_FILE}")
    model = joblib.load("${MODEL_PATH}")
    
    X_scaled = scaler.transform(X)
    proba = model.predict_proba(X_scaled)[:, 1]
    
    with open("${PROBA_TSV}", 'w') as f:
        f.write("variant_id\tprob_splice_altering\tlabel\n")
        for v, p, l in zip(vids, proba, y):
            f.write(f"{v}\t{p:.6f}\t{int(l)}\n")
    print(f"Saved all variant predictions: ${PROBA_TSV} ({len(vids)} variants)")
except Exception as e:
    print(f"ERROR: {e}", file=sys.stderr)
    sys.exit(1)
PY

if [ ! -f "$PROBA_TSV" ]; then
    echo "ERROR: Prediction failed"
    exit 1
fi

echo ""
echo "=== DNABERT-2 MLP Rat Pipeline Completed ==="
echo "Embeddings:  $EMB_NPZ"
echo "Model:       $MODEL_PATH"
echo "Scaler:      $SCALER_FILE"
echo "Confusion:   $CONF_IMG"
echo "Metrics:     $METRICS_JSON"
echo "Test Scores: $TEST_PROBA_TSV"
echo "All Scores:  $PROBA_TSV"
echo "End: $(date)"

if [ -f "$METRICS_JSON" ]; then
    echo ""
    cat "$METRICS_JSON"
fi
