#!/bin/bash
set -euo pipefail

echo "=== Evo2 Concatenate Embeddings + MLP Pipeline ==="
DATE=$(date)
echo "Start time: $DATE"

# Host paths (adjust to your server setup)
HOST_PROJECT_DIR="/mnt/userdata4/splicing/Evo2Splicevardb"
CONTAINER_PATH="$HOST_PROJECT_DIR/vortex.sif"
CONTAINER_PROJECT_DIR="/project"

# Files
INPUT_TSV="splicevardb_data/evo2_splicevardb_dataset_dedup.tsv"
EMB_OUTPUT="splicevardb_data/evo2_concat_embeddings.npz"
MODEL_FILE="splicevardb_data/evo2_mlp_model.pkl"
SCALER_FILE="splicevardb_data/evo2_mlp_scaler.pkl"
TRAIN_METRICS="splicevardb_data/evo2_mlp_train_test_metrics.json"
PREDICTIONS_FILE="splicevardb_data/evo2_mlp_variant_scores.tsv"
PRED_METRICS="splicevardb_data/evo2_mlp_full_metrics.json"

# Step 1: Extract embeddings
echo ""
echo "--- Step 1: Extracting Evo2 Concatenated Embeddings ---"
EMB_COMMAND="cd $CONTAINER_PROJECT_DIR && \
python extract_evo2_embeddings.py \
    --config_path vortex-main/configs/evo2-7b-8k.yml \
    --checkpoint_path model_checkpoint/evo2_7b_base.pt \
    --input_file $INPUT_TSV \
    --output_file $EMB_OUTPUT \
    --batch_size 1"

apptainer exec \
    --nv \
    --bind "$HOST_PROJECT_DIR:$CONTAINER_PROJECT_DIR" \
    "$CONTAINER_PATH" \
    bash -c "$EMB_COMMAND"

EMB_EXIT=$?
if [ $EMB_EXIT -ne 0 ]; then
    echo "ERROR: Embedding extraction failed with exit code $EMB_EXIT"
    exit $EMB_EXIT
fi

echo "Embedding extraction completed successfully."
echo ""

# Step 2: Train MLP
echo "--- Step 2: Training MLP Classifier ---"
MLP_COMMAND="cd $CONTAINER_PROJECT_DIR && \
python train_evo2_mlp.py \
    --input_file $EMB_OUTPUT \
    --model_file $MODEL_FILE \
    --scaler_file $SCALER_FILE \
    --metrics_file $TRAIN_METRICS"

apptainer exec \
    --nv \
    --bind "$HOST_PROJECT_DIR:$CONTAINER_PROJECT_DIR" \
    "$CONTAINER_PATH" \
    bash -c "$MLP_COMMAND"

MLP_EXIT=$?
if [ $MLP_EXIT -ne 0 ]; then
    echo "ERROR: MLP training failed with exit code $MLP_EXIT"
    exit $MLP_EXIT
fi

echo "MLP training completed successfully."
echo ""

# Step 3: Generate predictions for all variants
echo "--- Step 3: Generating Predictions for All Variants ---"
apptainer exec \
    --nv \
    --bind "$HOST_PROJECT_DIR:$CONTAINER_PROJECT_DIR" \
    "$CONTAINER_PATH" \
    bash -c "cd $CONTAINER_PROJECT_DIR && python - << \"PYSCRIPT\"
import numpy as np
import pandas as pd
import json
import joblib
from sklearn.metrics import roc_auc_score, average_precision_score, accuracy_score, confusion_matrix

# Load data
emb_path = '$EMB_OUTPUT'
model_path = '$MODEL_FILE'
scaler_path = '$SCALER_FILE'
out_tsv = '$PREDICTIONS_FILE'
out_metrics = '$PRED_METRICS'

print('Loading embeddings...')
data = np.load(emb_path, allow_pickle=True)
variant_ids = data['variant_ids']
X = data['embeddings']
y = data['labels']

print('Loading model and scaler...')
scaler = joblib.load(scaler_path)
model = joblib.load(model_path)

print('Generating predictions...')
X_scaled = scaler.transform(X)
proba = model.predict_proba(X_scaled)[:, 1]
pred = (proba >= 0.5).astype(int)

print('Saving predictions TSV...')
df = pd.DataFrame({
    'variant_id': variant_ids,
    'true_label': y,
    'prediction_probability': proba,
    'predicted_label': pred
})
df.to_csv(out_tsv, sep='\t', index=False)

print('Calculating metrics...')
auroc = float(roc_auc_score(y, proba))
auprc = float(average_precision_score(y, proba))
acc = float(accuracy_score(y, pred))
cm = confusion_matrix(y, pred).tolist()

metrics = {
    'auroc': auroc,
    'auprc': auprc,
    'accuracy': acc,
    'confusion_matrix': cm,
    'n_variants': int(len(y))
}

with open(out_metrics, 'w') as f:
    json.dump(metrics, f, indent=2)

print(f'Saved predictions: {out_tsv}')
print(f'Saved metrics: {out_metrics}')
print(f'AUROC: {auroc:.4f}, AUPRC: {auprc:.4f}, Accuracy: {acc:.4f}')
PYSCRIPT"

PRED_EXIT=$?
if [ $PRED_EXIT -ne 0 ]; then
    echo "ERROR: Prediction generation failed with exit code $PRED_EXIT"
    exit $PRED_EXIT
fi

echo "Prediction generation completed successfully."
echo ""

# Summary
echo "=== Pipeline Completed Successfully ==="
echo "Results saved to:"
echo "  - Embeddings: $HOST_PROJECT_DIR/$EMB_OUTPUT"
echo "  - Model: $HOST_PROJECT_DIR/$MODEL_FILE"
echo "  - Scaler: $HOST_PROJECT_DIR/$SCALER_FILE"
echo "  - Train/Test Metrics: $HOST_PROJECT_DIR/$TRAIN_METRICS"
echo "  - Full Predictions TSV: $HOST_PROJECT_DIR/$PREDICTIONS_FILE"
echo "  - Full Metrics: $HOST_PROJECT_DIR/$PRED_METRICS"

DATE=$(date)
echo "End time: $DATE"

