#!/bin/bash
set -euo pipefail

echo "Start time: $(date)"

HOST_CURRENT_PROJECT_DIR="/mnt/userdata4/splicing/Evo2_chicken"
HOST_BASE_PROJECT_DIR="/mnt/userdata4/splicing/Evo2Splicevardb"
HOST_INPUT_FILE="/mnt/userdata4/splicing/ChickenGTEx/processed_data/chickengtex_silver_benchmark_balanced.tsv"
OUTPUT_FILE="evo2_chicken_predictions.tsv"
BATCH_SIZE=1

if [ ! -d "$HOST_CURRENT_PROJECT_DIR" ]; then
    echo "ERROR: Current project directory not found: $HOST_CURRENT_PROJECT_DIR"
    exit 1
fi

if [ ! -d "$HOST_BASE_PROJECT_DIR" ]; then
    echo "ERROR: Base project directory not found: $HOST_BASE_PROJECT_DIR"
    exit 1
fi

if [ ! -f "$HOST_INPUT_FILE" ]; then
    echo "ERROR: Input file not found: $HOST_INPUT_FILE"
    exit 1
fi

if [ ! -f "$HOST_CURRENT_PROJECT_DIR/predict_evo2_chicken.py" ]; then
    echo "ERROR: Prediction script not found: $HOST_CURRENT_PROJECT_DIR/predict_evo2_chicken.py"
    exit 1
fi

HOST_DATA_DIR="/mnt/userdata4/splicing/ChickenGTEx/processed_data"

EVO2_COMMAND="python /app/current_project/predict_evo2_chicken.py \
    --config_path /app/base_project/vortex-main/configs/evo2-7b-8k.yml \
    --checkpoint_path /app/base_project/model_checkpoint/evo2_7b_base.pt \
    --input_file /app/data/chickengtex_silver_benchmark_balanced.tsv \
    --output_file /app/current_project/$OUTPUT_FILE \
    --batch_size $BATCH_SIZE"

apptainer exec \
    --nv \
    --bind "$HOST_CURRENT_PROJECT_DIR:/app/current_project" \
    --bind "$HOST_BASE_PROJECT_DIR:/app/base_project" \
    --bind "$HOST_DATA_DIR:/app/data" \
    "$HOST_BASE_PROJECT_DIR/vortex.sif" \
    bash -c "$EVO2_COMMAND"

EXIT_CODE=$?

if [ $EXIT_CODE -eq 0 ]; then
    echo "Evo2 prediction completed successfully!"
    echo "Output: $HOST_CURRENT_PROJECT_DIR/$OUTPUT_FILE"
    if [ -f "$HOST_CURRENT_PROJECT_DIR/$OUTPUT_FILE" ]; then
        echo "Lines: $(wc -l < "$HOST_CURRENT_PROJECT_DIR/$OUTPUT_FILE")"
        head -5 "$HOST_CURRENT_PROJECT_DIR/$OUTPUT_FILE"
    fi
else
    echo "ERROR: Evo2 prediction failed with exit code $EXIT_CODE"
fi

echo "End time: $(date)"

