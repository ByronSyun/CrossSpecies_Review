#!/bin/bash
set -euo pipefail

echo "Start: $(date)"

HOST_CURRENT_PROJECT_DIR="/mnt/userdata4/splicing/Evo2_pig"
HOST_BASE_PROJECT_DIR="/mnt/userdata4/splicing/Evo2Splicevardb"
HOST_DATA_DIR="/mnt/userdata4/splicing/PigGTEx/processed_data"
HOST_INPUT_FILE="/mnt/userdata4/splicing/PigGTEx/processed_data/piggtex_silver_benchmark_balanced.tsv"
OUTPUT_FILE="evo2_pig_embeddings.npz"
BATCH_SIZE=1

if [ ! -d "$HOST_CURRENT_PROJECT_DIR" ]; then
    echo "ERROR: Project directory not found: $HOST_CURRENT_PROJECT_DIR"
    exit 1
fi

if [ ! -f "$HOST_INPUT_FILE" ]; then
    echo "ERROR: Input file not found: $HOST_INPUT_FILE"
    exit 1
fi

if [ ! -f "$HOST_CURRENT_PROJECT_DIR/extract_evo2_embeddings.py" ]; then
    echo "ERROR: Script not found: extract_evo2_embeddings.py"
    exit 1
fi

EVO2_COMMAND="python /app/current_project/extract_evo2_embeddings.py \
    --config_path /app/base_project/vortex-main/configs/evo2-7b-8k.yml \
    --checkpoint_path /app/base_project/model_checkpoint/evo2_7b_base.pt \
    --input_file /app/data/piggtex_silver_benchmark_balanced.tsv \
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
    echo "Success! Output: $HOST_CURRENT_PROJECT_DIR/$OUTPUT_FILE"
    if [ -f "$HOST_CURRENT_PROJECT_DIR/$OUTPUT_FILE" ]; then
        echo "File size: $(du -h "$HOST_CURRENT_PROJECT_DIR/$OUTPUT_FILE" | cut -f1)"
    fi
else
    echo "ERROR: Extraction failed with exit code $EXIT_CODE"
fi

echo "End: $(date)"
exit $EXIT_CODE

