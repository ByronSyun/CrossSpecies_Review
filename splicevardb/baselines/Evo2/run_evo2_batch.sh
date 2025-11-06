#!/bin/bash
# Evo2 batch job launcher (Apptainer/Singularity)

echo "--- Starting Evo2 Batch Job ---"
DATE=$(date)
echo "Start time: $DATE"

# Host paths
HOST_PROJECT_DIR="/mnt/userdata4/splicing/Evo2Splicevardb"
CONTAINER_PATH="$HOST_PROJECT_DIR/vortex.sif"
CONTAINER_PROJECT_DIR="/project"

# Command executed inside the container
EVO2_COMMAND="cd $CONTAINER_PROJECT_DIR && \
python run_splicevardb_evo2.py \
    --config_path vortex-main/configs/evo2-7b-8k.yml \
    --checkpoint_path model_checkpoint/evo2_7b_base.pt \
    --input_file splicevardb_data/evo2_splicevardb_dataset_dedup.tsv \
    --output_file splicevardb_data/evo2_predictions_evo2_splicevardb_dataset_dedup.tsv \
    --batch_size 2"

# Execute
echo "Binding host '$HOST_PROJECT_DIR' to '$CONTAINER_PROJECT_DIR'"
echo "Executing Evo2..."
apptainer exec \
    --nv \
    --bind "$HOST_PROJECT_DIR:$CONTAINER_PROJECT_DIR" \
    "$CONTAINER_PATH" \
    bash -c "$EVO2_COMMAND"

# Exit status
EXIT_CODE=$?
if [ $EXIT_CODE -eq 0 ]; then
    echo "Command executed successfully."
else
    echo "Error: Exit code $EXIT_CODE."
fi

DATE=$(date)
echo "End time: $DATE"
echo "--- Evo2 Batch Job Finished ---"