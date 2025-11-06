#!/bin/bash
# Nucleotide Transformer Pipeline for SpliceVarDB and RatGTEx
# Run this script on the server after environment setup

set -e  # Exit on any error

# Configuration
ENV_PATH="/mnt/userdata4/splicing/conda_envs/nt-repo"
RESULTS_DIR="/mnt/userdata4/splicing/Nucleotide Transformer/results"
MODEL_ID="InstaDeepAI/nucleotide-transformer-v2-500m-multi-species"

# Data paths
SPLICEVARDB_DATA="/mnt/userdata4/splicing/Evo2Splicevardb/splicevardb_data/evo2_splicevardb_dataset_dedup.tsv"
RATGTEX_DATA="/mnt/userdata4/splicing/ratgetx/processed_data/ratgtex_silver_benchmark_balanced.tsv"

# Create results directory
mkdir -p "$RESULTS_DIR"

echo "=== Nucleotide Transformer Pipeline ==="
echo "Environment: $ENV_PATH"
echo "Results: $RESULTS_DIR"
echo "Model: $MODEL_ID"
echo ""

# Set HuggingFace cache
export HF_HOME="/mnt/userdata4/splicing/.cache/huggingface"
export TRANSFORMERS_CACHE="/mnt/userdata4/splicing/.cache/huggingface"

# 1. Human SpliceVarDB
echo "1. Running NT on Human SpliceVarDB..."
if [ -f "$SPLICEVARDB_DATA" ]; then
    conda run -p "$ENV_PATH" python nt_score_variants.py \
        --model_id "$MODEL_ID" \
        --input_tsv "$SPLICEVARDB_DATA" \
        --output_prefix "$RESULTS_DIR/nt_v2_splicevardb" \
        --window 8192 \
        --batch_size 2 \
        --num_workers 1
    
    echo "1a. Evaluating Human results..."
    conda run -p "$ENV_PATH" python evaluate_nt_scores.py \
        --labels "$SPLICEVARDB_DATA" \
        --scores "$RESULTS_DIR/nt_v2_splicevardb_scores.csv" \
        --out_dir "$RESULTS_DIR/splicevardb_eval"
else
    echo "Warning: SpliceVarDB data not found at $SPLICEVARDB_DATA"
fi

echo ""

# 2. Rat RatGTEx
echo "2. Running NT on Rat RatGTEx..."
if [ -f "$RATGTEX_DATA" ]; then
    conda run -p "$ENV_PATH" python nt_score_variants.py \
        --model_id "$MODEL_ID" \
        --input_tsv "$RATGTEX_DATA" \
        --output_prefix "$RESULTS_DIR/nt_v2_ratgtex" \
        --window 8192 \
        --batch_size 2 \
        --num_workers 1
    
    echo "2a. Evaluating Rat results..."
    conda run -p "$ENV_PATH" python evaluate_nt_scores.py \
        --labels "$RATGTEX_DATA" \
        --scores "$RESULTS_DIR/nt_v2_ratgtex_scores.csv" \
        --out_dir "$RESULTS_DIR/ratgtex_eval"
else
    echo "Warning: RatGTEx data not found at $RATGTEX_DATA"
fi

echo ""
echo "=== Pipeline Complete ==="
echo "Results saved to: $RESULTS_DIR"
echo "Check the following files:"
echo "- $RESULTS_DIR/splicevardb_eval/nt_evaluation_results.json"
echo "- $RESULTS_DIR/ratgtex_eval/nt_evaluation_results.json"
