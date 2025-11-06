#!/bin/bash
# Nucleotide Transformer Pipeline for PigGTEx
# Simplified version: direct command execution

set -euo pipefail

echo "=========================================="
echo "Running Nucleotide Transformer on PigGTEx"
echo "=========================================="

# Activate conda environment
conda run -p /mnt/userdata4/splicing/conda_envs/nt-repo \
  python -u nt_score_variants.py \
    --model_id "InstaDeepAI/nucleotide-transformer-500m-1000g" \
    --input_tsv piggtex/piggtex_with_header.tsv \
    --output_prefix results/nt_pig \
    --window 8192 \
    --batch_size 16 \
    --num_workers 2

echo ""
echo "=========================================="
echo "NT Pipeline Completed!"
echo "=========================================="
