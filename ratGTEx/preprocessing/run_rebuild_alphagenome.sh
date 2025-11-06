#!/bin/bash
set -e

# --- Configuration ---
PYTHON_EXECUTABLE="/mnt/userdata4/splicing/conda_envs/sptransformer-env/bin/python"
SAMTOOLS_EXECUTABLE="/mnt/userdata4/splicing/conda_envs/sptransformer-env/bin/samtools"
echo "Using Python executable: ${PYTHON_EXECUTABLE}"

# Resolve project root as the parent of this script's directory (…/splicing)
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# Paths
# FIX: reference genome is under ratgetx/reference_genome/rn6
REF_GENOME_DIR="${PROJECT_ROOT}/ratgetx/reference_genome/rn6"
FASTA_FILE="${REF_GENOME_DIR}/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa"

INPUT_TSV="${PROJECT_ROOT}/ratgetx/processed_data/ratgtex_silver_benchmark_balanced.tsv"
OUTPUT_TSV="${PROJECT_ROOT}/ratgetx/processed_data/ratgtex_silver_benchmark_balanced_len16384.tsv"

REBUILD_SCRIPT="${PROJECT_ROOT}/ratgetx/rebuild_sequences_for_alphagenome.py"

echo "========================================================================"
echo "Rebuild sequences for AlphaGenome (len=16384)"
echo "PROJECT_ROOT: ${PROJECT_ROOT}"
echo "------------------------------------------------------------------------"

# Ensure output dir exists
mkdir -p "$(dirname "${OUTPUT_TSV}")"

# FASTA presence check (skip index build)
echo "[Step 1/2] Using existing FASTA (and .fai if present); skipping index build."
if [ ! -f "${FASTA_FILE}" ]; then
  echo "ERROR: FASTA not found at: ${FASTA_FILE}" >&2
  exit 1
fi
echo "------------------------------------------------------------------------"

# Run rebuild script
echo "[Step 2/2] Rebuilding sequences from variant coordinates (16384 bp windows)…"
"${PYTHON_EXECUTABLE}" "${REBUILD_SCRIPT}" \
  --input_tsv "${INPUT_TSV}" \
  --fasta_file "${FASTA_FILE}" \
  --seq_len 16384 \
  --output_tsv "${OUTPUT_TSV}"

echo "------------------------------------------------------------------------"
echo "✅ Done."
echo "   Input : ${INPUT_TSV}"
echo "   Output: ${OUTPUT_TSV}"
echo "========================================================================"