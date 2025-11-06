#!/bin/bash
# SpliceAI preprocessing pipeline (server):
# 1) Convert TSV -> VCF
# 2) Run SpliceAI to annotate VCF (uses server reference genome)
# 3) Parse SpliceAI-annotated VCF -> TSV with scores

set -euo pipefail

# --- Configuration (absolute paths for robustness) ---
CONDA_ENV_PATH="/mnt/userdata4/splicing/conda_envs/spliceai_env"
ROOT="/mnt/userdata4/splicing"
PYTHON_EXE="$CONDA_ENV_PATH/bin/python"
SPLICEAI_BIN="$CONDA_ENV_PATH/bin/spliceai"

INPUT_TSV="$ROOT/Evo2Splicevardb/splicevardb_data/evo2_splicevardb_dataset_dedup.tsv"
OUT_DIR="$ROOT/SpliceAI/data"
OUT_VCF="$OUT_DIR/spliceai_input_splicevardb.vcf"
PRED_VCF="$OUT_DIR/spliceai_predictions_splicevardb.vcf"
OUT_TSV="$OUT_DIR/spliceai_parsed_scores.tsv"
REF_FASTA="$ROOT/Evo2Splicevardb/splicevardb_data/hg38.fa"
ANNOT_ALIAS="grch38"

mkdir -p "$OUT_DIR"

# 1) TSV -> VCF
echo "[1/3] Converting TSV -> VCF: $OUT_VCF"
"${PYTHON_EXE:-python}" "$ROOT/SpliceAI/tsv_to_vcf.py" \
    --input_tsv "$INPUT_TSV" \
    --output_vcf "$OUT_VCF"

# 2) Run SpliceAI
echo "[2/3] Running SpliceAI -> $PRED_VCF"
"${SPLICEAI_BIN:-spliceai}" \
    -I "$OUT_VCF" \
    -O "$PRED_VCF" \
    -R "$REF_FASTA" \
    -A "$ANNOT_ALIAS" \
    -D 4999

# Verify predictions VCF
if [ ! -s "$PRED_VCF" ]; then
    echo "ERROR: Predictions VCF not created: $PRED_VCF" >&2
    echo "Hint: check FASTA path and chromosome naming; ensure SpliceAI finished successfully." >&2
    exit 1
fi


# 3) Parse SpliceAI-annotated VCF -> TSV
echo "[3/3] Parsing SpliceAI VCF -> TSV: $OUT_TSV"
"${PYTHON_EXE:-python}" "$ROOT/SpliceAI/parse_spliceai_output.py" \
    -I "$PRED_VCF" \
    -O "$OUT_TSV"

# Verify parsed TSV
if [ ! -s "$OUT_TSV" ]; then
    echo "ERROR: Parsed TSV not created: $OUT_TSV" >&2
    exit 1
fi



