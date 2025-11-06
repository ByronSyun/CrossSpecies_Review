#!/bin/bash
# SpliceAI pipeline for chicken data (server):
# 1) Convert chicken TSV -> VCF
# 2) Create chicken annotation file from GTF (if needed)
# 3) Run SpliceAI to annotate VCF (uses chicken reference genome)
# 4) Parse SpliceAI-annotated VCF -> TSV with scores

set -euo pipefail

# --- Configuration (server paths) ---
CONDA_ENV_PATH="/mnt/userdata4/splicing/conda_envs/spliceai_env"
ROOT="/mnt/userdata4/splicing"
PYTHON_EXE="$CONDA_ENV_PATH/bin/python"
SPLICEAI_BIN="$CONDA_ENV_PATH/bin/spliceai"

# Input/Output paths
INPUT_TSV="$ROOT/ChickenGTEx/processed_data/chickengtex_silver_benchmark_balanced.tsv"
OUT_DIR="$ROOT/SpliceAI/chicken_data"
OUT_VCF="$OUT_DIR/spliceai_input_chicken.vcf"
PRED_VCF="$OUT_DIR/spliceai_predictions_chicken.vcf"
OUT_TSV="$OUT_DIR/spliceai_parsed_scores_chicken.tsv"

# Reference files (chicken) - use Ensembl FASTA to match GTF
REF_FASTA="$ROOT/ChickenGTEx/reference_genome/Gallus_gallus.GRCg6a.dna.toplevel.fa"
CHICKEN_GTF="$ROOT/ChickenGTEx/reference_genome/Gallus_gallus.GRCg6a.102.gtf"
CHICKEN_ANNOTATION="$OUT_DIR/chicken_grcg6a_annotation.txt"

# Script paths (upload these to server)
SCRIPT_DIR="$ROOT/SpliceAI"
TSV_TO_VCF_SCRIPT="$SCRIPT_DIR/tsv_to_vcf_chicken.py"
CREATE_ANNOTATION_SCRIPT="$SCRIPT_DIR/create_chicken_annotation.py"
PARSE_OUTPUT_SCRIPT="$SCRIPT_DIR/parse_spliceai_output.py"

mkdir -p "$OUT_DIR"
mkdir -p "$SCRIPT_DIR"

echo "=== SpliceAI Chicken Pipeline ==="
echo "Input TSV: $INPUT_TSV"
echo "Reference FASTA: $REF_FASTA"
echo "Reference GTF: $CHICKEN_GTF"
echo "Output directory: $OUT_DIR"

# Check if input files exist
if [ ! -f "$INPUT_TSV" ]; then
    echo "ERROR: Input TSV not found: $INPUT_TSV" >&2
    exit 1
fi

if [ ! -f "$REF_FASTA" ]; then
    echo "ERROR: Reference FASTA not found: $REF_FASTA" >&2
    echo "Please ensure chicken reference genome is available at the specified path." >&2
    exit 1
fi

# Step 1: TSV -> VCF conversion
echo "[1/4] Converting chicken TSV -> VCF: $OUT_VCF"
"$PYTHON_EXE" "$TSV_TO_VCF_SCRIPT" \
    --input_tsv "$INPUT_TSV" \
    --output_vcf "$OUT_VCF"

# Verify VCF creation
if [ ! -s "$OUT_VCF" ]; then
    echo "ERROR: VCF not created: $OUT_VCF" >&2
    exit 1
fi

# Step 2: Create chicken annotation file (if not exists)
if [ ! -f "$CHICKEN_ANNOTATION" ]; then
    echo "[2/4] Creating chicken annotation file from GTF: $CHICKEN_ANNOTATION"
    if [ ! -f "$CHICKEN_GTF" ]; then
        echo "ERROR: Chicken GTF not found: $CHICKEN_GTF" >&2
        echo "Please download chicken GTF from Ensembl or provide the correct path." >&2
        exit 1
    fi
    
    "$PYTHON_EXE" "$CREATE_ANNOTATION_SCRIPT" \
        --input_gtf "$CHICKEN_GTF" \
        --output_txt "$CHICKEN_ANNOTATION"
else
    echo "[2/4] Using existing chicken annotation file: $CHICKEN_ANNOTATION"
fi

# Verify annotation file
if [ ! -s "$CHICKEN_ANNOTATION" ]; then
    echo "ERROR: Chicken annotation file not created: $CHICKEN_ANNOTATION" >&2
    exit 1
fi

# Step 3: Run SpliceAI
echo "[3/4] Running SpliceAI on chicken data -> $PRED_VCF"
echo "This may take a while for large datasets..."

"$SPLICEAI_BIN" \
    -I "$OUT_VCF" \
    -O "$PRED_VCF" \
    -R "$REF_FASTA" \
    -A "$CHICKEN_ANNOTATION" \
    -D 4999

# Verify predictions VCF
if [ ! -s "$PRED_VCF" ]; then
    echo "ERROR: Predictions VCF not created: $PRED_VCF" >&2
    echo "Hints:" >&2
    echo "  - Check if chicken FASTA is indexed (should have .fai file)" >&2
    echo "  - Ensure chromosome naming consistency between VCF, FASTA, and annotation" >&2
    echo "  - Check SpliceAI logs for specific error messages" >&2
    exit 1
fi

# Step 4: Parse SpliceAI-annotated VCF -> TSV
echo "[4/4] Parsing SpliceAI VCF -> TSV: $OUT_TSV"
"$PYTHON_EXE" "$PARSE_OUTPUT_SCRIPT" \
    -I "$PRED_VCF" \
    -O "$OUT_TSV"

# Verify parsed TSV
if [ ! -s "$OUT_TSV" ]; then
    echo "ERROR: Parsed TSV not created: $OUT_TSV" >&2
    exit 1
fi

echo "=== Pipeline completed successfully! ==="
echo "Results saved to: $OUT_TSV"
echo ""
echo "Next steps:"
echo "1. Download the results file to local machine"
echo "2. Run evaluation script to compute metrics"
echo "3. Compare with other baseline models"

# Show sample results
echo ""
echo "Sample results:"
head -5 "$OUT_TSV"

