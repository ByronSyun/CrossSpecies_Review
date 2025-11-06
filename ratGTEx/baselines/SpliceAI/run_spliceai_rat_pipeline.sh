#!/bin/bash
# SpliceAI pipeline for rat data (server):
# 1) Convert rat TSV -> VCF
# 2) Create rat annotation file from GTF (if needed)
# 3) Run SpliceAI to annotate VCF (uses rat reference genome)
# 4) Parse SpliceAI-annotated VCF -> TSV with scores

set -euo pipefail

# --- Configuration (server paths) ---
CONDA_ENV_PATH="/mnt/userdata4/splicing/conda_envs/spliceai_env"
ROOT="/mnt/userdata4/splicing"
PYTHON_EXE="$CONDA_ENV_PATH/bin/python"
SPLICEAI_BIN="$CONDA_ENV_PATH/bin/spliceai"

# Input/Output paths
INPUT_TSV="$ROOT/ratgetx/processed_data/ratgtex_silver_benchmark_balanced.tsv"
OUT_DIR="$ROOT/SpliceAI/rat_data"
OUT_VCF="$OUT_DIR/spliceai_input_rat.vcf"
PRED_VCF="$OUT_DIR/spliceai_predictions_rat.vcf"
OUT_TSV="$OUT_DIR/spliceai_parsed_scores_rat.tsv"

# Reference files (rat)
REF_FASTA="$ROOT/SpliceAI/reference_genome/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa"
RAT_GTF="$ROOT/SpliceAI/reference_genome/Rattus_norvegicus.Rnor_6.0.84.final.gtf"
RAT_ANNOTATION="$OUT_DIR/rat_rn6_annotation.txt"

# Script paths (upload these to server)
SCRIPT_DIR="$ROOT/SpliceAI"
TSV_TO_VCF_SCRIPT="$SCRIPT_DIR/tsv_to_vcf_rat.py"
CREATE_ANNOTATION_SCRIPT="$SCRIPT_DIR/create_rat_annotation.py"
PARSE_OUTPUT_SCRIPT="$SCRIPT_DIR/parse_spliceai_output.py"

mkdir -p "$OUT_DIR"
mkdir -p "$SCRIPT_DIR"

echo "=== SpliceAI Rat Pipeline ==="
echo "Input TSV: $INPUT_TSV"
echo "Reference FASTA: $REF_FASTA"
echo "Reference GTF: $RAT_GTF"
echo "Output directory: $OUT_DIR"

# Check if input files exist
if [ ! -f "$INPUT_TSV" ]; then
    echo "ERROR: Input TSV not found: $INPUT_TSV" >&2
    exit 1
fi

if [ ! -f "$REF_FASTA" ]; then
    echo "ERROR: Reference FASTA not found: $REF_FASTA" >&2
    echo "Please ensure rat reference genome is available at the specified path." >&2
    exit 1
fi

# Step 1: TSV -> VCF conversion
echo "[1/4] Converting rat TSV -> VCF: $OUT_VCF"
"$PYTHON_EXE" "$TSV_TO_VCF_SCRIPT" \
    --input_tsv "$INPUT_TSV" \
    --output_vcf "$OUT_VCF"

# Verify VCF creation
if [ ! -s "$OUT_VCF" ]; then
    echo "ERROR: VCF not created: $OUT_VCF" >&2
    exit 1
fi

# Step 2: Create rat annotation file (if not exists)
if [ ! -f "$RAT_ANNOTATION" ]; then
    echo "[2/4] Creating rat annotation file from GTF: $RAT_ANNOTATION"
    if [ ! -f "$RAT_GTF" ]; then
        echo "ERROR: Rat GTF not found: $RAT_GTF" >&2
        echo "Please download rat GTF from Ensembl or provide the correct path." >&2
        exit 1
    fi
    
    "$PYTHON_EXE" "$CREATE_ANNOTATION_SCRIPT" \
        --input_gtf "$RAT_GTF" \
        --output_txt "$RAT_ANNOTATION"
else
    echo "[2/4] Using existing rat annotation file: $RAT_ANNOTATION"
fi

# Verify annotation file
if [ ! -s "$RAT_ANNOTATION" ]; then
    echo "ERROR: Rat annotation file not created: $RAT_ANNOTATION" >&2
    exit 1
fi

# Step 3: Run SpliceAI
echo "[3/4] Running SpliceAI on rat data -> $PRED_VCF"
echo "This may take a while for large datasets..."

"$SPLICEAI_BIN" \
    -I "$OUT_VCF" \
    -O "$PRED_VCF" \
    -R "$REF_FASTA" \
    -A "$RAT_ANNOTATION" \
    -D 4999

# Verify predictions VCF
if [ ! -s "$PRED_VCF" ]; then
    echo "ERROR: Predictions VCF not created: $PRED_VCF" >&2
    echo "Hints:" >&2
    echo "  - Check if rat FASTA is indexed (should have .fai file)" >&2
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
echo "3. Compare with AlphaGenome and Pangolin results"

# Show sample results
echo ""
echo "Sample results:"
head -5 "$OUT_TSV"
