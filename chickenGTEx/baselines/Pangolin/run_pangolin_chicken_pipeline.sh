#!/bin/bash
# Pangolin ChickenGTEx Pipeline (Server Only)
# Converts chicken TSV to VCF and runs Pangolin prediction

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(dirname "$(dirname "$(dirname "$(dirname "$SCRIPT_DIR")")")")"

# Data paths
CHICKEN_DATA="$BASE_DIR/data/processed_data/chickenGTEx/chickengtex_silver_benchmark_balanced.tsv"
SERVER_DATA="/mnt/userdata4/splicing/ChickenGTEx/processed_data/chickengtex_silver_benchmark_balanced.tsv"

# Output paths
VCF_OUTPUT="$SCRIPT_DIR/chickengtex_for_pangolin.vcf"
MAPPING_OUTPUT="$SCRIPT_DIR/chickengtex_for_pangolin_mapping.csv"
PREDICTIONS_OUTPUT="$SCRIPT_DIR/results/pangolin_chickengtex_results.json"

LOG_FILE="$SCRIPT_DIR/pangolin_pipeline_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "=== Pangolin ChickenGTEx Pipeline Started ==="
echo "Timestamp: $(date)"
echo "============================================="

# Environment Detection
if [[ -d "/mnt/userdata4/splicing/" ]]; then
    echo "🖥️  Server environment detected"
    IS_SERVER=true
    INPUT_DATA="$SERVER_DATA"
else
    echo "💻 Local environment detected"
    IS_SERVER=false
    INPUT_DATA="$CHICKEN_DATA"
fi

# Step 1: Validate Input Data
echo ""
echo "Step 1: Validating input data..."
if [[ ! -f "$INPUT_DATA" ]]; then
    echo "❌ Error: Input data file not found: $INPUT_DATA"
    exit 1
fi

echo "✅ Input data found: $INPUT_DATA"
echo "   File size: $(du -h "$INPUT_DATA" | cut -f1)"
echo "   Line count: $(wc -l < "$INPUT_DATA")"

# Step 2: Check Dependencies
echo ""
echo "Step 2: Checking Pangolin dependencies..."

if [[ "$IS_SERVER" == "true" ]]; then
    echo "Checking server environment..."
    
    PANGOLIN_ENV="/mnt/userdata4/splicing/conda_envs/pangolin-env"
    if [[ ! -d "$PANGOLIN_ENV" ]]; then
        echo "❌ Pangolin conda environment not found: $PANGOLIN_ENV"
        echo "Please create the environment first."
        exit 1
    fi
    
    CHICKEN_REFERENCE="/mnt/userdata4/splicing/ChickenGTEx/reference_genome/Gallus_gallus.GRCg6a.dna.toplevel.fa"
    if [[ ! -f "$CHICKEN_REFERENCE" ]]; then
        echo "❌ Chicken reference genome not found: $CHICKEN_REFERENCE"
        echo "Please download Gallus_gallus.GRCg6a.dna.toplevel.fa to this location."
        exit 1
    fi
    
    CHICKEN_ANNOTATION="/mnt/userdata4/splicing/ChickenGTEx/reference_genome/GRCg6a.annotation.db"
    if [[ ! -f "$CHICKEN_ANNOTATION" ]]; then
        echo "❌ Chicken annotation database not found: $CHICKEN_ANNOTATION"
        echo "Please create the annotation database using Pangolin's create_db.py"
        exit 1
    fi
    
    echo "✅ Server dependencies validated"
    PYTHON_CMD="conda run -p $PANGOLIN_ENV python3"
else
    echo "Local environment - dependency check..."
    echo "✅ Local environment check completed"
    PYTHON_CMD="python3"
fi

# Step 3: Convert TSV to VCF Format
echo ""
echo "Step 3: Converting ChickenGTEx TSV to VCF format..."

# Check if VCF and mapping already exist
if [[ -f "$VCF_OUTPUT" && -s "$VCF_OUTPUT" && -f "$MAPPING_OUTPUT" && -s "$MAPPING_OUTPUT" ]]; then
    VCF_LINES=$(wc -l < "$VCF_OUTPUT")
    echo "✅ VCF files already exist (skipping conversion)"
    echo "   VCF output: $VCF_OUTPUT ($VCF_LINES lines)"
    echo "   Mapping file: $MAPPING_OUTPUT"
else
    echo "   Running VCF conversion..."
    ${PYTHON_CMD} "$SCRIPT_DIR/convert_chickengtex_to_vcf.py" \
        --input "$INPUT_DATA" \
        --output "$VCF_OUTPUT" \
        --sequence_length 8192

    if [[ $? -ne 0 ]]; then
        echo "❌ VCF conversion failed!"
        exit 1
    fi

    echo "✅ VCF conversion completed"
    echo "   VCF output: $VCF_OUTPUT"
    echo "   Mapping file: $MAPPING_OUTPUT"
fi

# Step 4: Run Pangolin Prediction
echo ""
echo "Step 4: Running Pangolin predictions..."

# Check if predictions already exist
if [[ -f "$PREDICTIONS_OUTPUT" && -s "$PREDICTIONS_OUTPUT" ]]; then
    PRED_SIZE=$(du -h "$PREDICTIONS_OUTPUT" | cut -f1)
    echo "✅ Pangolin predictions already exist (skipping prediction)"
    echo "   Results file: $PREDICTIONS_OUTPUT ($PRED_SIZE)"
    echo "   To re-run, delete the file and restart the pipeline."
else
    echo "   Running Pangolin prediction (this may take up to 6 hours for ~25k variants)..."
    
    ${PYTHON_CMD} "$SCRIPT_DIR/predict_chickengtex_pangolin.py" \
        --input_vcf "$VCF_OUTPUT" \
        --mapping_csv "$MAPPING_OUTPUT" \
        --output_json "$PREDICTIONS_OUTPUT"

    PANGOLIN_EXIT_CODE=$?
    if [[ $PANGOLIN_EXIT_CODE -ne 0 ]]; then
        echo "❌ Pangolin prediction failed with exit code $PANGOLIN_EXIT_CODE"
        exit 1
    fi

    echo "✅ Pangolin predictions completed"
    echo "   Results saved to: $PREDICTIONS_OUTPUT"
fi

echo ""
echo "🎉 Pangolin ChickenGTEx Pipeline Completed Successfully!"
echo "Generated files:"
echo "  📄 VCF input: $VCF_OUTPUT"
echo "  📄 Mapping file: $MAPPING_OUTPUT"
echo "  📄 Predictions: $PREDICTIONS_OUTPUT"
echo ""
echo "📋 Log file: $LOG_FILE"
echo ""
echo "============================================="
echo "Pipeline completed at: $(date)"

