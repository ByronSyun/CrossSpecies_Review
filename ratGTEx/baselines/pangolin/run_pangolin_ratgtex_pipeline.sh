#!/bin/bash

# Comprehensive Pangolin RatGTEx Pipeline
# This script runs the complete pipeline from data conversion to evaluation

set -euo pipefail  # Exit on any error

# === Configuration ===
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(dirname "$(dirname "$(dirname "$(dirname "$SCRIPT_DIR")")")")"

# Data paths
RATGTEX_DATA="$BASE_DIR/data/processed_data/ratGTEx/ratgtex_silver_benchmark_balanced.tsv"

# Output paths
VCF_OUTPUT="$SCRIPT_DIR/ratgtex_for_pangolin.vcf"
MAPPING_OUTPUT="$SCRIPT_DIR/ratgtex_for_pangolin_mapping.csv"
PREDICTIONS_OUTPUT="$SCRIPT_DIR/pangolin_ratgtex_results.json"
EVALUATION_DIR="$SCRIPT_DIR"

# Server paths (will be used when running on server)
SERVER_DATA="/mnt/userdata4/splicing/ratgtex/processed_data/ratgtex_silver_benchmark_balanced.tsv"

# === Setup Logging ===
LOG_FILE="$SCRIPT_DIR/pangolin_pipeline_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "=== Pangolin RatGTEx Pipeline Started ==="
echo "Timestamp: $(date)"
echo "Script directory: $SCRIPT_DIR"
echo "Base directory: $BASE_DIR"
echo "Log file: $LOG_FILE"
echo "============================================="

# === Environment Detection ===
if [[ -d "/mnt/userdata4/splicing/" ]]; then
    echo "🖥️  Server environment detected"
    IS_SERVER=true
    INPUT_DATA="$SERVER_DATA"
else
    echo "💻 Local environment detected"
    IS_SERVER=false
    INPUT_DATA="$RATGTEX_DATA"
fi

# === Step 1: Validate Input Data ===
echo ""
echo "Step 1: Validating input data..."
if [[ ! -f "$INPUT_DATA" ]]; then
    echo "❌ Error: Input data file not found: $INPUT_DATA"
    if [[ "$IS_SERVER" == "false" ]]; then
        echo "For local testing, ensure the RatGTEx data has been preprocessed."
        echo "Expected location: $RATGTEX_DATA"
    fi
    exit 1
fi

echo "✅ Input data found: $INPUT_DATA"
echo "   File size: $(du -h "$INPUT_DATA" | cut -f1)"
echo "   Line count: $(wc -l < "$INPUT_DATA")"

# === Step 2: Convert TSV to VCF Format ===
echo ""
echo "Step 2: Converting RatGTEx TSV to VCF format..."
python3 "$SCRIPT_DIR/convert_ratgtex_to_vcf.py" \
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

# === Step 3: Environment-Specific Dependency Check ===
echo ""
echo "Step 3: Checking Pangolin dependencies..."

if [[ "$IS_SERVER" == "true" ]]; then
    echo "Checking server environment..."
    
    # Check conda environment
    PANGOLIN_ENV="/mnt/userdata4/splicing/conda_envs/pangolin-env"
    if [[ ! -d "$PANGOLIN_ENV" ]]; then
        echo "❌ Pangolin conda environment not found: $PANGOLIN_ENV"
        echo "Please create the environment first."
        exit 1
    fi
    
    # Check rat reference genome
    RAT_REFERENCE="/mnt/userdata4/splicing/references/rat/ratNor6.fa"
    if [[ ! -f "$RAT_REFERENCE" ]]; then
        echo "❌ Rat reference genome not found: $RAT_REFERENCE"
        echo "Please download ratNor6.fa to this location."
        exit 1
    fi
    
    # Check rat annotation database
    RAT_ANNOTATION="/mnt/userdata4/splicing/references/rat/ratNor6.annotation.db"
    if [[ ! -f "$RAT_ANNOTATION" ]]; then
        echo "❌ Rat annotation database not found: $RAT_ANNOTATION"
        echo "Please create the annotation database using Pangolin's create_db.py"
        exit 1
    fi
    
    echo "✅ Server dependencies validated"
    
else
    echo "Local environment - dependency check..."
    
    # Check if pangolin is installed
    if ! command -v pangolin &> /dev/null; then
        echo "❌ Pangolin not found in PATH"
        echo "Please install Pangolin: pip install pangolin"
        exit 1
    fi
    
    # Check for reference files (provide guidance)
    if [[ ! -f "ratNor6.fa" ]] || [[ ! -f "ratNor6.annotation.db" ]]; then
        echo "⚠️  Warning: Rat reference files not found in current directory"
        echo "   Required files:"
        echo "   - ratNor6.fa (rat reference genome)"
        echo "   - ratNor6.annotation.db (annotation database)"
        echo ""
        echo "To set up locally:"
        echo "1. Download ratNor6.fa from UCSC or Ensembl"
        echo "2. Create annotation database:"
        echo "   python -c \"from pangolin.create_db import create_db; create_db('ratNor6.fa', 'ratNor6.gtf', 'ratNor6.annotation.db')\""
        echo ""
        echo "For now, showing what the command would be..."
    fi
    
    echo "✅ Local environment check completed"
fi

# === Step 4: Run Pangolin Prediction ===
echo ""
echo "Step 4: Running Pangolin predictions..."

# Note: The Python script handles environment detection internally
python3 "$SCRIPT_DIR/predict_ratgtex_pangolin.py" \
    --input_vcf "$VCF_OUTPUT" \
    --mapping_csv "$MAPPING_OUTPUT" \
    --output_json "$PREDICTIONS_OUTPUT"

PANGOLIN_EXIT_CODE=$?
if [[ $PANGOLIN_EXIT_CODE -ne 0 ]]; then
    echo "❌ Pangolin prediction failed with exit code $PANGOLIN_EXIT_CODE"
    
    if [[ "$IS_SERVER" == "false" ]]; then
        echo ""
        echo "For local development, this step requires:"
        echo "1. Rat reference genome (ratNor6.fa)"
        echo "2. Rat annotation database (ratNor6.annotation.db)"
        echo "3. Pangolin installation"
        echo ""
        echo "Consider running this pipeline on the server where dependencies are set up."
    fi
    
    exit 1
fi

echo "✅ Pangolin predictions completed"
echo "   Results saved to: $PREDICTIONS_OUTPUT"

# === Step 5: Evaluate Results ===
echo ""
echo "Step 5: Evaluating Pangolin results..."

python3 "$SCRIPT_DIR/evaluate_pangolin_results.py" \
    --input_json "$PREDICTIONS_OUTPUT" \
    --output_dir "$EVALUATION_DIR"

if [[ $? -ne 0 ]]; then
    echo "❌ Evaluation failed!"
    exit 1
fi

echo "✅ Evaluation completed"

# === Step 6: Summary ===
echo ""
echo "🎉 Pangolin RatGTEx Pipeline Completed Successfully!"
echo ""
echo "Generated files:"
echo "  📄 VCF input: $VCF_OUTPUT"
echo "  📄 Mapping file: $MAPPING_OUTPUT"
echo "  📄 Predictions: $PREDICTIONS_OUTPUT"
echo "  📄 Performance metrics: $EVALUATION_DIR/pangolin_performance_metrics.json"
echo "  📊 Performance curves: $EVALUATION_DIR/pangolin_ratgtex_performance_curves.png"
echo "  📊 Score distribution: $EVALUATION_DIR/pangolin_score_distribution.png"
echo "  📋 Summary report: $EVALUATION_DIR/pangolin_evaluation_summary.txt"
echo ""
echo "📋 Log file: $LOG_FILE"
echo ""

# Display key results if evaluation was successful
if [[ -f "$EVALUATION_DIR/pangolin_performance_metrics.json" ]]; then
    echo "🎯 Key Results:"
    python3 -c "
import json
try:
    with open('$EVALUATION_DIR/pangolin_performance_metrics.json', 'r') as f:
        metrics = json.load(f)
    print(f\"   AUROC: {metrics['auroc']:.4f}\")
    print(f\"   AUPRC: {metrics['auprc']:.4f}\")
    print(f\"   Accuracy: {metrics['accuracy']:.4f}\")
    print(f\"   F1-score: {metrics['f1_score']:.4f}\")
    print(f\"   Total variants: {metrics['total_variants']:,}\")
except:
    print('   Could not load metrics')
"
fi

echo ""
echo "Next steps:"
echo "  - Compare with other baseline models (AlphaGenome, SpliceTransformer)"
echo "  - Add results to cross-species comparison plots"
echo "  - Include in manuscript results section"
echo ""
echo "============================================="
echo "Pipeline completed at: $(date)"
