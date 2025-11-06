#!/bin/bash
# Run SpliceTransformer on ChickenGTEx VCF (server script)

set -euo pipefail

# --- Server Config ---
ENV_PATH="/mnt/userdata4/splicing/conda_envs/sptransformer-env"
SPTR_ROOT="/mnt/userdata4/splicing/SpliceTransformer/SpliceTransformer-master"
INPUT_VCF="/mnt/userdata4/splicing/SpliceAI/chicken_data/spliceai_input_chicken.vcf"
OUTPUT_VCF="/mnt/userdata4/splicing/SpliceTransformer/data/SpliceTransformer_chickengtex_predictions.vcf"
REF_BUILD="hg38"

mkdir -p "$(dirname "$OUTPUT_VCF")"

echo "==============================================="
echo "SpliceTransformer ChickenGTEx Prediction"
echo "==============================================="
echo "Input:  $INPUT_VCF"
echo "Output: $OUTPUT_VCF"
echo "Reference: $REF_BUILD (Note: Using hg38 model for chicken - cross-species test)"
echo ""

# Check prerequisites
if [ ! -f "$INPUT_VCF" ]; then
    echo "ERROR: Input VCF not found: $INPUT_VCF"
    echo "Please run SpliceAI chicken pipeline first to generate the VCF"
    exit 1
fi

if [ ! -f "$SPTR_ROOT/sptransformer.py" ]; then
    echo "ERROR: SpliceTransformer script not found: $SPTR_ROOT/sptransformer.py"
    exit 1
fi

# Activate environment
PYTHON_EXE="$ENV_PATH/bin/python"
if command -v conda >/dev/null 2>&1; then
  source "$(conda info --base)/etc/profile.d/conda.sh" || true
  conda activate "$ENV_PATH" || true
  PYTHON_EXE="$(command -v python)"
fi

echo "Python: $PYTHON_EXE"
echo ""

cd "$SPTR_ROOT"

# Run SpliceTransformer
echo "Running SpliceTransformer..."
"$PYTHON_EXE" sptransformer.py \
  --input "$INPUT_VCF" \
  --output "$OUTPUT_VCF" \
  --reference "$REF_BUILD"

# Check output
if [ ! -s "$OUTPUT_VCF" ]; then
    echo "ERROR: SpliceTransformer prediction failed: $OUTPUT_VCF is empty or missing"
    exit 1
fi

echo ""
echo "==============================================="
echo "SpliceTransformer Prediction Complete"
echo "==============================================="
echo "Output: $OUTPUT_VCF"
echo "Variants predicted: $(grep -v "^#" "$OUTPUT_VCF" | wc -l)"
echo ""
echo "Next steps:"
echo "1. Copy results back to local machine:"
echo "   scp yinuos@mlerp-monash-node04:$OUTPUT_VCF \\"
echo "       /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/chickenGTEx/baselines/spliceTransformer/results/"
echo ""
echo "2. Run evaluation:"
echo "   cd /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/chickenGTEx/baselines/spliceTransformer"
echo "   python evaluate_splicetransformer_results.py"

