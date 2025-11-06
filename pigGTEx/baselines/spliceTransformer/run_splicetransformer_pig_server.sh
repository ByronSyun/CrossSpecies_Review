#!/bin/bash
# Run SpliceTransformer on PigGTEx VCF (server script)

set -euo pipefail

# --- Server Config ---
ENV_PATH="/mnt/userdata4/splicing/conda_envs/sptransformer-env"
SPTR_ROOT="/mnt/userdata4/splicing/SpliceTransformer/SpliceTransformer-master"
INPUT_VCF="/mnt/userdata4/splicing/SpliceAI/pig_data/spliceai_input_pig.vcf"
OUTPUT_VCF="/mnt/userdata4/splicing/SpliceTransformer/data/SpliceTransformer_piggtex_predictions.vcf"
REF_BUILD="hg38"

mkdir -p "$(dirname "$OUTPUT_VCF")"

echo "Input:  $INPUT_VCF"
echo "Output: $OUTPUT_VCF"

# Check prerequisites
[ ! -f "$INPUT_VCF" ] && { echo "Input VCF not found: $INPUT_VCF"; exit 1; }
[ ! -f "$SPTR_ROOT/sptransformer.py" ] && { echo "SpliceTransformer script not found"; exit 1; }

# Activate environment
PYTHON_EXE="$ENV_PATH/bin/python"
if command -v conda >/dev/null 2>&1; then
  source "$(conda info --base)/etc/profile.d/conda.sh" || true
  conda activate "$ENV_PATH" || true
  PYTHON_EXE="$(command -v python)"
fi

cd "$SPTR_ROOT"

# Run SpliceTransformer
"$PYTHON_EXE" sptransformer.py \
  --input "$INPUT_VCF" \
  --output "$OUTPUT_VCF" \
  --reference "$REF_BUILD"

# Check output
if [ ! -s "$OUTPUT_VCF" ]; then
    echo "SpliceTransformer prediction failed: $OUTPUT_VCF"
    exit 1
fi

echo "Completed: $OUTPUT_VCF"
echo "Variants: $(grep -v "^#" "$OUTPUT_VCF" | wc -l)"


