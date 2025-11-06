#!/bin/bash
# Run SpliceTransformer on SpliceVarDB VCF (server-friendly, concise)

set -euo pipefail

# --- Config ---
ENV_PATH="/mnt/userdata4/splicing/conda_envs/sptransformer-env"
SPTR_ROOT="/mnt/userdata4/splicing/SpliceTransformer/SpliceTransformer-master"
INPUT_VCF="/mnt/userdata4/splicing/SpliceAI/data/evo2_splicevardb_dataset_dedup.vcf"
OUTPUT_VCF="/mnt/userdata4/splicing/SpliceTransformer/data/SpliceTransformer_splicevardb_predictions.vcf"
REF_BUILD="hg38"

mkdir -p "$(dirname "$OUTPUT_VCF")"

echo "--- SpliceTransformer run ---"
echo "Env:     $ENV_PATH"
echo "Script:  $SPTR_ROOT/sptransformer.py"
echo "Input:   $INPUT_VCF"
echo "Output:  $OUTPUT_VCF"
echo "Reference: $REF_BUILD"

# Activate environment (prefer conda; fallback to direct python)
PYTHON_EXE="$ENV_PATH/bin/python"
if command -v conda >/dev/null 2>&1; then
  # shellcheck source=/dev/null
  source "$(conda info --base)/etc/profile.d/conda.sh" || true
  conda activate "$ENV_PATH" || true
  PYTHON_EXE="$(command -v python)"
fi

echo "Using Python: $PYTHON_EXE"

# Ensure working directory is the repository root so relative resource paths resolve
cd "$SPTR_ROOT"

"$PYTHON_EXE" sptransformer.py \
  --input "$INPUT_VCF" \
  --output "$OUTPUT_VCF" \
  --reference "$REF_BUILD"

echo "--- Done ---"


