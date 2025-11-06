#!/bin/bash
set -euo pipefail

# Activate MMSplice environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /mnt/userdata4/splicing/conda_envs/mmsplice-env

# Rat reference files (all on server)
GTF=/mnt/userdata4/splicing/SpliceAI/reference_genome/Rattus_norvegicus.Rnor_6.0.84.final.gtf
FASTA=/mnt/userdata4/splicing/SpliceAI/reference_genome/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa
VCF=/mnt/userdata4/splicing/SpliceAI/rat_data/spliceai_input_rat.vcf
OUTDIR=/mnt/userdata4/splicing/MMSplice/ratgtex
OUTPUT_CSV=$OUTDIR/mmsplice_ratgtex_scores.csv

mkdir -p "$OUTDIR"

echo "=========================================="
echo "MMSplice Rat (ratGTEx) Prediction Pipeline"
echo "=========================================="
echo "GTF:    $GTF"
echo "FASTA:  $FASTA"
echo "VCF:    $VCF"
echo "Output: $OUTPUT_CSV"
echo ""

# Check file existence
echo "Checking input files..."
for file in "$GTF" "$FASTA" "$VCF"; do
    if [ ! -f "$file" ]; then
        echo "ERROR: File not found: $file"
        exit 1
    fi
done
echo "✓ All input files exist"
echo ""

# Check FASTA index
if [ ! -f "${FASTA}.fai" ]; then
    echo "WARNING: FASTA index not found, creating..."
    samtools faidx "$FASTA"
    echo "✓ FASTA index created"
else
    echo "✓ FASTA index exists"
fi
echo ""

# Check chromosome naming consistency (optional diagnostics)
echo "Verifying chromosome format (rat uses 1,2,3 not chr1,chr2,chr3)..."
vcf_chr=$(grep -v "^#" "$VCF" 2>/dev/null | head -1 | cut -f1 || echo "unknown")
echo "  VCF first chromosome: $vcf_chr"
echo "  Expected format: numerical (1, 2, 3..., X) without chr prefix"
echo ""

# NOTE: Rat genome (rn6/Rnor_6.0) uses numerical chromosome naming (1, 2, 3)
# without 'chr' prefix (Ensembl standard), unlike human hg38 which uses chr1, chr2.
# The FASTA and GTF should also use numerical format for consistency.
# If chromosome naming mismatch occurs, the Python script will report the error.

# Run MMSplice following README example
echo "Running MMSplice predictions (this may take several hours)..."
echo "Start time: $(date)"
python - <<'PYTHON_CODE' "$GTF" "$FASTA" "$VCF" "$OUTPUT_CSV"
import sys
from mmsplice.vcf_dataloader import SplicingVCFDataloader
from mmsplice import MMSplice, predict_all_table

# Get arguments
gtf, fasta, vcf, output_csv = sys.argv[1:]

print(f"Loading data from VCF: {vcf}")
try:
    # Dataloader (tissue_specific=False, as ratGTEx is tissue-aggregated)
    dl = SplicingVCFDataloader(gtf, fasta, vcf, tissue_specific=False)
    print(f"✓ Data loaded successfully")
except ValueError as e:
    if "Fasta chrom names do not match" in str(e):
        print(f"ERROR: Chromosome naming mismatch detected!")
        print(f"  {e}")
        print(f"\nTroubleshooting:")
        print(f"  1. Check if FASTA uses 'chr1' while VCF uses '1' (or vice versa)")
        print(f"  2. Uncomment chromosome conversion in the shell script")
        print(f"  3. Or create VCF with matching chromosome names")
        sys.exit(1)
    else:
        raise

print("Initializing MMSplice model...")
model = MMSplice()

print("Running predictions...")
# Predict and return as dataframe (with pathogenicity and splicing_efficiency)
predictions = predict_all_table(model, dl, pathogenicity=True, splicing_efficiency=True)

# Save to CSV
predictions.to_csv(output_csv, index=False)
print(f"\n{'='*60}")
print(f"Predictions saved to: {output_csv}")
print(f"Total predictions: {len(predictions):,} rows")
print(f"Unique variants: {predictions['ID'].nunique():,}")
print(f"{'='*60}")
PYTHON_CODE

echo ""
echo "End time: $(date)"
echo "=========================================="
echo "MMSplice rat prediction completed!"
echo "Output: $OUTPUT_CSV"
echo "=========================================="

# Summary statistics
echo ""
echo "Quick summary:"
python - <<'PY' "$OUTPUT_CSV"
import sys
import pandas as pd

csv_file = sys.argv[1]
df = pd.read_csv(csv_file)

print(f"  Total predictions:     {len(df):,}")
print(f"  Unique variants (ID):  {df['ID'].nunique():,}")
print(f"  Columns: {', '.join(df.columns[:5])}...")
print(f"\nKey score distributions:")
for col in ['delta_logit_psi', 'pathogenicity', 'efficiency']:
    if col in df.columns:
        print(f"  {col}: mean={df[col].mean():.4f}, std={df[col].std():.4f}")
PY

echo ""
echo "Next steps:"
echo "  1. Download results: $OUTPUT_CSV"
echo "  2. Run evaluation script locally against ratGTEx labels"
echo "  3. Compare with SpliceAI, SpliceTransformer rat performance"
