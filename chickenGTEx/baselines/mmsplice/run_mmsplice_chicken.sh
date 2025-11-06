#!/bin/bash
set -euo pipefail

# Activate MMSplice environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /mnt/userdata4/splicing/conda_envs/mmsplice-env

# Chicken reference files (all on server)
GTF=/mnt/userdata4/splicing/ChickenGTEx/reference_genome/Gallus_gallus.GRCg6a.gtf
FASTA=/mnt/userdata4/splicing/ChickenGTEx/reference_genome/Gallus_gallus.GRCg6a.dna.toplevel.fa
VCF=/mnt/userdata4/splicing/MMSplice/chickengtex/chickengtex_input.vcf
OUTDIR=/mnt/userdata4/splicing/MMSplice/chickengtex
OUTPUT_CSV=$OUTDIR/mmsplice_chickengtex_scores.csv

# Input TSV for VCF creation
INPUT_TSV=/mnt/userdata4/splicing/ChickenGTEx/processed_data/chickengtex_silver_benchmark_balanced.tsv

mkdir -p "$OUTDIR"

echo "=========================================="
echo "MMSplice Chicken (ChickenGTEx) Prediction Pipeline"
echo "=========================================="
echo "GTF:    $GTF"
echo "FASTA:  $FASTA"
echo "VCF:    $VCF"
echo "Output: $OUTPUT_CSV"
echo ""

# Check file existence
echo "Checking reference files..."
for file in "$GTF" "$FASTA" "$INPUT_TSV"; do
    if [ ! -f "$file" ]; then
        echo "ERROR: File not found: $file"
        exit 1
    fi
done
echo "✓ All reference files exist"
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

# Create VCF if not exists
if [ ! -f "$VCF" ]; then
    echo "Creating VCF from TSV (Ensembl chromosome naming)..."
    python - <<'VCF_CREATE' "$INPUT_TSV" "$VCF"
import sys
import pandas as pd

input_tsv, output_vcf = sys.argv[1:]

# Read TSV (no header)
df = pd.read_csv(input_tsv, sep='\t', header=None, 
                 names=['variant_id', 'ref_seq', 'alt_seq', 'label', 'tissue_id'])

print(f"Loaded {len(df)} chicken variants from: {input_tsv}")

# VCF header for GRCg6a (Ensembl format)
vcf_header = [
    "##fileformat=VCFv4.2",
    "##source=chickengtex_mmsplice_converter",
    "##reference=GRCg6a",
    "##contig=<ID=1,length=195276750>",
    "##contig=<ID=2,length=149385596>",
    "##contig=<ID=3,length=111300296>",
    "##contig=<ID=4,length=90912899>",
    "##contig=<ID=5,length=60034740>",
    "##contig=<ID=6,length=34964064>",
    "##contig=<ID=7,length=36338314>",
    "##contig=<ID=8,length=29809473>",
    "##contig=<ID=9,length=23675360>",
    "##contig=<ID=10,length=20730876>",
    "##contig=<ID=11,length=20187453>",
    "##contig=<ID=12,length=19953944>",
    "##contig=<ID=13,length=18403236>",
    "##contig=<ID=14,length=15830388>",
    "##contig=<ID=15,length=12724206>",
    "##contig=<ID=16,length=9007986>",
    "##contig=<ID=17,length=11082006>",
    "##contig=<ID=18,length=10534285>",
    "##contig=<ID=19,length=10021710>",
    "##contig=<ID=20,length=13896712>",
    "##contig=<ID=21,length=6990948>",
    "##contig=<ID=22,length=5913949>",
    "##contig=<ID=23,length=6179846>",
    "##contig=<ID=24,length=6349559>",
    "##contig=<ID=25,length=2175668>",
    "##contig=<ID=26,length=5150254>",
    "##contig=<ID=27,length=5257854>",
    "##contig=<ID=28,length=5452033>",
    "##contig=<ID=29,length=5493311>",
    "##contig=<ID=30,length=5152642>",
    "##contig=<ID=31,length=3006013>",
    "##contig=<ID=32,length=3129577>",
    "##contig=<ID=33,length=3359810>",
    "##contig=<ID=Z,length=82528814>",
    "##contig=<ID=W,length=6159098>",
    "##contig=<ID=MT,length=16775>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
]

with open(output_vcf, 'w') as f:
    for line in vcf_header:
        f.write(line + '\n')
    
    for _, row in df.iterrows():
        vid = row['variant_id']
        # Format: 'chrom:pos_ref/alt' (after normalization)
        try:
            chrom_pos, ref_alt = vid.split('_', 1)  # Split only on first underscore
            chrom, pos = chrom_pos.split(':')
            ref, alt = ref_alt.split('/')
            
            # VCF line (Ensembl chromosome naming)
            f.write(f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t.\tPASS\t.\n")
        except Exception as e:
            print(f"Warning: Failed to parse variant_id '{vid}': {e}")
            continue

print(f"✓ VCF created: {output_vcf} ({len(df)} variants)")
VCF_CREATE
    
    if [ ! -f "$VCF" ]; then
        echo "ERROR: VCF creation failed"
        exit 1
    fi
else
    echo "✓ VCF already exists: $VCF"
fi
echo ""

# Check chromosome naming consistency
echo "Verifying chromosome format (Ensembl naming)..."
vcf_chr=$(grep -v "^#" "$VCF" 2>/dev/null | head -1 | cut -f1 || echo "unknown")
echo "  VCF first chromosome: $vcf_chr"
echo "  Expected format: numerical (1-33, Z, W, MT) without chr prefix"
echo ""

# Run MMSplice
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
    # Dataloader (tissue_specific=False, as ChickenGTEx is tissue-aggregated)
    dl = SplicingVCFDataloader(gtf, fasta, vcf, tissue_specific=False)
    print(f"✓ Data loaded successfully")
except ValueError as e:
    if "Fasta chrom names do not match" in str(e):
        print(f"ERROR: Chromosome naming mismatch detected!")
        print(f"  {e}")
        print(f"\nTroubleshooting:")
        print(f"  1. Check if FASTA uses 'chr1' while VCF uses '1' (or vice versa)")
        print(f"  2. Ensure all files use GRCg6a chromosome naming")
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
echo "MMSplice chicken prediction completed!"
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
echo "  2. Run evaluation script locally against ChickenGTEx labels"
echo "  3. Compare with other chicken baseline performance"

