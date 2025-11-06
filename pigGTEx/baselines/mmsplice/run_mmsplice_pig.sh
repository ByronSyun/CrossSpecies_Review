#!/bin/bash
set -euo pipefail

# Activate MMSplice environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /mnt/userdata4/splicing/conda_envs/mmsplice-env

# Pig reference files (all on server)
GTF=/mnt/userdata4/splicing/PigGTEx/reference_genome/Sus_scrofa.Sscrofa11.1.gtf
FASTA=/mnt/userdata4/splicing/PigGTEx/reference_genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa
VCF=/mnt/userdata4/splicing/MMSplice/piggtex/piggtex_input.vcf
OUTDIR=/mnt/userdata4/splicing/MMSplice/piggtex
OUTPUT_CSV=$OUTDIR/mmsplice_piggtex_scores.csv

# Input TSV for VCF creation
INPUT_TSV=/mnt/userdata4/splicing/PigGTEx/processed_data/piggtex_silver_benchmark_balanced.tsv

mkdir -p "$OUTDIR"

echo "=========================================="
echo "MMSplice Pig (PigGTEx) Prediction Pipeline"
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

print(f"Loaded {len(df)} pig variants from: {input_tsv}")

# VCF header for Sscrofa11.1 (Ensembl format)
vcf_header = [
    "##fileformat=VCFv4.2",
    "##source=piggtex_mmsplice_converter",
    "##reference=Sscrofa11.1",
    "##contig=<ID=1,length=274330532>",
    "##contig=<ID=2,length=151935994>",
    "##contig=<ID=3,length=132848913>",
    "##contig=<ID=4,length=130910915>",
    "##contig=<ID=5,length=104526007>",
    "##contig=<ID=6,length=170843587>",
    "##contig=<ID=7,length=121844099>",
    "##contig=<ID=8,length=138966237>",
    "##contig=<ID=9,length=139009144>",
    "##contig=<ID=10,length=69359453>",
    "##contig=<ID=11,length=79169978>",
    "##contig=<ID=12,length=61602749>",
    "##contig=<ID=13,length=208334590>",
    "##contig=<ID=14,length=141755446>",
    "##contig=<ID=15,length=140412725>",
    "##contig=<ID=16,length=79944280>",
    "##contig=<ID=17,length=63494081>",
    "##contig=<ID=18,length=55982971>",
    "##contig=<ID=X,length=125939595>",
    "##contig=<ID=Y,length=43547828>",
    "##contig=<ID=MT,length=16613>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
]

with open(output_vcf, 'w') as f:
    for line in vcf_header:
        f.write(line + '\n')
    
    for _, row in df.iterrows():
        vid = row['variant_id']
        # Format: '13:65362642_C/T' -> chrom=13, pos=65362642, ref=C, alt=T
        chrom_pos, ref_alt = vid.split('_')
        chrom, pos = chrom_pos.split(':')
        ref, alt = ref_alt.split('/')
        
        # VCF line (simple Ensembl chromosome naming)
        f.write(f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t.\tPASS\t.\n")

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
echo "  Expected format: numerical (1-18, X, Y, MT) without chr prefix"
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
    # Dataloader (tissue_specific=False, as PigGTEx is tissue-aggregated)
    dl = SplicingVCFDataloader(gtf, fasta, vcf, tissue_specific=False)
    print(f"✓ Data loaded successfully")
except ValueError as e:
    if "Fasta chrom names do not match" in str(e):
        print(f"ERROR: Chromosome naming mismatch detected!")
        print(f"  {e}")
        print(f"\nTroubleshooting:")
        print(f"  1. Check if FASTA uses 'chr1' while VCF uses '1' (or vice versa)")
        print(f"  2. Ensure all files use Sscrofa11.1 chromosome naming")
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
echo "MMSplice pig prediction completed!"
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
echo "  2. Run evaluation script locally against PigGTEx labels"
echo "  3. Compare with other pig baseline performance"

