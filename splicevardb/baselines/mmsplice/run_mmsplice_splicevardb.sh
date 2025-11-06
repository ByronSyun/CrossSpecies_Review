#!/bin/bash
set -euo pipefail

# Activate MMSplice environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /mnt/userdata4/splicing/conda_envs/mmsplice-env

# Paths (all required files ready on server)
GTF=/mnt/userdata4/splicing/SpliceAI/reference_genome/gencode.v38.annotation.gtf
FASTA=/mnt/userdata4/splicing/Evo2Splicevardb/splicevardb_data/hg38.fa
VCF_ORIG=/mnt/userdata4/splicing/SpliceAI/data/spliceai_input_splicevardb.vcf
OUTDIR=/mnt/userdata4/splicing/MMSplice/splicevardb
VCF_FIXED=$OUTDIR/spliceai_input_splicevardb_chr.vcf
OUTPUT_CSV=$OUTDIR/mmsplice_splicevardb_scores.csv

mkdir -p "$OUTDIR"

echo "Running MMSplice on SpliceVarDB..."
echo "GTF:    $GTF"
echo "FASTA:  $FASTA"
echo "VCF:    $VCF_ORIG"
echo "Output: $OUTPUT_CSV"

# Fix chromosome naming: Add 'chr' prefix to match FASTA format
echo "Fixing chromosome naming (adding 'chr' prefix to both header and data)..."
python - <<'FIX_CHR' "$VCF_ORIG" "$VCF_FIXED"
import sys
import re

vcf_in = sys.argv[1]
vcf_out = sys.argv[2]

with open(vcf_in, 'r') as fin, open(vcf_out, 'w') as fout:
    for line in fin:
        if line.startswith('##contig='):
            # Fix contig lines: <ID=1,length=...> -> <ID=chr1,length=...>
            line = re.sub(r'<ID=(\d+|X|Y|MT),', r'<ID=chr\1,', line)
            fout.write(line)
        elif line.startswith('#'):
            # Keep other header lines as-is
            fout.write(line)
        else:
            # Add 'chr' prefix to chromosome column in data lines
            parts = line.split('\t', 1)
            chrom = parts[0]
            # Only add prefix if not already present
            if not chrom.startswith('chr'):
                chrom = 'chr' + chrom
            fout.write(chrom + '\t' + parts[1])

print(f"Fixed VCF saved to: {vcf_out}")
FIX_CHR

# Run MMSplice following README example code
python - <<'PYTHON_CODE' "$GTF" "$FASTA" "$VCF_FIXED" "$OUTPUT_CSV"
import sys
from mmsplice.vcf_dataloader import SplicingVCFDataloader
from mmsplice import MMSplice, predict_all_table

# Get arguments
gtf, fasta, vcf, output_csv = sys.argv[1:]

print(f"Loading data from VCF: {vcf}")
# Dataloader (tissue_specific=False as per README example)
dl = SplicingVCFDataloader(gtf, fasta, vcf, tissue_specific=False)

print("Running MMSplice predictions...")
# Specify model
model = MMSplice()

# Predict and return as dataframe (with pathogenicity and splicing_efficiency)
predictions = predict_all_table(model, dl, pathogenicity=True, splicing_efficiency=True)

# Save to CSV
predictions.to_csv(output_csv, index=False)
print(f"Predictions saved to: {output_csv}")
print(f"Total predictions: {len(predictions):,}")
PYTHON_CODE

echo "MMSplice completed successfully!"
echo "Output saved to: $OUTPUT_CSV"