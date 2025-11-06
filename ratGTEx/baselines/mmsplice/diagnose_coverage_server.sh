#!/bin/bash
# Diagnose MMSplice rat coverage issue
# Run this on the server to investigate GTF quality

set -euo pipefail

echo "============================================================"
echo "MMSplice Rat Coverage Diagnostic"
echo "============================================================"
echo ""

# Activate environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /mnt/userdata4/splicing/conda_envs/mmsplice-env

# Paths
GTF_RAT=/mnt/userdata4/splicing/SpliceAI/reference_genome/Rattus_norvegicus.Rnor_6.0.84.final.gtf
GTF_HUMAN=/mnt/userdata4/splicing/SpliceAI/reference_genome/gencode.v38.annotation.gtf
FASTA_RAT=/mnt/userdata4/splicing/SpliceAI/reference_genome/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa
FASTA_HUMAN=/mnt/userdata4/splicing/Evo2Splicevardb/splicevardb_data/hg38.fa
VCF_RAT=/mnt/userdata4/splicing/SpliceAI/rat_data/spliceai_input_rat.vcf
VCF_HUMAN=/mnt/userdata4/splicing/MMSplice/splicevardb/spliceai_input_splicevardb_chr.vcf

echo "Step 1: GTF Exon Count Comparison"
echo "-----------------------------------"
RAT_EXONS=$(grep -c "^[^#].*\texon\t" "$GTF_RAT" || echo 0)
HUMAN_EXONS=$(grep -c "^[^#].*\texon\t" "$GTF_HUMAN" || echo 0)

echo "Rat GTF exons:   $RAT_EXONS"
echo "Human GTF exons: $HUMAN_EXONS"
echo "Ratio (Rat/Human): $(python3 -c "print(f'{$RAT_EXONS/$HUMAN_EXONS:.2%}')")"
echo ""

echo "Step 2: GTF Gene Count Comparison"
echo "-----------------------------------"
RAT_GENES=$(grep -c "^[^#].*\tgene\t" "$GTF_RAT" || echo 0)
HUMAN_GENES=$(grep -c "^[^#].*\tgene\t" "$GTF_HUMAN" || echo 0)

echo "Rat GTF genes:   $RAT_GENES"
echo "Human GTF genes: $HUMAN_GENES"
echo ""

echo "Step 3: VCF Variant Count"
echo "--------------------------"
RAT_VARIANTS=$(grep -cv "^#" "$VCF_RAT" || echo 0)
HUMAN_VARIANTS=$(grep -cv "^#" "$VCF_HUMAN" || echo 0)

echo "Rat VCF variants:   $RAT_VARIANTS"
echo "Human VCF variants: $HUMAN_VARIANTS"
echo ""

echo "Step 4: MMSplice Dataloader Test"
echo "----------------------------------"
python3 - <<'PYTHON_TEST'
from mmsplice.vcf_dataloader import SplicingVCFDataloader
import sys

# Test Rat
print("Testing Rat dataloader...")
try:
    gtf_rat = "/mnt/userdata4/splicing/SpliceAI/reference_genome/Rattus_norvegicus.Rnor_6.0.84.final.gtf"
    fasta_rat = "/mnt/userdata4/splicing/SpliceAI/reference_genome/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa"
    vcf_rat = "/mnt/userdata4/splicing/SpliceAI/rat_data/spliceai_input_rat.vcf"
    
    dl_rat = SplicingVCFDataloader(gtf_rat, fasta_rat, vcf_rat, tissue_specific=False)
    
    # Count loaded variants
    variant_count = 0
    for batch in dl_rat:
        variant_count += len(batch['metadata'])
    
    print(f"✓ Rat dataloader loaded successfully")
    print(f"  Variants that MMSplice will process: {variant_count}")
    
except Exception as e:
    print(f"✗ Rat dataloader failed: {e}")
    sys.exit(1)

print("")

# Test Human for comparison
print("Testing Human dataloader...")
try:
    gtf_human = "/mnt/userdata4/splicing/SpliceAI/reference_genome/gencode.v38.annotation.gtf"
    fasta_human = "/mnt/userdata4/splicing/Evo2Splicevardb/splicevardb_data/hg38.fa"
    vcf_human = "/mnt/userdata4/splicing/MMSplice/splicevardb/spliceai_input_splicevardb_chr.vcf"
    
    dl_human = SplicingVCFDataloader(gtf_human, fasta_human, vcf_human, tissue_specific=False)
    
    # Count loaded variants
    variant_count = 0
    for batch in dl_human:
        variant_count += len(batch['metadata'])
    
    print(f"✓ Human dataloader loaded successfully")
    print(f"  Variants that MMSplice will process: {variant_count}")
    
except Exception as e:
    print(f"✗ Human dataloader failed: {e}")

PYTHON_TEST

echo ""
echo "Step 5: Sample Variant Position Analysis"
echo "------------------------------------------"
echo "Checking if rat variants fall within annotated gene regions..."

python3 - <<'POSITION_CHECK'
import pandas as pd

# Parse rat GTF to get gene regions
print("Loading rat GTF gene regions...")
gene_regions = {}
with open("/mnt/userdata4/splicing/SpliceAI/reference_genome/Rattus_norvegicus.Rnor_6.0.84.final.gtf") as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 9:
            continue
        if parts[2] == "gene":
            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            if chrom not in gene_regions:
                gene_regions[chrom] = []
            gene_regions[chrom].append((start, end))

print(f"Loaded gene regions for {len(gene_regions)} chromosomes")

# Sample rat variants
print("\nChecking first 100 rat variants...")
vcf_file = "/mnt/userdata4/splicing/SpliceAI/rat_data/spliceai_input_rat.vcf"
inside_gene = 0
outside_gene = 0
checked = 0

with open(vcf_file) as f:
    for line in f:
        if line.startswith("#"):
            continue
        if checked >= 100:
            break
        
        parts = line.strip().split("\t")
        chrom = parts[0]
        pos = int(parts[1])
        
        # Check if position is within any gene region
        is_in_gene = False
        if chrom in gene_regions:
            for start, end in gene_regions[chrom]:
                if start <= pos <= end:
                    is_in_gene = True
                    break
        
        if is_in_gene:
            inside_gene += 1
        else:
            outside_gene += 1
        
        checked += 1

print(f"Results from first 100 variants:")
print(f"  Inside gene regions:  {inside_gene} ({inside_gene/checked*100:.1f}%)")
print(f"  Outside gene regions: {outside_gene} ({outside_gene/checked*100:.1f}%)")

POSITION_CHECK

echo ""
echo "============================================================"
echo "Diagnostic Complete"
echo "============================================================"
echo ""
echo "Summary:"
echo "  1. If rat exon count << human exon count → GTF quality issue"
echo "  2. If most variants are outside gene regions → explains low coverage"
echo "  3. Compare with actual MMSplice prediction count (should match Step 4)"
