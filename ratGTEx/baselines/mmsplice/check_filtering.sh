#!/bin/bash
# Check why 28120 variants → 4184 loaded
# Run on server

set -euo pipefail

GTF=/mnt/userdata4/splicing/SpliceAI/reference_genome/Rattus_norvegicus.Rnor_6.0.84.final.gtf
VCF=/mnt/userdata4/splicing/SpliceAI/rat_data/spliceai_input_rat.vcf

echo "============================================================"
echo "MMSplice Rat Filtering Analysis"
echo "============================================================"
echo ""

echo "1. GTF Statistics"
echo "-----------------"
echo -n "Total genes: "
awk '$3=="gene"' "$GTF" | wc -l
echo -n "Total exons: "
awk '$3=="exon"' "$GTF" | wc -l
echo -n "Total transcripts: "
awk '$3=="transcript"' "$GTF" | wc -l

echo ""
echo "2. VCF Statistics"
echo "-----------------"
echo -n "Total variants: "
grep -cv "^#" "$VCF"

echo ""
echo "3. Sample GTF (first 3 genes)"
echo "------------------------------"
awk '$3=="gene"' "$GTF" | head -3 | cut -f1,3,4,5,9 | column -t

echo ""
echo "4. Sample VCF (first 5 variants)"
echo "---------------------------------"
grep -v "^#" "$VCF" | head -5 | cut -f1-5

echo ""
echo "5. Chromosome Distribution"
echo "---------------------------"
echo "GTF chromosomes:"
awk '$3=="gene"' "$GTF" | cut -f1 | sort -u | head -10

echo ""
echo "VCF chromosomes:"
grep -v "^#" "$VCF" | cut -f1 | sort -u | head -10

echo ""
echo "============================================================"
echo "Run Python script for detailed overlap analysis..."
echo "============================================================"
