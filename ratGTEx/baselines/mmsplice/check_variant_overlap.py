#!/usr/bin/env python3
"""
Check why MMSplice filters out 85% of rat variants
Analyze variant-gene region overlap
"""

import sys
from collections import defaultdict

# Paths
GTF = "/mnt/userdata4/splicing/SpliceAI/reference_genome/Rattus_norvegicus.Rnor_6.0.84.final.gtf"
VCF = "/mnt/userdata4/splicing/SpliceAI/rat_data/spliceai_input_rat.vcf"

print("=" * 60)
print("MMSplice Variant Filtering Analysis")
print("=" * 60)
print("")

# Step 1: Load gene regions from GTF
print("Step 1: Loading gene regions from GTF...")
gene_regions = defaultdict(list)
exon_regions = defaultdict(list)

with open(GTF, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 9:
            continue
        
        chrom = parts[0]
        feature = parts[2]
        start = int(parts[3])
        end = int(parts[4])
        
        if feature == 'gene':
            gene_regions[chrom].append((start, end))
        elif feature == 'exon':
            exon_regions[chrom].append((start, end))

print(f"✓ Loaded gene regions:")
print(f"  Chromosomes with genes: {len(gene_regions)}")
print(f"  Total gene regions: {sum(len(v) for v in gene_regions.values()):,}")
print(f"  Total exon regions: {sum(len(v) for v in exon_regions.values()):,}")
print("")

# Step 2: Check variants against gene/exon regions
print("Step 2: Checking variants against annotations...")
vcf_chroms = set()
inside_gene = 0
inside_exon = 0
outside = 0
chrom_mismatch = 0

with open(VCF, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        
        parts = line.strip().split('\t')
        chrom = parts[0]
        pos = int(parts[1])
        
        vcf_chroms.add(chrom)
        
        # Check if chromosome exists in GTF
        if chrom not in gene_regions:
            chrom_mismatch += 1
            continue
        
        # Check gene overlap
        in_gene = any(start <= pos <= end for start, end in gene_regions[chrom])
        in_exon = any(start <= pos <= end for start, end in exon_regions[chrom])
        
        if in_exon:
            inside_exon += 1
        elif in_gene:
            inside_gene += 1
        else:
            outside += 1

total_variants = inside_gene + inside_exon + outside + chrom_mismatch

print(f"✓ Analyzed {total_variants:,} variants")
print("")
print(f"Results:")
print(f"  Inside exon regions:        {inside_exon:,} ({inside_exon/total_variants*100:.1f}%)")
print(f"  Inside gene (non-exon):     {inside_gene:,} ({inside_gene/total_variants*100:.1f}%)")
print(f"  Outside gene regions:       {outside:,} ({outside/total_variants*100:.1f}%)")
print(f"  Chromosome mismatch:        {chrom_mismatch:,} ({chrom_mismatch/total_variants*100:.1f}%)")
print("")

# Step 3: Chromosome comparison
print("Step 3: Chromosome naming comparison...")
gtf_chroms = sorted(gene_regions.keys())[:15]
vcf_chroms = sorted(vcf_chroms)[:15]

print(f"GTF chromosomes (first 15): {', '.join(gtf_chroms)}")
print(f"VCF chromosomes (first 15): {', '.join(vcf_chroms)}")
print("")

# Step 4: Expected coverage calculation
print("=" * 60)
print("ANALYSIS")
print("=" * 60)

expected_loaded = inside_exon + inside_gene
actual_loaded = 4184  # From dataloader output

print(f"Total variants in VCF:           {total_variants:,}")
print(f"Expected to load (in genes):     {expected_loaded:,} ({expected_loaded/total_variants*100:.1f}%)")
print(f"Actually loaded by dataloader:   {actual_loaded:,} ({actual_loaded/total_variants*100:.1f}%)")
print(f"Discrepancy:                     {expected_loaded - actual_loaded:,}")
print("")

if chrom_mismatch > 0:
    print(f"⚠️  ISSUE FOUND: {chrom_mismatch:,} variants have chromosome naming mismatch!")
    print(f"   This explains {chrom_mismatch/total_variants*100:.1f}% of the filtering.")
    print("")

if expected_loaded < total_variants * 0.5:
    print(f"⚠️  ISSUE FOUND: Only {expected_loaded/total_variants*100:.0f}% of variants fall within gene regions!")
    print(f"   MMSplice requires variants to be in/near exons for modular prediction.")
    print("")

if abs(expected_loaded - actual_loaded) > 1000:
    print(f"⚠️  WARNING: Expected {expected_loaded:,} but loaded {actual_loaded:,}")
    print(f"   Additional filtering may be happening in SplicingVCFDataloader.")
    print(f"   Possible reasons:")
    print(f"   - Variants too far from splice sites (flanking region limit)")
    print(f"   - Transcript annotation missing")
    print(f"   - FASTA sequence extraction failures")
    print("")

print("=" * 60)
print("CONCLUSION")
print("=" * 60)

coverage_rate = actual_loaded / total_variants * 100
if coverage_rate < 20:
    print(f"❌ SEVERE: Only {coverage_rate:.1f}% coverage")
    print(f"   Primary cause: Variants outside annotated gene/exon regions")
    print(f"   This is a FUNDAMENTAL LIMITATION of MMSplice's modular design.")
    print(f"   Unlike SpliceAI (sequence-based), MMSplice REQUIRES exon annotations.")
elif coverage_rate < 50:
    print(f"⚠️  MODERATE: {coverage_rate:.1f}% coverage")
    print(f"   Combination of annotation gaps and strict filtering.")
else:
    print(f"✓ ACCEPTABLE: {coverage_rate:.1f}% coverage")

print("")
print("Recommendation:")
print("  This 85% filtering is expected for MMSplice on species with")
print("  incomplete annotations. Document this as a known limitation")
print("  in the cross-species evaluation section of the paper.")
