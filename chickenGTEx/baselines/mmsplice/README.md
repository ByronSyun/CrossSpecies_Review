# MMSplice Evaluation on ChickenGTEx

## Overview
MMSplice predicts variant effects on splicing using a mechanistic deep learning model.

## Pipeline

### 1. Prerequisites (Server)

**IMPORTANT**: Download Ensembl FASTA first (to match GTF chromosome naming):
```bash
cd /mnt/userdata4/splicing/ChickenGTEx/reference_genome
wget -c http://ftp.ensembl.org/pub/release-102/fasta/gallus_gallus/dna/Gallus_gallus.GRCg6a.dna.toplevel.fa.gz
gunzip Gallus_gallus.GRCg6a.dna.toplevel.fa.gz
samtools faidx Gallus_gallus.GRCg6a.dna.toplevel.fa
```

Ensure these files exist:
- **GTF**: `/mnt/userdata4/splicing/ChickenGTEx/reference_genome/Gallus_gallus.GRCg6a.gtf`
- **FASTA**: `/mnt/userdata4/splicing/ChickenGTEx/reference_genome/Gallus_gallus.GRCg6a.dna.toplevel.fa` (Ensembl)
- **FASTA index**: `.fa.fai` (created by samtools faidx)
- **Input TSV**: `/mnt/userdata4/splicing/ChickenGTEx/processed_data/chickengtex_silver_benchmark_balanced.tsv`

### 2. Upload Script to Server
```bash
scp run_mmsplice_chicken.sh yinuos@mlerp-monash-node04:/mnt/userdata4/splicing/MMSplice/
```

### 3. Run on Server
```bash
ssh yinuos@mlerp-monash-node04
cd /mnt/userdata4/splicing/MMSplice
bash run_mmsplice_chicken.sh
```

The pipeline will:
1. Check all reference files exist
2. Create FASTA index if needed
3. Create VCF from TSV (if not exists)
4. Run MMSplice predictions with pathogenicity and splicing efficiency
5. Save results to CSV

### 4. Download Results
```bash
scp yinuos@mlerp-monash-node04:/mnt/userdata4/splicing/MMSplice/chickengtex/mmsplice_chickengtex_scores.csv \
  /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/chickenGTEx/baselines/mmsplice/results/
```

## Output Format
`mmsplice_chickengtex_scores.csv` contains:
- `ID`: variant_id (format: 1:32375826_T/A)
- `delta_logit_psi`: Change in exon inclusion logit
- `pathogenicity`: Pathogenicity score (0-1)
- `efficiency`: Splicing efficiency score
- Additional columns for 5 MMSplice modules

## Model Details
- **Model**: MMSplice (VEP architecture)
- **Scores**: Delta logit PSI, pathogenicity, efficiency
- **Tissue-specific**: False (aggregated across tissues)
- **Reference genome**: GRCg6a (Gallus gallus reference genome assembly 6a)

## Expected Performance
- Coverage: ~95-100% (depends on annotation quality)
- Runtime: ~2-3 hours for 25,000 variants
- AUROC: TBD (compare with rat and pig)

## Troubleshooting

### Chromosome naming mismatch
If you see "GTF chrom names do not match with vcf chrom names":
- **Solution**: Use Ensembl FASTA (see Prerequisites above)
- FASTA must match GTF: both use Ensembl format (1, 2, 3..., Z, W, MT)
- **Do NOT use** NCBI FASTA format - it won't match GTF
- Delete old VCF if it exists: `rm /mnt/userdata4/splicing/MMSplice/chickengtex/chickengtex_input.vcf`

### Memory issues
MMSplice can be memory-intensive:
- Use a node with at least 32GB RAM
- Or split VCF into chunks and merge results

### GTF file location
If GTF is not found, download from Ensembl:
```bash
cd /mnt/userdata4/splicing/ChickenGTEx/reference_genome/
wget http://ftp.ensembl.org/pub/release-102/gtf/gallus_gallus/Gallus_gallus.GRCg6a.102.gtf.gz
gunzip Gallus_gallus.GRCg6a.102.gtf.gz
mv Gallus_gallus.GRCg6a.102.gtf Gallus_gallus.GRCg6a.gtf
```

