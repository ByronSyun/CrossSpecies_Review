# MMSplice Evaluation on PigGTEx

## Overview
MMSplice predicts variant effects on splicing using a mechanistic deep learning model.

## Pipeline

### 1. Prerequisites (Server)

**IMPORTANT**: Download Ensembl FASTA first (to match GTF chromosome naming):
```bash
cd /mnt/userdata4/splicing/PigGTEx/reference_genome
wget -c http://ftp.ensembl.org/pub/release-110/fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz
gunzip Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz
samtools faidx Sus_scrofa.Sscrofa11.1.dna.toplevel.fa
```

Ensure these files exist:
- **GTF**: `/mnt/userdata4/splicing/PigGTEx/reference_genome/Sus_scrofa.Sscrofa11.1.gtf`
- **FASTA**: `/mnt/userdata4/splicing/PigGTEx/reference_genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa` (Ensembl)
- **FASTA index**: `.fa.fai` (created by samtools faidx)
- **Input TSV**: `/mnt/userdata4/splicing/PigGTEx/processed_data/piggtex_silver_benchmark_balanced.tsv`

### 2. Upload Script to Server
```bash
scp run_mmsplice_pig.sh YOUR_SERVER:/mnt/userdata4/splicing/MMSplice/
```

### 3. Run on Server
```bash
ssh YOUR_SERVER
cd /mnt/userdata4/splicing/MMSplice
bash run_mmsplice_pig.sh
```

The pipeline will:
1. Check all reference files exist
2. Create FASTA index if needed
3. Create VCF from TSV (if not exists)
4. Run MMSplice predictions with pathogenicity and splicing efficiency
5. Save results to CSV

### 4. Download Results
```bash
scp YOUR_SERVER:/mnt/userdata4/splicing/MMSplice/piggtex/mmsplice_piggtex_scores.csv \
  /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/pigGTEx/baselines/mmsplice/results/
```

## Output Format
`mmsplice_piggtex_scores.csv` contains:
- `ID`: variant_id (format: 13:65362642_C/T)
- `delta_logit_psi`: Change in exon inclusion logit
- `pathogenicity`: Pathogenicity score (0-1)
- `efficiency`: Splicing efficiency score
- Additional columns for 5 MMSplice modules

## Model Details
- **Model**: MMSplice (VEP architecture)
- **Scores**: Delta logit PSI, pathogenicity, efficiency
- **Tissue-specific**: False (aggregated across tissues)
- **Reference genome**: Sscrofa11.1

## Expected Performance
- Coverage: ~95-100% (depends on annotation quality)
- Runtime: ~2-4 hours for 26,358 variants
- AUROC: TBD (compare with rat)

## Troubleshooting

### Chromosome naming mismatch
If you see "GTF chrom names do not match with vcf chrom names":
- **Solution**: Use Ensembl FASTA (see Prerequisites above)
- FASTA must match GTF: both use Ensembl format (1, 2, 3..., X, Y)
- **Do NOT use** RefSeq FASTA (NC_XXXXXX format) - it won't match GTF
- Delete old VCF if it exists: `rm /mnt/userdata4/splicing/MMSplice/piggtex/piggtex_input.vcf`

### Memory issues
MMSplice can be memory-intensive:
- Use a node with at least 32GB RAM
- Or split VCF into chunks and merge results

