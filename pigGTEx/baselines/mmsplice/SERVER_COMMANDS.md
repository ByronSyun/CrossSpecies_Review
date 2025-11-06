# MMSplice PigGTEx Server Commands

## Upload Script
```bash
scp run_mmsplice_pig.sh YOUR_SERVER:/mnt/userdata4/splicing/MMSplice/
```

## Run on Server
```bash
ssh YOUR_SERVER
cd /mnt/userdata4/splicing/MMSplice
bash run_mmsplice_pig.sh
```

## Download Results
```bash
# Main results
scp YOUR_SERVER:/mnt/userdata4/splicing/MMSplice/piggtex/mmsplice_piggtex_scores.csv \
  /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/pigGTEx/baselines/mmsplice/results/

# VCF (optional, for debugging)
scp YOUR_SERVER:/mnt/userdata4/splicing/MMSplice/piggtex/piggtex_input.vcf \
  /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/pigGTEx/baselines/mmsplice/results/
```

## Expected Runtime
- VCF creation: ~1 minute
- MMSplice predictions: ~2-4 hours (26,358 variants)
- Total: ~2-4 hours

## Expected Output
`mmsplice_piggtex_scores.csv`:
- Rows: ~26,358 (one per variant)
- Key columns: ID, delta_logit_psi, pathogenicity, efficiency
- Size: ~5-10 MB

## Troubleshooting

### If GTF not found
```bash
# Download GTF (if not already done for SpliceAI)
cd /mnt/userdata4/splicing/PigGTEx/reference_genome/
wget http://ftp.ensembl.org/pub/release-110/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.110.gtf.gz
gunzip Sus_scrofa.Sscrofa11.1.110.gtf.gz
mv Sus_scrofa.Sscrofa11.1.110.gtf Sus_scrofa.Sscrofa11.1.gtf
```

### If FASTA index missing
The script will create it automatically, or manually:
```bash
samtools faidx /mnt/userdata4/splicing/PigGTEx/reference_genome/Sscrofa11.1_genomic.fna
```

