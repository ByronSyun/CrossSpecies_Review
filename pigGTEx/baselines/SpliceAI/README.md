# SpliceAI Evaluation on PigGTEx

## Overview
SpliceAI evaluation on PigGTEx requires VCF format input and pig reference genome annotation.

## Pipeline Steps

### 1. Prepare Reference Files (Server)
The reference files are already available on server at:
- **Reference FASTA**: `/mnt/userdata4/splicing/PigGTEx/reference_genome/Sscrofa11.1_genomic.fna`
- **FASTA index**: `/mnt/userdata4/splicing/PigGTEx/reference_genome/Sscrofa11.1_genomic.fna.fai`
- **GTF Annotation**: Check if exists in `/mnt/userdata4/splicing/PigGTEx/reference_genome/`

If GTF is missing, download from Ensembl:
```bash
cd /mnt/userdata4/splicing/PigGTEx/reference_genome/
wget http://ftp.ensembl.org/pub/release-110/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.110.gtf.gz
gunzip Sus_scrofa.Sscrofa11.1.110.gtf.gz
```

### 2. Upload Scripts to Server
```bash
scp tsv_to_vcf_pig.py YOUR_SERVER:/mnt/userdata4/splicing/SpliceAI/
scp create_pig_annotation.py YOUR_SERVER:/mnt/userdata4/splicing/SpliceAI/
scp parse_spliceai_output.py YOUR_SERVER:/mnt/userdata4/splicing/SpliceAI/
scp run_spliceai_pig_pipeline.sh YOUR_SERVER:/mnt/userdata4/splicing/SpliceAI/
```

### 3. Run Pipeline on Server
```bash
ssh YOUR_SERVER
cd /mnt/userdata4/splicing/SpliceAI
bash run_spliceai_pig_pipeline.sh
```

The pipeline will:
1. Convert pig TSV → VCF
2. Create pig annotation file from GTF (if not exists)
3. Run SpliceAI to annotate variants
4. Parse SpliceAI output to TSV with scores

### 4. Download Results
```bash
scp YOUR_SERVER:/mnt/userdata4/splicing/SpliceAI/pig_data/spliceai_parsed_scores_pig.tsv \
  /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/pigGTEx/baselines/SpliceAI/results/
```

## Output Format
`spliceai_parsed_scores_pig.tsv`:
- Columns: CHROM, POS, REF, ALT, SYMBOL, DS_AG, DS_AL, DS_DG, DS_DL, DP_AG, DP_AL, DP_DG, DP_DL
- DS_* = delta scores (0-1, probability of splice site alteration)
- DP_* = delta positions (relative to variant)

## Model Details
- **Model**: SpliceAI v1.3
- **Distance**: 4,999 bp context window
- **Scoring**: Max delta score across all 4 splice types

## Expected Performance
- Coverage: ~95-100% (depends on annotation quality)
- AUROC: TBD (compare with rat)

