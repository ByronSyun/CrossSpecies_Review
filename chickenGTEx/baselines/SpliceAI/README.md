# SpliceAI Chicken Baseline

This directory contains scripts to run SpliceAI predictions on the ChickenGTEx silver-standard benchmark.

## Overview

SpliceAI is a deep neural network that predicts splice junctions from pre-mRNA sequences. This pipeline adapts SpliceAI for chicken (Gallus gallus) data using the GRCg6a reference genome.

## Pipeline Steps

1. **TSV → VCF Conversion**: Convert ChickenGTEx benchmark to VCF format
2. **Annotation Creation**: Generate SpliceAI-compatible annotation from chicken GTF
3. **SpliceAI Prediction**: Run SpliceAI on chicken variants
4. **Result Parsing**: Extract delta scores from SpliceAI output

## Files

- `tsv_to_vcf_chicken.py`: Convert TSV to VCF format
- `create_chicken_annotation.py`: Generate SpliceAI annotation from GTF
- `parse_spliceai_output.py`: Parse SpliceAI VCF output to TSV
- `run_spliceai_chicken_pipeline.sh`: Complete server-side pipeline

## Requirements

### Server Environment

- **SpliceAI**: Installed in `/mnt/userdata4/splicing/conda_envs/spliceai_env`
- **Python packages**: pandas, gffutils
- **Reference files**:
  - FASTA: `Gallus_gallus.GRCg6a.dna.toplevel.fa`
  - GTF: `Gallus_gallus.GRCg6a.102.gtf`

## Usage

### On Server

```bash
# 1. Upload scripts to server
scp -r * your_server:/mnt/userdata4/splicing/SpliceAI/

# 2. SSH to server
ssh your_server

# 3. Run complete pipeline
cd /mnt/userdata4/splicing/SpliceAI/
bash run_spliceai_chicken_pipeline.sh
```

### Pipeline Output

```
/mnt/userdata4/splicing/SpliceAI/chicken_data/
├── spliceai_input_chicken.vcf           # Input VCF (25,000 variants)
├── chicken_grcg6a_annotation.txt        # SpliceAI annotation file
├── spliceai_predictions_chicken.vcf     # SpliceAI predictions
└── spliceai_parsed_scores_chicken.tsv   # Parsed scores (CHROM, POS, REF, ALT, MAX_DS)
```

## Data Format

### Input TSV (ChickenGTEx)
```
variant_id    ref_seq    alt_seq    label    tissue_id
1:15450992_C/T    ATGC...    ATGT...    1    liver
```

### Output TSV (SpliceAI Scores)
```
CHROM    POS          ID    REF    ALT    GENE      DS_AG    DS_AL    DS_DG    DS_DL    MAX_DS
1        15450992     .     C      T      ENSGALG... 0.12     0.05     0.89     0.03     0.89
```

## Reference Genome

- **Assembly**: GRCg6a (Gallus gallus 6a)
- **Source**: Ensembl Release 102
- **Chromosomes**: 1-28, W, Z, MT
- **FASTA**: ~1.2 GB
- **GTF**: ~50 MB

## SpliceAI Parameters

```bash
spliceai -I input.vcf \
         -O output.vcf \
         -R reference.fa \
         -A annotation.txt \
         -D 4999  # Distance parameter (default 50, max 10000)
```

### Distance Parameter (-D)
- Default: 50 (looks ±50bp around splice sites)
- Used: 4999 (looks ±5kb for long-range effects)
- Max: 10000

## Performance Notes

- **Runtime**: ~2-4 hours for 25,000 variants
- **Memory**: ~8-16 GB RAM
- **Disk**: ~500 MB for outputs
- **GPU**: Optional (speeds up prediction)

## Troubleshooting

### Issue: VCF conversion fails
```bash
# Check variant ID format
head chickengtex_silver_benchmark_balanced.tsv | cut -f1

# Expected format: '1:15450992_C/T'
```

### Issue: Annotation creation hangs
```bash
# Check GTF file
head -50 Gallus_gallus.GRCg6a.102.gtf

# Remove existing database and retry
rm Gallus_gallus.GRCg6a.102.gtf.db
python create_chicken_annotation.py --input_gtf ... --output_txt ...
```

### Issue: SpliceAI prediction fails
```bash
# Check FASTA index exists
ls -lh Gallus_gallus.GRCg6a.dna.toplevel.fa.fai

# Create index if missing
samtools faidx Gallus_gallus.GRCg6a.dna.toplevel.fa

# Check chromosome naming consistency
# VCF should use: 1, 2, 3, ... (not chr1, chr2, ...)
# FASTA and GTF should match
```

### Issue: No SpliceAI annotations in output
```bash
# Check if variants are within annotated genes
# SpliceAI only predicts for variants near genes

# Verify annotation file format
head chicken_grcg6a_annotation.txt

# Expected format:
# #NAME  CHROM  STRAND  TX_START  TX_END  EXON_START  EXON_END
```

## Expected Results

Based on cross-species evaluation patterns:

- **Coverage**: ~50-60% (lower than human due to annotation differences)
- **AUROC**: ~0.55-0.65 (performance degradation from human)
- **AUPRC**: ~0.60-0.70

Performance degradation expected due to:
1. Human-centric training of SpliceAI
2. Chicken-specific regulatory elements
3. Annotation quality differences

## Citations

```bibtex
@article{jaganathan2019spliceai,
  title={Predicting splicing from primary sequence with deep learning},
  author={Jaganathan, Kishore and others},
  journal={Cell},
  year={2019}
}
```

## Contact

For issues or questions about this chicken-specific implementation, refer to the main BIB_review repository documentation.

