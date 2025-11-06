# MMSplice Cross-Species Evaluation: Rat (ratGTEx)

## Overview

This directory contains scripts to evaluate MMSplice's cross-species transfer performance on rat genetic variants from the ratGTEx dataset. Following the successful human (SpliceVarDB) evaluation where MMSplice achieved AUROC 0.9230 (pathogenicity) and 0.9011 (delta_logit_psi magnitude), this evaluation tests whether these human-trained modules can generalize to rat genomic contexts.

## Background

**MMSplice Performance on Human:**
- Pathogenicity score: AUROC 0.9230, AUPRC 0.9429
- Delta_logit_psi (absolute): AUROC 0.9011, AUPRC 0.9317  
- 23,588 variants evaluated

**Cross-Species Challenge:**
MMSplice's modular design, trained on human GTEx RNA-seq and experimental data (MaPSy), embeds human-specific splicing regulatory patterns. The model assumes hg38/hg19 reference genomes and GENCODE/Ensembl human annotations, creating barriers for cross-species application similar to those observed with:
- SpliceAI: AUROC 0.9583 (human) → 0.5696 (rat)
- SpliceTransformer: AUROC 0.9431 (human) → 0.5011 (rat)
- Pangolin: AUROC 0.9604 (human) → 0.5100 (rat)

## Data Requirements

### Server Paths (Already Available)

```bash
# Rat reference files (Ensembl Rnor_6.0, rn6 assembly)
GTF=/mnt/userdata4/splicing/SpliceAI/reference_genome/Rattus_norvegicus.Rnor_6.0.84.final.gtf
FASTA=/mnt/userdata4/splicing/SpliceAI/reference_genome/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa
VCF=/mnt/userdata4/splicing/SpliceAI/rat_data/spliceai_input_rat.vcf

# FASTA index (will be created if not exists)
FASTA_INDEX=${FASTA}.fai
```

### Chromosome Naming Format

**Critical:** Rat genome (rn6/Rnor_6.0) uses **numerical chromosome naming** (1, 2, 3..., X) **without 'chr' prefix** (Ensembl standard), unlike human hg38 which uses chr1, chr2, chr3.

Expected format consistency:
- **VCF**: `1`, `2`, `3` (confirmed)
- **FASTA**: `>1`, `>2`, `>3` (Ensembl standard)
- **GTF**: `1`, `2`, `3` (Ensembl standard)

If chromosome naming mismatch occurs ("Fasta chrom names do not match with vcf chrom names"), the VCF would need conversion. However, since all Ensembl rat files use numerical format, they should already match.

## Usage

### 1. Server Execution

```bash
# On server (with MMSplice environment activated)
bash /mnt/userdata4/splicing/MMSplice/ratgtex/run_mmsplice_rat.sh
```

**Expected Runtime:** ~3-20 hours for 28,120 variants (depending on server load)

**Output:**
- `/mnt/userdata4/splicing/MMSplice/ratgtex/mmsplice_ratgtex_scores.csv`

### 2. Local Evaluation

```bash
# Download results from server
scp user@server:/mnt/userdata4/splicing/MMSplice/ratgtex/mmsplice_ratgtex_scores.csv \
    /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/ratGTEx/baselines/mmsplice/results/

# Run evaluation locally
cd /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/ratGTEx/baselines/mmsplice
python3 evaluate_mmsplice_rat.py
```

**Outputs:**
- `results/mmsplice_ratgtex_merged.csv`: Merged predictions and labels
- `results/mmsplice_rat_pathogenicity_curves.png`: ROC/PR curves
- `results/mmsplice_rat_pathogenicity_results.json`: Detailed metrics
- `results/mmsplice_rat_delta_logit_psi_abs_curves.png`
- `results/mmsplice_rat_delta_logit_psi_abs_results.json`
- `results/mmsplice_rat_efficiency_abs_curves.png`
- `results/mmsplice_rat_efficiency_abs_results.json`
- `results/mmsplice_rat_scores_comparison.json`: Summary comparison

## Key Technical Details

### Score Interpretation

Following the human evaluation findings, we use:

1. **Pathogenicity score** (primary): 
   - Pre-calibrated during training
   - Direct interpretation: higher = more pathogenic
   - Human AUROC: 0.9230

2. **Delta_logit_psi (absolute value)**:
   - Captures magnitude of splicing change
   - Direction-magnitude trade-off: negative values (exon exclusion) more common for splice-altering
   - Using absolute value improves AUROC from 0.08 → 0.90 on human
   - Human AUROC (abs): 0.9011

3. **Efficiency (absolute value)**:
   - Similar direction-magnitude pattern
   - Human AUROC (abs): 0.8960

### Aggregation Strategy

MMSplice outputs multiple predictions per variant (one per affected exon/transcript). Following best practices:
- **Aggregation method**: Maximum absolute effect across all exons
- **Rationale**: Most severe disruption likely determines pathogenicity

## Expected Outcomes

### Hypothesis 1: Human-Specific Degradation
If MMSplice's modular components are human-specific:
- **Expected:** AUROC ~0.50-0.55 (near-random), similar to SpliceAI/SpliceTransformer
- **Interpretation:** Module weights trained on human GTEx don't capture rat splicing code

### Hypothesis 2: Partial Transferability
If core splicing signals are conserved:
- **Expected:** AUROC 0.60-0.75 (moderate)
- **Interpretation:** Basic splice site recognition transfers, but context-specific regulation doesn't

### Hypothesis 3: Robust Generalization (Unlikely)
If modular design enables cross-species application:
- **Expected:** AUROC >0.80
- **Interpretation:** Modular decomposition captures universal splicing principles

## Comparison Context

### Human Performance (SpliceVarDB)
| Model | AUROC | AUPRC | Notes |
|-------|-------|-------|-------|
| Pangolin | 0.9604 | 0.9687 | Best human performance |
| SpliceAI | 0.9583 | 0.9652 | Clinical gold standard |
| SpliceTransformer | 0.9431 | 0.9625 | Tissue-aware |
| **MMSplice (pathogenicity)** | **0.9230** | **0.9429** | **Modular, interpretable** |
| MMSplice (delta_logit_psi, abs) | 0.9011 | 0.9317 | Magnitude-based |

### Rat Cross-Species Performance
| Model | AUROC | Notes |
|-------|-------|-------|
| SpliceAI | 0.5696 | "Genome substitution" bypass |
| Pangolin | 0.5100 | Multi-species training |
| SpliceTransformer | 0.5011 | "Reference approximation" bypass |
| **MMSplice** | **TBD** | **Modular human-trained** |

## Troubleshooting

### Common Issues

**1. Chromosome Naming Mismatch**
```
ValueError: Fasta chrom names do not match with vcf chrom names
```
**Solution:** Check if FASTA uses 'chr1' while VCF uses '1'. Rat Ensembl files should all use numerical format, so this error suggests file inconsistency. Verify:
```bash
# Check formats
grep "^>" $FASTA | head -3
grep -v "^#" $VCF | head -3 | cut -f1
grep -v "^#" $GTF | cut -f1 | sort -u | head -3
```

**2. Missing FASTA Index**
```
ERROR: .fai file not found
```
**Solution:** Script auto-creates it, but you can manually:
```bash
samtools faidx /mnt/userdata4/splicing/SpliceAI/reference_genome/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa
```

**3. Memory Issues**
If prediction fails with memory error, consider:
- Run on a node with more RAM
- Process variants in batches (requires script modification)

## Files

```
mmsplice/
├── README.md                      # This file
├── run_mmsplice_rat.sh           # Server execution script
├── evaluate_mmsplice_rat.py      # Local evaluation script
└── results/
    ├── mmsplice_ratgtex_scores.csv          # Raw predictions (from server)
    ├── mmsplice_ratgtex_merged.csv          # Merged with labels
    ├── mmsplice_rat_pathogenicity_curves.png
    ├── mmsplice_rat_pathogenicity_results.json
    ├── mmsplice_rat_delta_logit_psi_abs_curves.png
    ├── mmsplice_rat_delta_logit_psi_abs_results.json
    ├── mmsplice_rat_efficiency_abs_curves.png
    ├── mmsplice_rat_efficiency_abs_results.json
    └── mmsplice_rat_scores_comparison.json  # Summary
```

## References

- **MMSplice Paper:** Cheng J, et al. MMSplice: modular modeling of RNA splicing. *Genome Biol* 2019;20:46.
- **Human Evaluation:** SpliceVarDB benchmark (23,588 variants, AUROC 0.9230)
- **Rat Dataset:** ratGTEx silver-standard (28,120 variants, balanced)
