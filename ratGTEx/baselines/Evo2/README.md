# Evo2 Cross-Species Evaluation: Rat (ratGTEx)

This directory contains scripts for evaluating Evo2 (7B parameter, 8kb context window) on rat genomic variants from the ratGTEx silver-standard benchmark.

## Overview

Evo2 is a large-scale genomic foundation model pre-trained on sequences from 1000+ species. Unlike task-specific models (SpliceAI, Pangolin) or specialized foundation models (SpliceBERT), Evo2 uses a sequence-driven approach with delta log-probability ($\Delta \log p$) scoring to quantify mutation effects without requiring task-specific fine-tuning.

### Hardcode Bypass Strategy: Sequence-Driven Prediction

Evo2's architecture enables direct cross-species application through its sequence-driven paradigm:
- **No coordinate dependencies**: Operates purely on DNA sequences without genomic annotations
- **Species-agnostic tokenization**: Character-level tokenizer processes any nucleotide sequence
- **Universal pre-training**: Trained on diverse genomes, capturing fundamental sequence grammar

This makes Evo2 inherently cross-species compatible, requiring only properly formatted sequence windows as input.

## Directory Structure

```
Evo2/
├── README.md                        # This file
├── predict_evo2_rat.py              # Prediction script (run on server)
├── evaluate_evo2_rat.py             # Evaluation script (run locally)
├── run_evo2_rat_server.sh           # Server execution wrapper
└── results/                         # Output directory (created automatically)
    ├── evo2_rat_predictions.tsv     # Prediction results
    ├── evo2_rat_results.json        # Performance metrics
    ├── evo2_rat_curves.png          # ROC/PR curves
    └── evo2_rat_score_distribution.png  # Score distributions
```

## Input Data Format

The input file (`ratgtex_silver_benchmark_balanced.tsv`) contains:
- **Format**: Space-separated, no header
- **Columns**: 
  1. `variant_id` (e.g., "15:55290237_A/G")
  2. `ref_seq` (8192bp reference sequence centered on variant)
  3. `alt_seq` (8192bp alternative sequence centered on variant)
  4. `label` (0=normal, 1=splice-altering)
  5. `chrom` (chromosome identifier)
- **Size**: 28,120 variants (50/50 balanced)

## Server Setup

### Prerequisites

1. **Server paths** (already configured):
   - Working directory: `/mnt/userdata4/splicing/Evo2_rat`
   - Base project: `/mnt/userdata4/splicing/Evo2Splicevardb`
   - Container: `vortex.sif` (includes Evo2 model and dependencies)

2. **Files to upload to server**:
   ```bash
   # From local machine
   cd /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/ratGTEx/baselines/Evo2
   
   # Upload prediction script
   scp predict_evo2_rat.py user@server:/mnt/userdata4/splicing/Evo2_rat/
   
   # Upload run script
   scp run_evo2_rat_server.sh user@server:/mnt/userdata4/splicing/Evo2_rat/
   
   # Upload input data (if not already there)
   scp /Users/byronsun/Desktop/AS_复现模型/BIB_review/data/processed_data/ratGTEx/ratgtex_silver_benchmark_balanced.tsv \
       user@server:/mnt/userdata4/splicing/Evo2_rat/
   ```

## Execution Workflow

### Step 1: Run Prediction on Server

```bash
# SSH to server
ssh user@server

# Navigate to working directory
cd /mnt/userdata4/splicing/Evo2_rat

# Make script executable
chmod +x run_evo2_rat_server.sh

# Run prediction (expect ~3-6 hours for 28k variants with batch_size=1)
bash run_evo2_rat_server.sh
```

**Expected runtime**: 3-6 hours (depends on GPU availability)
- Batch size: 1 (conservative for 8kb sequences to avoid OOM)
- GPU memory: ~40GB for Evo2-7B model

### Step 2: Download Results

```bash
# From local machine
cd /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/ratGTEx/baselines/Evo2
mkdir -p results

# Download predictions
scp user@server:/mnt/userdata4/splicing/Evo2_rat/evo2_rat_predictions.tsv results/
```

### Step 3: Run Evaluation Locally

```bash
# On local machine
cd /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/ratGTEx/baselines/Evo2

python3 evaluate_evo2_rat.py \
    --predictions results/evo2_rat_predictions.tsv \
    --output_dir results
```

**Output files**:
- `evo2_rat_results.json`: Performance metrics (AUROC, AUPRC, etc.)
- `evo2_rat_curves.png`: ROC and PR curves
- `evo2_rat_score_distribution.png`: Score distributions by class

## Expected Results

### Hypothesis

Based on Evo2's performance on human SpliceVarDB (reported in draft.tex), we expect:

**Human performance** (from literature/previous evaluation):
- AUROC: ~0.65-0.75 (sequence-driven, no task-specific training)
- AUPRC: ~0.70-0.75

**Rat performance (hypothesis)**:
- **Scenario 1 (Strong cross-species transfer)**: AUROC 0.60-0.70
  - Suggests Evo2's multi-species pre-training captures universal splicing patterns
  - Would represent the best cross-species transfer among evaluated models

- **Scenario 2 (Moderate degradation)**: AUROC 0.55-0.60
  - Indicates some species-specific limitations despite broad pre-training
  - Still better than complete failure (AUROC ~0.50)

- **Scenario 3 (Severe degradation)**: AUROC 0.50-0.55
  - Would align with other models' failure patterns (SpliceAI, Pangolin, SpliceTransformer)
  - Suggests fundamental cross-species prediction gap persists even for large-scale GFMs

### Scientific Significance

The rat evaluation will test whether:
1. **Scale matters**: Does 1000+ species pre-training overcome the cross-species gap?
2. **Architecture matters**: Do sequence-driven models generalize better than coordinate-driven ones?
3. **Task alignment matters**: Can general-purpose pre-training substitute for splicing-specific training?

## Comparison with Other Models

### Cross-Species Performance on ratGTEx

| Model               | Architecture | Species Training | Human AUROC | Rat AUROC | Performance Drop |
|:--------------------|:-------------|:-----------------|:------------|:----------|:-----------------|
| Pangolin            | CNN          | Human+3          | 0.9604      | 0.5100    | -46.9%           |
| SpliceAI            | CNN          | Human            | 0.9583      | 0.5696    | -40.6%           |
| SpliceTransformer   | Transformer  | Human+4 partial  | 0.9431      | 0.5011    | -46.9%           |
| SpliceBERT          | Transformer  | 72 vertebrates   | 0.5163      | 0.5044    | -2.3%            |
| Nucleotide Transformer | Transformer | 850+ species | 0.5073      | 0.5063    | -0.2%            |
| **Evo2**            | Hyena        | 1000+ species    | **TBD**     | **TBD**   | **TBD**          |

### Key Questions

1. **Will Evo2 outperform SpliceBERT and NT on human data?**
   - SpliceBERT (splicing-specialized, 72 species): AUROC 0.5163
   - NT (general-purpose, 850+ species): AUROC 0.5073
   - Evo2 (general-purpose, 1000+ species, longer context): ?

2. **Will Evo2 maintain performance on rat data?**
   - If yes: Strongest evidence yet for cross-species transferability via foundation models
   - If no: Further confirmation that scale alone cannot overcome architectural limitations

## Troubleshooting

### Common Issues

1. **Out of Memory (OOM) on GPU**
   - Solution: Reduce batch size to 1 in `run_evo2_rat_server.sh`
   - Alternative: Use CPU (very slow, not recommended)

2. **Container not found**
   - Check that `vortex.sif` exists in base project directory
   - Verify apptainer is available: `which apptainer`

3. **Tokenization errors**
   - Ensure sequences contain only A, C, G, T, N characters
   - Check input file format (space-separated, no header)

4. **Slow prediction speed**
   - Expected: ~1-2 seconds per variant with batch_size=1
   - Monitor GPU usage: `nvidia-smi` (should show ~40GB memory)

## References

- **Evo2 paper**: Brandes et al. (2022). Genome-wide prediction of disease-variant function. bioRxiv.
- **SpliceVarDB**: Sullivan et al. (2024). Am J Hum Genet 111:2164–2175.
- **ratGTEx**: Derived from FarmGTEx project (www.farmgtex.org)

## Citation

If you use these scripts, please cite:
- Evo2 model: Brandes et al. (2022)
- This review: Sun Y. et al. (2025). Cross-Species Transfer Learning for Splicing Variant Prediction. Brief Bioinform (in preparation)

