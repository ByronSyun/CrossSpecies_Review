# Pangolin RatGTEx Zero-Shot Evaluation

This directory contains scripts to evaluate Pangolin's zero-shot performance on the RatGTEx silver-standard dataset.

## Overview

**Pangolin** is a dilated CNN model designed for predicting splicing variants across multiple species. Unlike other models, Pangolin was actually trained on rat data (among other species), making it an interesting test case for cross-species transfer learning.

### Key Features:
- **Multi-species training**: Human + Rhesus macaque + Rat + Mouse
- **10kb context window**: Can capture long-range regulatory elements
- **Tissue-specific predictions**: Outputs scores for multiple tissues
- **Zero-shot capability**: Designed for cross-species application

## Pipeline Overview

The complete pipeline consists of:

1. **Data Conversion**: Convert RatGTEx TSV → VCF format
2. **Pangolin Prediction**: Run zero-shot predictions on rat variants
3. **Results Processing**: Parse Pangolin output and merge with labels
4. **Performance Evaluation**: Compute AUROC, AUPRC, and generate plots

## Files

### Scripts
- `convert_ratgtex_to_vcf.py`: Convert RatGTEx data to VCF format
- `predict_ratgtex_pangolin.py`: Run Pangolin predictions
- `evaluate_pangolin_results.py`: Evaluate performance metrics
- `run_pangolin_ratgtex_pipeline.sh`: Complete automated pipeline

### Expected Outputs
- `ratgtex_for_pangolin.vcf`: Input VCF for Pangolin
- `ratgtex_for_pangolin_mapping.csv`: Mapping file with labels
- `pangolin_ratgtex_results.json`: Raw prediction results
- `pangolin_performance_metrics.json`: Performance metrics
- `pangolin_ratgtex_performance_curves.png`: ROC/PR curves
- `pangolin_score_distribution.png`: Score distribution plots

## Requirements

### Server Environment (Recommended)
```bash
# Conda environment with Pangolin
/mnt/userdata4/splicing/conda_envs/pangolin-env

# Reference files
/mnt/userdata4/splicing/references/rat/ratNor6.fa
/mnt/userdata4/splicing/references/rat/ratNor6.annotation.db

# Data
/mnt/userdata4/splicing/ratgtex/processed_data/ratgtex_silver_benchmark_balanced.tsv
```

### Local Environment (Development)
```bash
# Install Pangolin
pip install pangolin

# Download rat reference genome (UCSC ratNor6)
wget http://hgdownload.cse.ucsc.edu/goldenPath/rn6/bigZips/rn6.fa.gz
gunzip rn6.fa.gz
mv rn6.fa ratNor6.fa

# Create annotation database
# (Requires rat GTF annotation file)
python -c "from pangolin.create_db import create_db; create_db('ratNor6.fa', 'ratNor6.gtf', 'ratNor6.annotation.db')"
```

## Usage

### Quick Start (Server)
```bash
# Run complete pipeline on server
./run_pangolin_ratgtex_pipeline.sh
```

### Step-by-Step (Local Development)

1. **Convert data format:**
```bash
python3 convert_ratgtex_to_vcf.py \
    --input /path/to/ratgtex_silver_benchmark_balanced.tsv \
    --output ratgtex_for_pangolin.vcf
```

2. **Run predictions:**
```bash
python3 predict_ratgtex_pangolin.py \
    --input_vcf ratgtex_for_pangolin.vcf \
    --mapping_csv ratgtex_for_pangolin_mapping.csv \
    --output_json pangolin_ratgtex_results.json
```

3. **Evaluate results:**
```bash
python3 evaluate_pangolin_results.py \
    --input_json pangolin_ratgtex_results.json \
    --output_dir ./
```

## Expected Performance

Based on Pangolin's design characteristics:

### Strengths:
- **Multi-species training**: Including rat data should provide good zero-shot performance
- **Large context window**: Can capture long-range regulatory elements (unlike SpliceAI's 10kb limit)
- **Tissue specificity**: Can identify tissue-specific effects

### Potential Limitations:
- **Training data balance**: Rat was only 28 samples vs 32 human samples
- **Limited validation**: No explicit rat performance validation in original paper
- **Cross-species annotations**: May rely on human-centric annotation features

## Technical Details

### Input Format
Pangolin expects standard VCF format:
```
##fileformat=VCFv4.2
##reference=ratNor6
#CHROM  POS     ID              REF ALT QUAL FILTER INFO
chr15   55290237 15:55290237_A/G A   G   .    PASS   LABEL=1;TISSUE=15
```

### Output Processing
Pangolin outputs tissue-specific predictions in format:
```
gene|pos:largest_increase|pos:largest_decrease|
```

We extract the maximum absolute change across tissues as the variant effect score.

### Sequence Length
- **RatGTEx data**: 8192bp context (±4096bp from variant)
- **Pangolin model**: Up to 10kb context (should be sufficient)

## Integration with Cross-Species Analysis

Results from this pipeline will be integrated with:
- AlphaGenome predictions (`../alphagenome/`)
- SpliceTransformer predictions (`../spliceTransformer/`)
- Cross-species comparison plots (`../plot_cross_species_comparison.py`)

This allows direct comparison of:
1. **Pure human models** (SpliceAI, SpliceTransformer)
2. **Multi-species augmented models** (Pangolin)
3. **Universal foundation models** (AlphaGenome, DNABERT-2, Evo2)

## Troubleshooting

### Common Issues

1. **Missing reference files**:
   - Ensure ratNor6.fa and ratNor6.annotation.db are available
   - Check file permissions and paths

2. **Memory issues**:
   - Pangolin can be memory-intensive for large datasets
   - Consider batch processing if needed

3. **Conda environment activation**:
   - Ensure pangolin-env is properly set up on server
   - Test with: `conda run -p /mnt/userdata4/splicing/conda_envs/pangolin-env pangolin --help`

4. **VCF parsing errors**:
   - Check variant ID format in RatGTEx data
   - Ensure chromosome naming consistency (chr15 vs 15)

### Performance Debugging

If performance is unexpectedly poor:
1. Check score distribution plots
2. Verify variant parsing accuracy
3. Compare with human baseline performance
4. Examine tissue-specific score patterns

## Citation

When using these results, cite:
- **Pangolin**: Zeng, H. & Li, Y.I. Predicting RNA splicing from DNA sequence using Pangolin. *Genome Biology* 23, 103 (2022).
- **RatGTEx**: (Our silver-standard dataset methodology)

## Next Steps

After running this evaluation:
1. Compare AUROC/AUPRC with other models
2. Analyze failure cases and score distributions
3. Include results in manuscript cross-species comparison
4. Consider tissue-specific analysis if relevant
