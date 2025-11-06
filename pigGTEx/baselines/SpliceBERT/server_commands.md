# SpliceBERT PigGTEx Server Commands

## Step 1: Upload Script to Server

From your local machine:

```bash
scp /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/pigGTEx/baselines/SpliceBERT/splicebert_score_pig_variants.py \
    yinuos@mlerp-monash-node04:/mnt/userdata4/splicing/SpliceBERT/pig_analysis/
```

## Step 2: Setup Environment

On the server:

```bash
# Create working directory
mkdir -p /mnt/userdata4/splicing/SpliceBERT/pig_analysis/results

# Activate SpliceBERT environment
conda activate /mnt/userdata4/splicing/conda_envs/splicebert-env

# Navigate to working directory
cd /mnt/userdata4/splicing/SpliceBERT/pig_analysis
```

## Step 3: Run SpliceBERT Scoring

```bash
# Run SpliceBERT variant scoring on pig data
python splicebert_score_pig_variants.py \
    --input_tsv /mnt/userdata4/splicing/PigGTEx/processed_data/piggtex_silver_benchmark_balanced.tsv \
    --model_path /mnt/userdata4/splicing/SpliceBERT/SpliceBERT-main/models/models/SpliceBERT.1024nt \
    --output_prefix results/splicebert_pig \
    --flanking_window 100 \
    --batch_size 4
```

**Expected output files:**
- `results/splicebert_pig_scores.jsonl`
- `results/splicebert_pig_scores.csv`

**Expected runtime:** 2-4 hours for ~26,358 variants

## Step 4: Download Results

To analyze results locally:

```bash
# From your local machine, download the results
scp -r yinuos@mlerp-monash-node04:/mnt/userdata4/splicing/SpliceBERT/pig_analysis/results/* \
    /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/pigGTEx/baselines/SpliceBERT/results/
```

## Expected Performance

Based on SpliceBERT's cross-species characteristics:
- **Coverage**: Should achieve 100% coverage (no annotation dependency)
- **RNA-centric advantage**: Uses RNA sequences, which may be more conserved across species
- **Zero-shot**: No fine-tuning on pig data

## Troubleshooting

### Common Issues:

1. **CUDA memory error**: Reduce batch size to 2 or 1
2. **Model loading error**: Verify model path exists at `/mnt/userdata4/splicing/SpliceBERT/SpliceBERT-main/models/models/SpliceBERT.1024nt`
3. **Tokenization issues**: Ensure transformers library is compatible

### Resource Requirements:

- **GPU memory**: ~8-12GB for batch size 4
- **Time**: ~2-4 hours for 26,358 variants
- **Disk space**: ~500MB for results

## Notes

- SpliceBERT requires RNA format (T->U conversion) and space-separated tokens
- The model uses KL-divergence scoring in a ±100nt flanking window
- Lower KL-context scores typically indicate more splice-altering variants (inverted scoring)
- This is a zero-shot cross-species evaluation (no fine-tuning on pig data)

