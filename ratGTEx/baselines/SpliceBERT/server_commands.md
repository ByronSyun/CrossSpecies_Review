
# Activate SpliceBERT environment
conda activate /mnt/userdata4/splicing/conda_envs/splicebert-env

# Navigate to working directory
cd /mnt/userdata4/splicing/SpliceBERT/rat_analysis
```

## Step 3: Run SpliceBERT Scoring

```bash
# Run SpliceBERT variant scoring on rat data
python splicebert_score_rat_variants.py \
    --input_tsv /mnt/userdata4/splicing/ratgetx/processed_data/ratgtex_silver_benchmark_balanced.tsv \
    --model_path /mnt/userdata4/splicing/SpliceBERT/SpliceBERT-main/models/models/SpliceBERT.1024nt \
    --output_prefix results/splicebert_rat \
    --flanking_window 100 \
    --batch_size 4
```

**Expected output files:**
- `results/splicebert_rat_scores.jsonl`
- `results/splicebert_rat_scores.csv`

## Step 4: Evaluate Performance

```bash
# Create results directory if it doesn't exist
mkdir -p results/evaluation

python evaluate_splicebert_rat.py \
    --labels /mnt/userdata4/splicing/ratgetx/processed_data/ratgtex_silver_benchmark_balanced.tsv \
    --scores results/splicebert_rat_scores.csv \
    --output_dir results/evaluation
	

**Expected output files:**
- `results/evaluation/splicebert_rat_evaluation_results.csv`
- `results/evaluation/splicebert_rat_evaluation_summary.json`
- `results/evaluation/roc_curve.png`
- `results/evaluation/precision_recall_curve.png`
- `results/evaluation/score_distribution.png`
- `results/evaluation/merged_data.csv`

## Step 5: Download Results (Optional)

To analyze results locally:

```bash
# From your local machine, download the results
scp -r user@server:/mnt/userdata4/splicing/SpliceBERT/rat_analysis/results /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/ratGTEx/baselines/SpliceBERT/
```

## Expected Performance

Based on SpliceBERT's characteristics:
- **Better than SpliceAI**: Should perform better than SpliceAI's AUROC 0.57 due to multi-species pre-training
- **RNA-centric advantage**: Uses RNA sequences, which may be more conserved across species
- **Coverage**: Should achieve 100% coverage (no annotation dependency like SpliceAI)

## Troubleshooting

### Common Issues:

1. **CUDA memory error**: Reduce batch size to 2 or 1
2. **Model loading error**: Verify model path exists
3. **Tokenization issues**: Check if transformers library is compatible

### Resource Requirements:

- **GPU memory**: ~8-12GB for batch size 4
- **Time**: ~2-4 hours for 28,120 variants
- **Disk space**: ~500MB for results

## Notes

- SpliceBERT requires RNA format (T->U conversion) and space-separated tokens
- The model uses KL-divergence scoring in a ±100nt flanking window
- Lower KL-context scores typically indicate more splice-altering variants (inverted scoring)
- This is a zero-shot cross-species evaluation (no fine-tuning on rat data)
