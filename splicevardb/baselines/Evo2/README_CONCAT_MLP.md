# Evo2 Concatenate Embeddings + MLP Pipeline

## Overview

This pipeline implements the **embedding-based approach** for Evo2, similar to DNABERT-2's successful concatenate + MLP method. Instead of using Evo2's direct `delta_logp` scoring (which suffers from quantization issues), we extract hidden state embeddings and train a downstream MLP classifier.

## Workflow

1. **Extract Embeddings** (`extract_evo2_embeddings.py`)
   - Extract hidden state embeddings from Evo2's penultimate layer for both WT and MT sequences
   - Concatenate: `[WT_embedding, MT_embedding, diff_embedding]`
   - Save as compressed NPZ file

2. **Train MLP Classifier** (`train_evo2_mlp.py`)
   - Train a 2-layer MLP [256, 128] on concatenated embeddings
   - 80/20 train/test split with stratification
   - StandardScaler normalization
   - Early stopping to prevent overfitting
   - Save model, scaler, metrics, and test predictions

## Files

- `extract_evo2_embeddings.py` - Embedding extraction script
- `train_evo2_mlp.py` - MLP training script
- `run_evo2_concat_mlp.sh` - Complete pipeline for server (Apptainer)

## Usage

### On Server (Apptainer/Singularity)

```bash
cd /mnt/userdata4/splicing/Evo2Splicevardb
bash run_evo2_concat_mlp.sh
```

### Local Testing (if you have Evo2 model locally)

```bash
# Step 1: Extract embeddings
python extract_evo2_embeddings.py \
    --config_path /path/to/evo2-7b-8k.yml \
    --checkpoint_path /path/to/evo2_7b_base.pt \
    --input_file /path/to/input.tsv \
    --output_file evo2_concat_embeddings.npz \
    --batch_size 2

# Step 2: Train MLP
python train_evo2_mlp.py \
    --input_file evo2_concat_embeddings.npz \
    --model_file evo2_mlp_model.pkl \
    --scaler_file evo2_mlp_scaler.pkl \
    --metrics_file evo2_mlp_train_test_metrics.json

# Step 3: Generate predictions for all variants (using Python script)
# See run_evo2_concat_mlp.sh for the inline Python code
```

## Expected Output

After successful execution:
- `evo2_concat_embeddings.npz` - Concatenated embeddings (WT + MT + diff) for all variants
- `evo2_mlp_model.pkl` - Trained MLP classifier
- `evo2_mlp_scaler.pkl` - StandardScaler for feature normalization
- `evo2_mlp_train_test_metrics.json` - Train/test split performance metrics
- **`evo2_mlp_variant_scores.tsv`** - **Full predictions TSV for all variants** (can be used for downstream evaluation)
- `evo2_mlp_full_metrics.json` - Performance metrics on all variants

## Key Differences from Direct Scoring

| Method | Pros | Cons |
|--------|------|------|
| **Direct `delta_logp`** | No training data needed; immediate deployment | Quantization issues; poor discrimination (AUROC ~0.48 on rat) |
| **Embedding + MLP** | Better performance; bypasses quantization | Requires labeled training data; more complex pipeline |

## Notes

- **Consistency with rat implementation**: This pipeline uses the same method as successfully deployed on ratGTEx data (AUROC 0.7761)
- **Data format handling**: SpliceVarDB has a single `sequence` column (WT), so we manually construct mutant sequences. RatGTEx has separate `ref_seq` and `alt_seq` columns
- **MLP hyperparameters**: Aligned with rat implementation (`batch_size=64`, `max_iter=1000`, `n_iter_no_change=20`)
- Batch size for embedding extraction is kept small (default=2) due to Evo2's large model size and long sequence support (8192 bp)
- Embedding extraction is the bottleneck - expect ~1-2 hours for 23K variants
- MLP training is fast (~5-10 minutes)
- This approach should work better for cross-species transfer than direct scoring

