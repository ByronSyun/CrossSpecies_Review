# DNABERT-2 Logistic Regression Baseline for Rat Data

## Overview

This directory contains scripts to run the **DNABERT-2 Logistic Regression** baseline on rat (RatGTEx) data. This method uses **differential embeddings (ALT - REF)** with a simple Logistic Regression classifier, representing a zero-shot approach with minimal supervision.

## Files

### Logistic Regression Method (Difference Embeddings)
- `predict_dnabert2_logistic_rat.py` - Extract differential embeddings (ALT - REF)
- `train_classifier_rat.py` - Train Logistic Regression on differential embeddings
- `run_dnabert2_logistic_rat_server.sh` - Complete pipeline for server

### MLP Method (Concatenate Embeddings) - Separate Implementation
- `predict_dnabert2_rat.py` - Extract concatenated embeddings (REF + ALT + diff)
- `train_mlp_rat.py` - Train MLP on concatenated embeddings
- `run_dnabert2_rat_server.sh` - Complete MLP pipeline for server

## Running on Server

### 1. Upload scripts to server

```bash
# On local machine
scp predict_dnabert2_logistic_rat.py train_classifier_rat.py run_dnabert2_logistic_rat_server.sh \
    username@server:/mnt/userdata4/splicing/DNABert_2_rat/
```

### 2. Run the pipeline on server

```bash
# On server
cd /mnt/userdata4/splicing/DNABert_2_rat/
bash run_dnabert2_logistic_rat_server.sh
```

## Expected Outputs

After successful execution, you will find in `/mnt/userdata4/splicing/DNABert_2_rat/results/`:

- `dnabert2_rat_logistic_embeddings.npz` - Differential embeddings (ALT - REF)
- `dnabert2_rat_logistic_classifier.joblib` - Trained Logistic Regression model
- `dnabert2_rat_logistic_confusion_matrix.png` - Confusion matrix visualization
- **`dnabert2_rat_logistic_scores.tsv`** - Per-variant predictions (for benchmarking)

The final TSV file has the following format:
```
variant_id	prob_splice_altering	label
15:55290237_A/G	0.523456	1
7:3089418_A/T	0.487123	1
...
```

## Method Comparison

| Method | Embedding Strategy | Classifier | Use Case |
|--------|-------------------|------------|----------|
| **Logistic (this)** | Differential (ALT - REF) | Logistic Regression | Zero-shot baseline |
| **MLP** | Concatenate (REF + ALT + diff) | Multi-layer Perceptron | Higher capacity |

## Key Differences from Human Pipeline

1. **Input Format**: Rat data has no header and uses format `variant_id, ref_seq, alt_seq, label, tissue_id`
2. **Sequences**: Rat data already has separate `ref_seq` and `alt_seq` columns (no need to generate mutants)
3. **Labels**: Already in binary format (0/1) instead of strings
4. **Genome**: Rnor_6.0 (rn6) instead of GRCh38 (hg38)

## Integration with Summary Scripts

After downloading the results, update the path in `/Users/byronsun/Desktop/AS_复现模型/BIB_review/code/ratGTEx/baselines/rat_summary_results.py`:

```python
dnabert2_log_file = base_dir / 'code/ratGTEx/baselines/DNABert_2/results/dnabert2_rat_logistic_scores.tsv'
```

Then run:
```bash
python rat_summary_results.py
python rat_evaluate_results.py
```

## Notes

- Batch size is set to 16 (can be adjusted based on GPU memory)
- Uses balanced class weights in Logistic Regression
- Train/test split is 80/20 with stratification
- All file naming uses `_logistic_` to distinguish from MLP method

