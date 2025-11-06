# DNABERT-2 Baseline for Rat Data

## Overview
Apply successful DNABERT-2 MLP + Concatenate Embeddings strategy (proven on human SpliceVarDB) to rat data for cross-species evaluation.

**Strategy**: Concatenate `[WT, MT, diff]` embeddings (2304D) + MLP [256, 128]

---

## Scripts

1. **predict_dnabert2_rat.py** - Generate concatenated embeddings
   - Parse variant_id to extract ref/alt (format: `chr:pos_ref/alt`)
   - Generate WT and MT embeddings
   - Compute difference and concatenate

2. **train_mlp_rat.py** - Train MLP classifier
   - Architecture: [256, 128] hidden layers
   - Preprocessing: StandardScaler
   - Output: AUROC, AUPRC, confusion matrix

3. **run_dnabert2_rat_server.sh** - Server execution script

---

## Server Execution

### Setup
```bash
# On server, create project directory
mkdir -p /mnt/userdata4/splicing/DNABert_2_rat
cd /mnt/userdata4/splicing/DNABert_2_rat

# Copy scripts
scp predict_dnabert2_rat.py user@server:/mnt/userdata4/splicing/DNABert_2_rat/
scp train_mlp_rat.py user@server:/mnt/userdata4/splicing/DNABert_2_rat/
scp run_dnabert2_rat_server.sh user@server:/mnt/userdata4/splicing/DNABert_2_rat/
```

### Run
```bash
# Execute pipeline
bash run_dnabert2_rat_server.sh
```

### Expected Output
```
/mnt/userdata4/splicing/DNABert_2_rat/results/
├── dnabert2_rat_concat_embeddings.npz          # Concatenated embeddings
├── dnabert2_rat_mlp_classifier.joblib          # Trained MLP model
├── dnabert2_rat_mlp_classifier_scaler.joblib   # StandardScaler
├── dnabert2_rat_mlp_metrics.json               # Performance metrics
├── dnabert2_rat_mlp_variant_scores.tsv         # Variant predictions
└── dnabert2_rat_mlp_confusion_matrix.png       # Confusion matrix plot
```

---

## Download Results to Local

After server execution, download results:

```bash
cd /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/ratGTEx/baselines/DNABert_2
mkdir -p results

SERVER="user@your_server"
DIR="/mnt/userdata4/splicing/DNABert_2_rat/results"

# Essential files
scp $SERVER:$DIR/dnabert2_rat_concat_embeddings.npz results/
scp $SERVER:$DIR/dnabert2_rat_mlp_metrics.json results/
scp $SERVER:$DIR/dnabert2_rat_mlp_variant_scores.tsv results/
scp $SERVER:$DIR/dnabert2_rat_mlp_confusion_matrix.png results/

# Optional (models can be regenerated)
scp $SERVER:$DIR/dnabert2_rat_mlp_classifier.joblib results/
scp $SERVER:$DIR/dnabert2_rat_mlp_classifier_scaler.joblib results/
```

---

## Data Format

### Input (ratgtex_silver_benchmark_balanced.tsv)
No header, columns:
```
variant_id    ref_seq    alt_seq    label    chrom
```

Example:
```
15:55290237_A/G    GGCCGC...    GGCCGC...    1    15
```

### Output (dnabert2_rat_mlp_variant_scores.tsv)
```
variant_id    true_label    predicted_label    predicted_proba
```

---

## Expected Performance

**Human (SpliceVarDB)**: AUROC 0.9603, AUPRC 0.9561

**Rat (ratGTEx)**: TBD (cross-species evaluation)

---

## Comparison to Other Models

| Model | Strategy | Human AUROC | Rat AUROC | Notes |
|-------|----------|-------------|-----------|-------|
| SpliceAI | Direct prediction | ~0.95 | ~0.58 | Hardcode bypass |
| Evo2 | $\Delta \log p$ | 0.8916 | 0.7365 | Generative model |
| DNABERT-2 | MLP + concat | 0.9603 | **TBD** | This work |

---

## Notes

- Rat data uses different chromosome naming (numeric vs chr prefix)
- Variant parsing adapted for `chr:pos_ref/alt` format
- Same MLP architecture and hyperparameters as human data for fair comparison

