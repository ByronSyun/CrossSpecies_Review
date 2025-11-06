# DNABERT-2 Evaluation on PigGTEx

## Overview
DNABERT-2 (117M parameters) evaluation on PigGTEx with two embedding-based methods.

## Methods

### Method 1: Logistic Regression (Differential Embeddings)
- **Embeddings**: ALT - REF (difference)
- **Classifier**: Logistic Regression
- **Purpose**: Simpler baseline, close to zero-shot

**Server execution:**
```bash
# Upload
scp predict_dnabert2_logistic_pig.py YOUR_SERVER:/mnt/userdata4/splicing/DNABert_2_pig/
scp train_classifier_pig.py YOUR_SERVER:/mnt/userdata4/splicing/DNABert_2_pig/
scp run_dnabert2_logistic_pig_server.sh YOUR_SERVER:/mnt/userdata4/splicing/DNABert_2_pig/

# Run
ssh YOUR_SERVER
cd /mnt/userdata4/splicing/DNABert_2_pig
bash run_dnabert2_logistic_pig_server.sh

# Download
scp YOUR_SERVER:/mnt/userdata4/splicing/DNABert_2_pig/results/dnabert2_pig_logistic_scores.tsv ./results/
```

**Output:**
- `dnabert2_pig_logistic_scores.tsv`: Predictions (prob_splice_altering, label)

---

### Method 2: MLP (Concatenated Embeddings)
- **Embeddings**: [REF, ALT, ALT-REF] (concatenated)
- **Classifier**: Multi-layer Perceptron (768 → 128 → 1)
- **Purpose**: Within-model strategy comparison (not zero-shot)

**Server execution:**
```bash
# Upload
scp predict_dnabert2_pig.py YOUR_SERVER:/mnt/userdata4/splicing/DNABert_2_pig/
scp train_mlp_pig.py YOUR_SERVER:/mnt/userdata4/splicing/DNABert_2_pig/
scp run_dnabert2_pig_server.sh YOUR_SERVER:/mnt/userdata4/splicing/DNABert_2_pig/

# Run
ssh YOUR_SERVER
cd /mnt/userdata4/splicing/DNABert_2_pig
bash run_dnabert2_pig_server.sh

# Download
scp YOUR_SERVER:/mnt/userdata4/splicing/DNABert_2_pig/results/dnabert2_pig_mlp_variant_scores.tsv ./results/
scp YOUR_SERVER:/mnt/userdata4/splicing/DNABert_2_pig/results/dnabert2_pig_mlp_metrics.json ./results/
```

**Output:**
- `dnabert2_pig_mlp_variant_scores.tsv`: Predictions
- `dnabert2_pig_mlp_metrics.json`: Performance metrics

## Model Details
- **Model**: zhihan1996/DNABERT-2-117M
- **Embedding dim**: 768
- **Tokenization**: 6-mer k-mer tokenization

## Expected Performance
- Coverage: 100%
- AUROC: TBD (compare with rat)

