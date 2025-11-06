# Evo2 Zero-Shot Evaluation on ChickenGTEx

## Overview
Evo2 (7B parameters, 8K context) zero-shot splicing variant prediction on ChickenGTEx balanced benchmark.

## Method
- **Model**: Evo2-7B-8K (generative DNA foundation model)
- **Scoring**: Δlog p = log p(alt) - log p(ref)
- **Higher Δlog p** → more splice-altering

## Server Execution

### 1. Upload files to server
```bash
# Create directory
ssh yinuos@mlerp-monash-node04 "mkdir -p /mnt/userdata4/splicing/Evo2_chicken"

# Upload scripts
scp predict_evo2_chicken.py yinuos@mlerp-monash-node04:/mnt/userdata4/splicing/Evo2_chicken/
scp run_evo2_chicken_server.sh yinuos@mlerp-monash-node04:/mnt/userdata4/splicing/Evo2_chicken/
```

### 2. Run on server
```bash
ssh yinuos@mlerp-monash-node04
cd /mnt/userdata4/splicing/Evo2_chicken
bash run_evo2_chicken_server.sh
```

### 3. Download results
```bash
scp yinuos@mlerp-monash-node04:/mnt/userdata4/splicing/Evo2_chicken/evo2_chicken_predictions.tsv ./results/
```

## Output
- `results/evo2_chicken_predictions.tsv`: Variant-level predictions
  - Columns: variant_id, logp_ref, logp_alt, delta_logp, label

## Expected Performance
- Coverage: 100% (all variants scored)
- AUROC: TBD (compare with rat ~0.48, pig TBD)

---

## Alternative Method: Embedding + MLP (Within-Model Comparison)

### Overview
Extract Evo2 embeddings and train MLP classifier (for strategy comparison, not zero-shot).

### Steps

#### 1. Extract embeddings
```bash
scp extract_evo2_embeddings.py yinuos@mlerp-monash-node04:/mnt/userdata4/splicing/Evo2_chicken/
scp run_evo2_embeddings_chicken_server.sh yinuos@mlerp-monash-node04:/mnt/userdata4/splicing/Evo2_chicken/

ssh yinuos@mlerp-monash-node04
cd /mnt/userdata4/splicing/Evo2_chicken
bash run_evo2_embeddings_chicken_server.sh
```

#### 2. Train MLP
```bash
scp train_evo2_mlp.py yinuos@mlerp-monash-node04:/mnt/userdata4/splicing/Evo2_chicken/
scp run_evo2_mlp_chicken_server.sh yinuos@mlerp-monash-node04:/mnt/userdata4/splicing/Evo2_chicken/

ssh yinuos@mlerp-monash-node04
cd /mnt/userdata4/splicing/Evo2_chicken
bash run_evo2_mlp_chicken_server.sh
```

#### 3. Download results
```bash
scp yinuos@mlerp-monash-node04:/mnt/userdata4/splicing/Evo2_chicken/results/evo2_mlp_chicken_predictions.tsv ./results/
scp yinuos@mlerp-monash-node04:/mnt/userdata4/splicing/Evo2_chicken/results/evo2_mlp_chicken_metrics_noplot.json ./results/
```

### Output Files
- `evo2_chicken_embeddings.npz`: Raw embeddings
- `results/evo2_mlp_chicken_model.pkl`: Trained MLP model
- `results/evo2_mlp_chicken_predictions.tsv`: Predictions
- `results/evo2_mlp_chicken_metrics_noplot.json`: Performance metrics

