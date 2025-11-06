# Evo2 Zero-Shot Evaluation on PigGTEx

## Overview
Evo2 (7B parameters, 8K context) zero-shot splicing variant prediction on PigGTEx balanced benchmark.

## Method
- **Model**: Evo2-7B-8K (generative DNA foundation model)
- **Scoring**: Δlog p = log p(alt) - log p(ref)
- **Higher Δlog p** → more splice-altering

## Server Execution

### 1. Upload files to server
```bash
# Create directory
ssh YOUR_SERVER "mkdir -p /mnt/userdata4/splicing/Evo2_pig"

# Upload scripts
scp predict_evo2_pig.py YOUR_SERVER:/mnt/userdata4/splicing/Evo2_pig/
scp run_evo2_pig_server.sh YOUR_SERVER:/mnt/userdata4/splicing/Evo2_pig/
```

### 2. Run on server
```bash
ssh YOUR_SERVER
cd /mnt/userdata4/splicing/Evo2_pig
bash run_evo2_pig_server.sh
```

### 3. Download results
```bash
scp YOUR_SERVER:/mnt/userdata4/splicing/Evo2_pig/evo2_pig_predictions.tsv ./results/
```

## Output
- `results/evo2_pig_predictions.tsv`: Variant-level predictions
  - Columns: variant_id, logp_ref, logp_alt, delta_logp, label

## Expected Performance
- Coverage: 100% (all variants scored)
- AUROC: TBD (compare with rat ~0.48)

---

## Alternative Method: Embedding + MLP (Within-Model Comparison)

### Overview
Extract Evo2 embeddings and train MLP classifier (for strategy comparison, not zero-shot).

### Steps

#### 1. Extract embeddings
```bash
scp extract_evo2_embeddings.py YOUR_SERVER:/mnt/userdata4/splicing/Evo2_pig/
scp run_evo2_embeddings_pig_server.sh YOUR_SERVER:/mnt/userdata4/splicing/Evo2_pig/

ssh YOUR_SERVER
cd /mnt/userdata4/splicing/Evo2_pig
bash run_evo2_embeddings_pig_server.sh
```

#### 2. Train MLP
```bash
scp train_evo2_mlp.py YOUR_SERVER:/mnt/userdata4/splicing/Evo2_pig/
scp run_evo2_mlp_pig_server.sh YOUR_SERVER:/mnt/userdata4/splicing/Evo2_pig/

ssh YOUR_SERVER
cd /mnt/userdata4/splicing/Evo2_pig
bash run_evo2_mlp_pig_server.sh
```

#### 3. Download results
```bash
scp YOUR_SERVER:/mnt/userdata4/splicing/Evo2_pig/results/evo2_mlp_pig_predictions.tsv ./results/
scp YOUR_SERVER:/mnt/userdata4/splicing/Evo2_pig/results/evo2_mlp_pig_metrics_noplot.json ./results/
```

### Output Files
- `evo2_pig_embeddings.npz`: Raw embeddings
- `results/evo2_mlp_pig_model.pkl`: Trained MLP model
- `results/evo2_mlp_pig_predictions.tsv`: Predictions
- `results/evo2_mlp_pig_metrics_noplot.json`: Performance metrics

