# Evo2 Embedding + MLP Method

## Overview

This directory implements an **alternative approach** to Evo2 variant effect prediction:
- **Method 1 (Failed)**: Direct `Δlog p` scoring → AUROC 0.4848 on rat (quantization issues)
- **Method 2 (This)**: Extract embeddings + train MLP → Expected to match/exceed DNABERT-2 (AUROC 0.71)

## Motivation

The direct `Δlog p` scoring method failed catastrophically on rat data (AUROC 0.48, worse than random) because:
1. **95.58% of variants have Δlog p = 0** (model cannot distinguish them)
2. **Coarse quantization** (all logp values are multiples of 64)
3. **Architectural insensitivity** to single-nucleotide changes

However, Evo2's **hidden state representations** (embeddings) may preserve fine-grained information lost in the quantized probability outputs. This approach mirrors DNABERT-2's success:
- DNABERT-2 with concatenate+MLP: human AUROC 0.96, rat 0.71
- Evo2 has **1000+ species pre-training** vs. DNABERT-2's 135 species
- Hypothesis: Evo2 embeddings could **match or exceed DNABERT-2** performance

## Pipeline Overview

```
Human Data (SpliceVarDB)          Rat Data (ratGTEx)
        ↓                               ↓
[1] Extract Embeddings        [3] Extract Embeddings
        ↓                               ↓
[2] Train MLP                          ↓
        ↓                               ↓
        └───────── [4] Evaluate ←──────┘
                (Zero-Shot Cross-Species)
```

## File Structure

```
Evo2/
├── extract_evo2_embeddings.py          # Extract hidden states from Evo2
├── train_evo2_mlp.py                   # Train MLP on embeddings
├── evaluate_evo2_mlp_rat.py            # Evaluate on rat (zero-shot)
├── run_evo2_embeddings_rat_server.sh   # Server script for rat embedding extraction
├── README_EMBEDDING_METHOD.md          # This file
└── results/                            # Output directory
```

## Prerequisites

### Server Environment
- Evo2 model checkpoint: `/mnt/userdata4/splicing/Evo2Splicevardb/model_checkpoint/evo2_7b_base.pt`
- Container: `/mnt/userdata4/splicing/Evo2Splicevardb/vortex.sif`
- Config: `/mnt/userdata4/splicing/Evo2Splicevardb/vortex-main/configs/evo2-7b-8k.yml`

### Local Environment
```bash
pip install numpy pandas scikit-learn joblib matplotlib seaborn tqdm
```

## Complete Workflow

### Step 1: Extract Human Embeddings (Server)

**Note**: For human data, you need to first prepare SpliceVarDB in the same format as rat data.

```bash
# On local machine
cd /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/ratGTEx/baselines/Evo2

# Upload embedding extraction script
scp extract_evo2_embeddings.py user@server:/mnt/userdata4/splicing/Evo2_rat/

# SSH to server
ssh user@server

# Run embedding extraction for human data (if available)
cd /mnt/userdata4/splicing/Evo2_rat

python extract_evo2_embeddings.py \
    --config_path /mnt/userdata4/splicing/Evo2Splicevardb/vortex-main/configs/evo2-7b-8k.yml \
    --checkpoint_path /mnt/userdata4/splicing/Evo2Splicevardb/model_checkpoint/evo2_7b_base.pt \
    --input_file /path/to/human/splicevardb_8192bp.tsv \
    --output_file evo2_human_embeddings.npz \
    --batch_size 1

# Download embeddings
# On local machine
scp user@server:/mnt/userdata4/splicing/Evo2_rat/evo2_human_embeddings.npz \
    /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/ratGTEx/baselines/Evo2/results/
```

**Expected Output**:
- File: `evo2_human_embeddings.npz`
- Size: ~500MB-1GB (depends on embedding dimension)
- Contains: variant_ids, embeddings (concatenated: ref+alt+diff), labels

**Estimated Time**: 6-12 hours for ~24k variants

### Step 2: Train MLP on Human Embeddings (Local)

```bash
cd /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/ratGTEx/baselines/Evo2

python train_evo2_mlp.py \
    --input_file results/evo2_human_embeddings.npz \
    --model_file results/evo2_mlp_model.pkl \
    --scaler_file results/evo2_mlp_scaler.pkl \
    --metrics_file results/evo2_mlp_human_metrics.json
```

**Expected Output**:
- Model: `evo2_mlp_model.pkl` (~1MB)
- Scaler: `evo2_mlp_scaler.pkl` (~100KB)
- Metrics: `evo2_mlp_human_metrics.json`
  - Expected AUROC: **0.85-0.95** (should match or exceed Evo2 delta_logp's 0.71)

**Training Time**: 5-15 minutes on laptop

### Step 3: Extract Rat Embeddings (Server)

```bash
# On local machine
cd /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/ratGTEx/baselines/Evo2

# Upload scripts
scp extract_evo2_embeddings.py user@server:/mnt/userdata4/splicing/Evo2_rat/
scp run_evo2_embeddings_rat_server.sh user@server:/mnt/userdata4/splicing/Evo2_rat/

# SSH to server
ssh user@server

# Navigate to working directory
cd /mnt/userdata4/splicing/Evo2_rat

# Make script executable
chmod +x run_evo2_embeddings_rat_server.sh

# Run embedding extraction
bash run_evo2_embeddings_rat_server.sh

# OR run in background
nohup bash run_evo2_embeddings_rat_server.sh > evo2_embeddings_rat.log 2>&1 &

# Monitor progress
tail -f evo2_embeddings_rat.log
watch -n 10 'ls -lh evo2_rat_embeddings.npz 2>/dev/null || echo "Not created yet"'
```

**Expected Output**:
- File: `evo2_rat_embeddings.npz`
- Size: ~600MB-1.2GB (28,120 variants)
- Contains: variant_ids, embeddings (concatenated), labels

**Estimated Time**: 8-16 hours (28k variants, batch_size=1)

**Download Results**:
```bash
# On local machine
scp user@server:/mnt/userdata4/splicing/Evo2_rat/evo2_rat_embeddings.npz \
    /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/ratGTEx/baselines/Evo2/results/
```

### Step 4: Evaluate on Rat Data (Local - Zero-Shot)

```bash
cd /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/ratGTEx/baselines/Evo2

python evaluate_evo2_mlp_rat.py \
    --input_file results/evo2_rat_embeddings.npz \
    --model_file results/evo2_mlp_model.pkl \
    --scaler_file results/evo2_mlp_scaler.pkl \
    --output_dir results
```

**Expected Output**:
- Metrics: `results/evo2_mlp_rat_metrics.json`
  - **Target**: AUROC **0.70-0.75** (match or exceed DNABERT-2's 0.71)
  - **Baseline**: Evo2 delta_logp = 0.48 (failed)
- Predictions: `results/evo2_mlp_rat_predictions.tsv`
- Plots:
  - `results/evo2_mlp_rat_curves.png` (ROC + PR curves)
  - `results/evo2_mlp_rat_confusion_matrix.png`

**Evaluation Time**: 1-2 minutes

## Expected Performance Comparison

| Method | Human AUROC | Rat AUROC | Cross-Species Drop |
|--------|-------------|-----------|-------------------|
| **Evo2 Δlog p** | 0.7129 (AUPRC) | **0.4848** | -31% (failed) |
| **DNABERT-2 MLP** | 0.9603 | 0.7112 | -25.9% |
| **Evo2 MLP (this)** | 0.85-0.95 (target) | **0.70-0.75 (target)** | ~20-25% (expected) |

**Hypothesis**: Evo2 embeddings should perform **at least as well as DNABERT-2** because:
1. Evo2 pre-trained on **1000+ species** vs. DNABERT-2's 135
2. Embeddings bypass the quantization bottleneck affecting delta_logp
3. 8,192bp context vs. DNABERT-2's typically shorter windows

## Troubleshooting

### If Embedding Extraction Fails

**Check model forward signature**:
```python
# In extract_evo2_embeddings.py, around line 50-60
# Try different extraction strategies:

# Strategy 1: output_hidden_states parameter
outputs = model(input_tensor, output_hidden_states=True)

# Strategy 2: Access model internals
# May need to inspect vortex.model.model.StripedHyena source code
```

**OOM (Out of Memory)**:
```bash
# Reduce batch size to 1 (already default)
# Check GPU memory: nvidia-smi
# Kill other processes if needed
```

### If Embeddings Are Wrong Dimension

```python
# Check in extract_evo2_embeddings.py
# Line ~73: print(f"Hidden states shape: {hidden_states.shape}")
# Expected: (batch_size, seq_len, hidden_dim)
# After pooling: (batch_size, hidden_dim)
```

### If MLP Training Fails

```bash
# Check embedding file integrity
python -c "import numpy as np; d = np.load('results/evo2_human_embeddings.npz'); print(d.files); print({k: d[k].shape for k in d.files})"

# Expected output:
# ['variant_ids', 'embeddings', 'labels', 'embedding_dim']
# {'variant_ids': (N,), 'embeddings': (N, D*3), 'labels': (N,)}
```

## Next Steps After Results

1. **Compare to DNABERT-2**:
   - If Evo2 MLP ≥ DNABERT-2 (0.71): Evo2 embeddings are superior (larger pre-training pays off)
   - If Evo2 MLP < DNABERT-2: Embedding quality or architecture differences matter

2. **Update Draft.tex**:
   - Add Evo2 embedding results to Section 2.3.2
   - Update summary paragraph with embedding vs. probability-based scoring comparison
   - Add to Key Insights: "Generative foundation models' value lies in representations, not direct probability outputs"

3. **Consider NT and SpliceBERT**:
   - Apply same embedding+MLP strategy to NT (if embeddings accessible)
   - May rescue their near-random zero-shot performance

## Scientific Implications

This experiment tests a fundamental hypothesis:
> **For variant effect prediction, are foundation models' intermediate representations (embeddings) more valuable than their task-specific outputs (probabilities, distances)?**

**If successful**, this demonstrates:
1. Generative and discriminative foundation models **converge** via embeddings
2. Pre-training scale matters for embeddings, not just model capacity
3. Task-specific fine-tuning can be replaced by **frozen embeddings + lightweight classifiers**

This would fundamentally change how we use foundation models for genomic variant prediction.

## Citation

If you use this method, cite:
- Evo2: Nguyen et al., 2024 (generative genomic foundation model)
- DNABERT-2: Zhou et al., 2023 (embedding extraction + downstream classifier paradigm)
- Your review paper: Sun et al., 2025 (systematic comparison of foundation model utilization strategies)

