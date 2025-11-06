# Evo2 Embedding Method - Quick Start Commands

## 🚀 TL;DR

Test if Evo2 **embeddings** (not delta_logp) can match DNABERT-2's success (rat AUROC 0.71).

---

## 📋 Prerequisites Check

```bash
# Server files exist
ssh user@server "ls -lh /mnt/userdata4/splicing/Evo2Splicevardb/model_checkpoint/evo2_7b_base.pt"
ssh user@server "ls -lh /mnt/userdata4/splicing/ratgetx/processed_data/ratgtex_silver_benchmark_balanced.tsv"

# Local directory
cd /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/ratGTEx/baselines/Evo2
mkdir -p results
```

---

## 🔄 Complete Pipeline (Copy-Paste)

### Step 1: Upload Scripts to Server

```bash
cd /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/ratGTEx/baselines/Evo2

scp extract_evo2_embeddings.py user@server:/mnt/userdata4/splicing/Evo2_rat/
scp run_evo2_embeddings_rat_server.sh user@server:/mnt/userdata4/splicing/Evo2_rat/
```

### Step 2: Run Embedding Extraction on Server

```bash
ssh user@server

cd /mnt/userdata4/splicing/Evo2_rat
chmod +x run_evo2_embeddings_rat_server.sh

# Run in background (8-16 hours)
nohup bash run_evo2_embeddings_rat_server.sh > evo2_emb_rat.log 2>&1 &

# Get process ID
echo $!

# Monitor (Ctrl+C to exit monitoring)
tail -f evo2_emb_rat.log

# Check progress
watch -n 30 'ls -lh evo2_rat_embeddings.npz 2>/dev/null || echo "Not ready yet"; tail -3 evo2_emb_rat.log'
```

### Step 3: Download Rat Embeddings

```bash
# On local machine
cd /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/ratGTEx/baselines/Evo2

scp user@server:/mnt/userdata4/splicing/Evo2_rat/evo2_rat_embeddings.npz results/

# Verify download
ls -lh results/evo2_rat_embeddings.npz
python3 -c "import numpy as np; d=np.load('results/evo2_rat_embeddings.npz', allow_pickle=True); print(f'Variants: {len(d[\"variant_ids\"]):,}, Embedding dim: {d[\"embeddings\"].shape[1]:,}D')"
```

### Step 4: Train MLP on Human Data

**Option A: If human embeddings already exist**
```bash
# Use existing human embeddings (if available from Evo2Splicevardb)
python3 train_evo2_mlp.py \
    --input_file /path/to/evo2_human_embeddings.npz \
    --model_file results/evo2_mlp_model.pkl \
    --scaler_file results/evo2_mlp_scaler.pkl \
    --metrics_file results/evo2_mlp_human_metrics.json
```

**Option B: Extract human embeddings first** *(if not available)*
```bash
# TODO: Need to prepare human SpliceVarDB in TSV format first
# Then run extract_evo2_embeddings.py on server for human data
```

### Step 5: Evaluate on Rat (Zero-Shot)

```bash
python3 evaluate_evo2_mlp_rat.py \
    --input_file results/evo2_rat_embeddings.npz \
    --model_file results/evo2_mlp_model.pkl \
    --scaler_file results/evo2_mlp_scaler.pkl \
    --output_dir results

# Check results
cat results/evo2_mlp_rat_metrics.json
open results/evo2_mlp_rat_curves.png  # macOS
```

---

## 📊 Expected Results

```json
{
  "auroc": 0.70-0.75,  // Target: match or exceed DNABERT-2 (0.7112)
  "auprc": 0.68-0.73,  // Should be similar to AUROC
  "accuracy": 0.65-0.70,
  "n_variants": 28120
}
```

**Success Criteria**:
- ✅ AUROC ≥ 0.71 → Evo2 embeddings ≥ DNABERT-2
- ✅ AUROC ≥ 0.48 → Better than Evo2 delta_logp (failed baseline)
- ❌ AUROC < 0.55 → Embeddings don't help, architecture issue

---

## ⏱️ Time Estimates

| Step | Time | Note |
|------|------|------|
| 1. Upload | 1 min | Fast |
| 2. Extract embeddings (server) | **8-16 hrs** | 28k variants, batch_size=1 |
| 3. Download | 5-10 min | ~600MB file |
| 4. Train MLP | 5-15 min | On laptop CPU |
| 5. Evaluate | 1-2 min | Fast |
| **Total** | **~8-16 hrs** | Mostly waiting for Step 2 |

---

## 🔍 Quick Checks

### Monitor Server Job

```bash
# Check if running
ssh user@server "ps aux | grep extract_evo2_embeddings.py"

# GPU usage
ssh user@server "nvidia-smi"

# Log tail
ssh user@server "tail -20 /mnt/userdata4/splicing/Evo2_rat/evo2_emb_rat.log"

# Output file growth
ssh user@server "ls -lh /mnt/userdata4/splicing/Evo2_rat/evo2_rat_embeddings.npz"
```

### Verify NPZ Files

```bash
python3 << 'EOF'
import numpy as np

# Check rat embeddings
data = np.load('results/evo2_rat_embeddings.npz', allow_pickle=True)
print("Rat Embeddings:")
print(f"  Files: {data.files}")
print(f"  Variants: {len(data['variant_ids']):,}")
print(f"  Embedding shape: {data['embeddings'].shape}")
print(f"  Labels: 0={np.sum(data['labels']==0):,}, 1={np.sum(data['labels']==1):,}")
EOF
```

---

## 🚨 Troubleshooting

### Embedding extraction stuck/slow
```bash
# Check GPU is being used
ssh user@server "nvidia-smi | grep python"

# If not using GPU, check container CUDA setup
ssh user@server "apptainer exec --nv /mnt/userdata4/splicing/Evo2Splicevardb/vortex.sif nvidia-smi"
```

### MLP training fails
```bash
# Check embedding file
python3 -c "import numpy as np; d=np.load('results/evo2_rat_embeddings.npz', allow_pickle=True); print({k: type(d[k]) for k in d.files})"

# Should show all arrays, not objects (unless variant_ids)
```

### Results look wrong
```bash
# Compare to baselines
cat << 'BASELINE'
Model                   Rat AUROC
---------------------------------
Evo2 delta_logp:        0.4848    ← Failed baseline
DNABERT-2 MLP:          0.7112    ← Target to match/exceed
Evo2 MLP (this):        0.7X?     ← Your result
BASELINE
```

---

## 📝 Next Actions After Success

1. **Update draft.tex** with Evo2 embedding results
2. **Compare**: Evo2 MLP vs. DNABERT-2 MLP vs. Evo2 delta_logp
3. **Conclude**: Embeddings vs. direct probability scoring for foundation models

---

## 🎯 Scientific Question

**Does Evo2's 1000+ species pre-training beat DNABERT-2's 135 species when both use embeddings+MLP?**

- If YES → Pre-training scale matters
- If NO → Embedding architecture/quality matters more than scale
- Either way → Embeddings > probability-based scoring for variants

---

## 📧 Contact

If extraction fails or results are unexpected, check:
1. `evo2_emb_rat.log` on server for errors
2. StripedHyena model API for hidden state extraction
3. README_EMBEDDING_METHOD.md for detailed troubleshooting

