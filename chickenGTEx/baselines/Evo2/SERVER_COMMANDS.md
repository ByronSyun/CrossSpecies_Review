# ChickenGTEx Evo2 - Complete Server Workflow

## Part 1: Zero-Shot Prediction (Δlog p)

### Step 1: Upload scripts
```bash
cd /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/chickenGTEx/baselines/Evo2

scp predict_evo2_chicken.py run_evo2_chicken_server.sh \
    yinuos@mlerp-monash-node04:/mnt/userdata4/splicing/Evo2_chicken/
```

### Step 2: Run prediction on server
```bash
ssh yinuos@mlerp-monash-node04
cd /mnt/userdata4/splicing/Evo2_chicken
bash run_evo2_chicken_server.sh
```

### Step 3: Download results
```bash
scp yinuos@mlerp-monash-node04:/mnt/userdata4/splicing/Evo2_chicken/evo2_chicken_predictions.tsv \
    /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/chickenGTEx/baselines/Evo2/results/
```

### Step 4: Evaluate locally
```bash
cd /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/chickenGTEx/baselines/Evo2

python3 - << 'EOF'
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score, accuracy_score

df = pd.read_csv('results/evo2_chicken_predictions.tsv', sep='\t')
y_true = df['label']
y_score = df['delta_logp']
y_pred = (y_score > 0).astype(int)

print(f"ChickenGTEx Evo2 Zero-Shot Results:")
print(f"  AUROC: {roc_auc_score(y_true, y_score):.4f}")
print(f"  AUPRC: {average_precision_score(y_true, y_score):.4f}")
print(f"  Accuracy: {accuracy_score(y_true, y_pred):.4f}")
print(f"  Total variants: {len(df):,}")
EOF
```

---

## Part 2: Embedding + MLP

### Step 1: Upload embedding extraction scripts
```bash
cd /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/chickenGTEx/baselines/Evo2

scp extract_evo2_embeddings.py run_evo2_embeddings_chicken_server.sh \
    yinuos@mlerp-monash-node04:/mnt/userdata4/splicing/Evo2_chicken/
```

### Step 2: Extract embeddings on server
```bash
ssh yinuos@mlerp-monash-node04
cd /mnt/userdata4/splicing/Evo2_chicken
bash run_evo2_embeddings_chicken_server.sh
```

### Step 3: Upload MLP training scripts
```bash
cd /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/chickenGTEx/baselines/Evo2

scp train_evo2_mlp.py run_evo2_mlp_chicken_server.sh \
    yinuos@mlerp-monash-node04:/mnt/userdata4/splicing/Evo2_chicken/
```

### Step 4: Train MLP on server
```bash
ssh yinuos@mlerp-monash-node04
cd /mnt/userdata4/splicing/Evo2_chicken
bash run_evo2_mlp_chicken_server.sh
```

### Step 5: Download results
```bash
scp yinuos@mlerp-monash-node04:/mnt/userdata4/splicing/Evo2_chicken/results/evo2_mlp_chicken_predictions.tsv \
    /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/chickenGTEx/baselines/Evo2/results/

scp yinuos@mlerp-monash-node04:/mnt/userdata4/splicing/Evo2_chicken/results/evo2_mlp_chicken_metrics_noplot.json \
    /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/chickenGTEx/baselines/Evo2/results/
```

### Step 6: View results
```bash
cd /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/chickenGTEx/baselines/Evo2

cat results/evo2_mlp_chicken_metrics_noplot.json
```

---

## Quick Reference

### All upload commands (run from local)
```bash
cd /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/chickenGTEx/baselines/Evo2

# Create directory
ssh yinuos@mlerp-monash-node04 "mkdir -p /mnt/userdata4/splicing/Evo2_chicken/results"

# Upload all scripts
scp predict_evo2_chicken.py extract_evo2_embeddings.py train_evo2_mlp.py \
    run_evo2_chicken_server.sh run_evo2_embeddings_chicken_server.sh run_evo2_mlp_chicken_server.sh \
    yinuos@mlerp-monash-node04:/mnt/userdata4/splicing/Evo2_chicken/
```

### All server commands (run on server)
```bash
ssh yinuos@mlerp-monash-node04
cd /mnt/userdata4/splicing/Evo2_chicken

# 1. Zero-shot prediction
bash run_evo2_chicken_server.sh

# 2. Extract embeddings
bash run_evo2_embeddings_chicken_server.sh

# 3. Train MLP
bash run_evo2_mlp_chicken_server.sh
```

### All download commands (run from local)
```bash
cd /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/chickenGTEx/baselines/Evo2

scp yinuos@mlerp-monash-node04:/mnt/userdata4/splicing/Evo2_chicken/evo2_chicken_predictions.tsv ./results/
scp yinuos@mlerp-monash-node04:/mnt/userdata4/splicing/Evo2_chicken/results/evo2_mlp_chicken_*.tsv ./results/
scp yinuos@mlerp-monash-node04:/mnt/userdata4/splicing/Evo2_chicken/results/evo2_mlp_chicken_*.json ./results/
```

---

## Expected Runtime
- **Zero-shot prediction**: ~8-12 hours (25,000 variants, batch_size=1)
- **Embedding extraction**: ~8-12 hours
- **MLP training**: ~5-10 minutes

## Expected Results
Based on rat/pig performance:
- **Zero-shot**: AUROC ~0.45-0.50 (baseline)
- **Embedding+MLP**: AUROC ~0.70-0.75 (significant improvement)

