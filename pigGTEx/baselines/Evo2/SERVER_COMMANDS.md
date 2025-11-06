# Evo2 PigGTEx Server Commands

## Method 1: Zero-Shot (Direct Scoring)

### Uploadi ge
```bash
scp predict_evo2_pig.py YOUR_SERVER:/mnt/userdata4/splicing/Evo2_pig/
scp run_evo2_pig_server.sh YOUR_SERVER:/mnt/userdata4/splicing/Evo2_pig/
```

### Run
```bash
ssh YOUR_SERVER
cd /mnt/userdata4/splicing/Evo2_pig
bash run_evo2_pig_server.sh
```

### Download
```bash
scp YOUR_SERVER:/mnt/userdata4/splicing/Evo2_pig/evo2_pig_predictions.tsv \
  /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/pigGTEx/baselines/Evo2/results/
```

---

## Method 2: Embedding + MLP (Within-Model Strategy Comparison)

### Step 1: Extract Embeddings

#### Upload
```bash
scp extract_evo2_embeddings.py YOUR_SERVER:/mnt/userdata4/splicing/Evo2_pig/
scp run_evo2_embeddings_pig_server.sh YOUR_SERVER:/mnt/userdata4/splicing/Evo2_pig/
```

#### Run
```bash
ssh YOUR_SERVER
cd /mnt/userdata4/splicing/Evo2_pig
bash run_evo2_embeddings_pig_server.sh
```

### Step 2: Train MLP

#### Upload
```bash
scp train_evo2_mlp.py YOUR_SERVER:/mnt/userdata4/splicing/Evo2_pig/
scp run_evo2_mlp_pig_server.sh YOUR_SERVER:/mnt/userdata4/splicing/Evo2_pig/
```

#### Run
```bash
ssh YOUR_SERVER
cd /mnt/userdata4/splicing/Evo2_pig
bash run_evo2_mlp_pig_server.sh
```

### Download Results
```bash
# Predictions
scp YOUR_SERVER:/mnt/userdata4/splicing/Evo2_pig/results/evo2_mlp_pig_predictions.tsv \
  /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/pigGTEx/baselines/Evo2/results/

# Metrics
scp YOUR_SERVER:/mnt/userdata4/splicing/Evo2_pig/results/evo2_mlp_pig_metrics_noplot.json \
  /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/pigGTEx/baselines/Evo2/results/
```
