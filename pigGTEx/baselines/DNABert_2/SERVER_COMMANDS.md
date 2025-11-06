# DNABERT-2 PigGTEx Server Commands

## Method 1: Logistic Regression (Differential Embeddings)

### Upload
```bash
scp predict_dnabert2_logistic_pig.py YOUR_SERVER:/mnt/userdata4/splicing/DNABert_2_pig/
scp train_classifier_pig.py YOUR_SERVER:/mnt/userdata4/splicing/DNABert_2_pig/
scp run_dnabert2_logistic_pig_server.sh YOUR_SERVER:/mnt/userdata4/splicing/DNABert_2_pig/
```

### Run
```bash
ssh YOUR_SERVER
cd /mnt/userdata4/splicing/DNABert_2_pig
bash run_dnabert2_logistic_pig_server.sh
```

### Download
```bash
scp YOUR_SERVER:/mnt/userdata4/splicing/DNABert_2_pig/results/dnabert2_pig_logistic_scores.tsv \
  /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/pigGTEx/baselines/DNABert_2/results/
```

---

## Method 2: MLP (Concatenated Embeddings)

### Upload
```bash
scp predict_dnabert2_pig.py YOUR_SERVER:/mnt/userdata4/splicing/DNABert_2_pig/
scp train_mlp_pig.py YOUR_SERVER:/mnt/userdata4/splicing/DNABert_2_pig/
scp run_dnabert2_pig_server.sh YOUR_SERVER:/mnt/userdata4/splicing/DNABert_2_pig/
```

### Run
```bash
ssh YOUR_SERVER
cd /mnt/userdata4/splicing/DNABert_2_pig
bash run_dnabert2_pig_server.sh
```

### Download
```bash
# Predictions
scp YOUR_SERVER:/mnt/userdata4/splicing/DNABert_2_pig/results/dnabert2_pig_mlp_variant_scores.tsv \
  /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/pigGTEx/baselines/DNABert_2/results/

# Metrics
scp YOUR_SERVER:/mnt/userdata4/splicing/DNABert_2_pig/results/dnabert2_pig_mlp_metrics.json \
  /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/pigGTEx/baselines/DNABert_2/results/
```

