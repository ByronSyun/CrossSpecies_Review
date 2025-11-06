# DNABERT-2 Baselines for ChickenGTEx

Two methods for evaluating DNABERT-2 on chicken splicing variants:

## 1. Logistic Regression (Zero-shot)
Uses differential embeddings (ALT - REF) with Logistic Regression classifier.

**Run on server:**
```bash
cd /mnt/userdata4/splicing/DNABert_2_chicken
bash run_dnabert2_logistic_chicken_server.sh
```

**Output:** `results/dnabert2_chicken_logistic_scores.tsv`

## 2. MLP (Embedding + Classifier)
Uses concatenated embeddings [REF, ALT, diff] (2304D) with MLP classifier.

**Run on server:**
```bash
cd /mnt/userdata4/splicing/DNABert_2_chicken
bash run_dnabert2_chicken_server_mlp.sh
```

**Output:** `results/dnabert2_chicken_mlp_variant_scores.tsv`

## Setup (First time only)
```bash
# Upload scripts to server
scp -r /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/chickenGTEx/baselines/DNABert_2/* \
    username@server:/mnt/userdata4/splicing/DNABert_2_chicken/

# Verify benchmark file exists
ls -lh /mnt/userdata4/splicing/ChickenGTEx/processed_data/chickengtex_silver_benchmark_balanced.tsv
```

## Notes
- Both methods use DNABERT-2 (117M parameters): `zhihan1996/DNABERT-2-117M`
- Input: 25,000 chicken variants (12,500 positive + 12,500 negative)
- Logistic method is faster but MLP typically achieves better performance
- MLP trains on 80% data, tests on 20% for performance evaluation

