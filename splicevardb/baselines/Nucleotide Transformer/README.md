# Nucleotide Transformer Baseline (Human: SpliceVarDB; Rat: Silver Benchmark)

This folder contains the zero-shot variant-effect baseline using the Nucleotide Transformer (NT) Genomic Foundation Model. We score variants by comparing sequence embeddings between the reference and the alternative alleles within a fixed window.

Important notes
- We only run on the server (GPU). Do NOT run deep learning locally.
- All environments live under `splicing/conda_envs/` on the server, per our convention.
- Outputs are stored under `splicing/Nucleotide Transformer/results/`.

## 1) Server Environment Setup

Create a dedicated conda environment (PyTorch + Hugging Face Transformers). We use PyTorch to avoid JAX installation complexity.


cd "/mnt/userdata4/splicing/Nucleotide Transformer"
git clone https://github.com/instadeepai/nucleotide-transformer.git
cd nucleotide-transformer

```bash
# Paths (server)
ENV_DIR=/mnt/userdata4/splicing/conda_envs/nt-v2
HF_CACHE=/mnt/userdata4/splicing/.cache/huggingface
RESULTS_DIR="/mnt/userdata4/splicing/Nucleotide Transformer/results"

mkdir -p $(dirname "$ENV_DIR")
mkdir -p "$HF_CACHE" "$RESULTS_DIR"

# 1. Create the environment
conda create -p "$ENV_DIR" python=3.10 -y

# 2. Install CUDA-enabled PyTorch (adjust CUDA wheel if needed)
conda run -p "$ENV_DIR" pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121

# 3. Install runtime deps
conda run -p "$ENV_DIR" pip install transformers tokenizers numpy pandas tqdm pyfaidx einops scikit-learn matplotlib seaborn

# 4. (Optional but recommended) set HF cache to shared path
export HF_HOME="$HF_CACHE"
```

GPU check:
```bash
conda run -p "$ENV_DIR" python - << 'PY'
import torch
print('CUDA available:', torch.cuda.is_available())
PY
```

## 2) Required Inputs

- Human (SpliceVarDB): TSV file with at least the following columns:
  - `variant_id, chrom, pos_0based, ref, alt`
  - Preferred: `ref_seq, alt_seq` (8,192 bp window centered at the variant). If absent, provide a reference FASTA to let the script build sequences on-the-fly.
- Rat (RatGTEx silver): same as above; FASTA must be rn6.

Recommended sequence window: 8,192 bp (fits NT-v2 context 12kb).

## 3) Model and Scoring Protocol

- Model: NT-v2 multispecies (e.g., 250M or 500M). We load from Hugging Face by id.
- Tokenisation: 6-mer tokenizer provided by the model.
- Pooling: mean-pool the last hidden state across tokens.
- Variant score: distance between embeddings of `alt_seq` and `ref_seq`.
  - Default: cosine distance `1 - cos(ref, alt)`.
  - Also outputs L2 distance.

Outputs
- JSON Lines (`.jsonl`) and CSV with columns: `variant_id, score_cosine, score_l2`.

## 4) Run Commands (Server)

Set paths (edit as needed):
```bash
ENV_DIR=/mnt/userdata4/splicing/conda_envs/nt-v2
RESULTS_DIR="/mnt/userdata4/splicing/Nucleotide Transformer/results"
HF_HOME=/mnt/userdata4/splicing/.cache/huggingface

mkdir -p "$RESULTS_DIR"
```

### 4.1 Human: SpliceVarDB
```bash
conda run -p "$ENV_DIR" python nt_score_variants.py \
  --model_id InstaDeepAI/nucleotide-transformer-v2-500m-multi-species \
  --input_tsv "/mnt/userdata4/splicing/Evo2Splicevardb/splicevardb_data/evo2_splicevardb_dataset_dedup.tsv" \
  --output_prefix "$RESULTS_DIR/nt_v2_splicevardb" \
  --window 8192 \
  --batch_size 4 \
  --num_workers 2
```

### 4.2 Rat: RatGTEx Silver
```bash
conda run -p "$ENV_DIR" python nt_score_variants.py \
  --model_id InstaDeepAI/nucleotide-transformer-v2-500m-multi-species \
  --input_tsv /mnt/userdata4/splicing/ratgetx/processed_data/ratgtex_silver_benchmark_balanced.tsv \
  --output_prefix "$RESULTS_DIR/nt_v2_ratgtex" \
  --window 8192 \
  --batch_size 4 \
  --num_workers 2
```

Notes
- If your TSV lacks `ref_seq/alt_seq`, add `--ref_fasta` and `--center_from_columns chrom pos_0based ref alt` to let the script generate sequences.

## 5) Evaluation

Reuse our generic evaluator to compute AUROC, AUPRC, confusion matrix and PR/ROC curves at fixed threshold 0.01.
```bash
conda run -p "$ENV_DIR" python evaluate_nt_scores.py \
  --labels "/mnt/userdata4/splicing/Evo2Splicevardb/splicevardb_data/evo2_splicevardb_dataset_dedup.tsv" \
  --scores "$RESULTS_DIR/nt_v2_splicevardb_scores.csv" \
  --out_dir "$RESULTS_DIR/splicevardb_eval"

conda run -p "$ENV_DIR" python evaluate_nt_scores.py \
  --labels /mnt/userdata4/splicing/ratgetx/processed_data/ratgtex_silver_benchmark_balanced.tsv \
  --scores "$RESULTS_DIR/nt_v2_ratgtex_scores.csv" \
  --out_dir "$RESULTS_DIR/ratgtex_eval"
```

## 6) Repro and Logging

- All scripts write a log header with model id, dataset, total variants and timing.
- Intermediate outputs:
  - `*_scores.jsonl` (streaming friendly)
  - `*_scores.csv` (merged view for downstream evaluation)

## 7) Caveats

- Sequence length must be ≤ 12,282 nt (NT‑v2). Our default 8,192 bp window is safe.
- When sequences contain many `N`, tokenisation splits into single letters; performance may degrade.
- Loss/Δloss zero-shot scoring is supported in principle but not implemented here for speed; embedding-distance works well and matches our cross-model protocol.
