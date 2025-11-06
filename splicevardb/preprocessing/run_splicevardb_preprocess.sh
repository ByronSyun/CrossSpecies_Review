#!/bin/bash
set -euo pipefail

# Absolute project root
ROOT="/Users/byronsun/Desktop/AS_复现模型"

# Paths
SCRIPT="$ROOT/BIB_review/code/splicevardb/preprocessing/prepare_dataset_for_evo2.py"
DATA_TSV="$ROOT/BIB_review/data/datasets/splicevardb.download_all.tsv"
FASTA="$ROOT/BIB_review/data/reference_genome/hg38.fa"
OUT_DIR="$ROOT/BIB_review/data/processed_data/splicevardb"
MAIN_OUT="$OUT_DIR/evo2_splicevardb_dataset_dedup.tsv"
VARLIST_OUT="$OUT_DIR/splicevardb_filter_benchmark.tsv"

echo "--- SpliceVarDB preprocessing (deduplicated, positives+negatives, with variant list) ---"

# Ensure directories
mkdir -p "$(dirname "$DATA_TSV")" "$(dirname "$FASTA")" "$OUT_DIR"

# Optional: create dataset symlink if missing
if [ ! -f "$DATA_TSV" ]; then
  SRC_DATA="$ROOT/Evo2Splicevardb/splicevardb_data/splicevardb.download_all.tsv"
  if [ -f "$SRC_DATA" ]; then
    ln -s "$SRC_DATA" "$DATA_TSV"
    echo "Linked dataset TSV -> $DATA_TSV"
  else
    echo "ERROR: Dataset TSV not found at $DATA_TSV or $SRC_DATA" >&2
    exit 1
  fi
fi

# FASTA should be symlinked/placed beforehand
if [ ! -f "$FASTA" ]; then
  echo "ERROR: FASTA not found at $FASTA" >&2
  echo "Hint: ln -s \"$ROOT/reference_genome/hg38/Homo_sapiens.GRCh38.dna.primary_assembly.fa\" \"$FASTA\"" >&2
  exit 1
fi

# Run preprocessing with deduplication and variant list output
python "$SCRIPT" \
  --hg38_info "$DATA_TSV" \
  --fasta "$FASTA" \
  --output "$MAIN_OUT" \
  --dedup_by_hg38 \
  --variant_list_out "$VARLIST_OUT" | cat

# Line count summary
if [ -f "$VARLIST_OUT" ]; then
  echo "\n--- Variant list ---"
  echo "$VARLIST_OUT"
  echo "Lines (incl header): $(wc -l < "$VARLIST_OUT")"
fi

if [ -f "$MAIN_OUT" ]; then
  echo "\n--- Main output ---"
  echo "$MAIN_OUT"
  echo "Lines (incl header): $(wc -l < "$MAIN_OUT")"
else
  echo "ERROR: Main output file not created: $MAIN_OUT" >&2
  exit 1
fi

echo "--- Done ---"
