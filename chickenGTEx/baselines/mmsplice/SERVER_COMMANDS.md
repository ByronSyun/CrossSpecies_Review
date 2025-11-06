# MMSplice ChickenGTEx Server Commands

## Upload Script
```bash
cd /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/chickenGTEx/baselines/mmsplice

scp run_mmsplice_chicken.sh yinuos@mlerp-monash-node04:/mnt/userdata4/splicing/MMSplice/
```

## Run on Server
```bash
ssh yinuos@mlerp-monash-node04
cd /mnt/userdata4/splicing/MMSplice
bash run_mmsplice_chicken.sh
```

## Download Results
```bash
# Main results
scp yinuos@mlerp-monash-node04:/mnt/userdata4/splicing/MMSplice/chickengtex/mmsplice_chickengtex_scores.csv \
  /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/chickenGTEx/baselines/mmsplice/results/

# VCF (optional, for debugging)
scp yinuos@mlerp-monash-node04:/mnt/userdata4/splicing/MMSplice/chickengtex/chickengtex_input.vcf \
  /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/chickenGTEx/baselines/mmsplice/results/
```

## Expected Runtime
- VCF creation: ~30 seconds
- MMSplice predictions: ~2-3 hours (25,000 variants)
- Total: ~2-3 hours

## Expected Output
`mmsplice_chickengtex_scores.csv`:
- Rows: ~25,000 (one per variant)
- Key columns: ID, delta_logit_psi, pathogenicity, efficiency
- Size: ~5-10 MB

## Evaluate Results Locally
```bash
cd /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/chickenGTEx/baselines/mmsplice

python3 - << 'EOF'
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score, accuracy_score

# Load MMSplice predictions
mmsplice_df = pd.read_csv('results/mmsplice_chickengtex_scores.csv')

# Load ground truth labels
benchmark_df = pd.read_csv(
    '/mnt/userdata4/splicing/ChickenGTEx/processed_data/chickengtex_silver_benchmark_balanced.tsv',
    sep='\t', header=None,
    names=['variant_id', 'ref_seq', 'alt_seq', 'label', 'tissue_id']
)

# Merge on variant ID
merged = mmsplice_df.merge(benchmark_df[['variant_id', 'label']], 
                            left_on='ID', right_on='variant_id', how='inner')

print(f"ChickenGTEx MMSplice Results:")
print(f"  Total variants: {len(merged):,}")
print(f"  Coverage: {len(merged) / len(benchmark_df) * 100:.2f}%")

# Evaluate using delta_logit_psi (primary score)
y_true = merged['label']
y_score = merged['delta_logit_psi'].abs()  # Use absolute value for effect magnitude
y_pred = (y_score > y_score.median()).astype(int)

print(f"\nUsing delta_logit_psi (absolute value):")
print(f"  AUROC: {roc_auc_score(y_true, y_score):.4f}")
print(f"  AUPRC: {average_precision_score(y_true, y_score):.4f}")
print(f"  Accuracy: {accuracy_score(y_true, y_pred):.4f}")

# Also try pathogenicity score
if 'pathogenicity' in merged.columns:
    y_score_path = merged['pathogenicity']
    print(f"\nUsing pathogenicity score:")
    print(f"  AUROC: {roc_auc_score(y_true, y_score_path):.4f}")
    print(f"  AUPRC: {average_precision_score(y_true, y_score_path):.4f}")
EOF
```

## Troubleshooting

### If GTF not found
```bash
# Download GTF from Ensembl
cd /mnt/userdata4/splicing/ChickenGTEx/reference_genome/
wget http://ftp.ensembl.org/pub/release-102/gtf/gallus_gallus/Gallus_gallus.GRCg6a.102.gtf.gz
gunzip Gallus_gallus.GRCg6a.102.gtf.gz
mv Gallus_gallus.GRCg6a.102.gtf Gallus_gallus.GRCg6a.gtf
```

### If FASTA not found or mismatched
```bash
# Download Ensembl FASTA (matches GTF chromosome naming)
cd /mnt/userdata4/splicing/ChickenGTEx/reference_genome/
wget -c http://ftp.ensembl.org/pub/release-102/fasta/gallus_gallus/dna/Gallus_gallus.GRCg6a.dna.toplevel.fa.gz
gunzip Gallus_gallus.GRCg6a.dna.toplevel.fa.gz
samtools faidx Gallus_gallus.GRCg6a.dna.toplevel.fa
```

### If FASTA index missing
The script will create it automatically, or manually:
```bash
samtools faidx /mnt/userdata4/splicing/ChickenGTEx/reference_genome/Gallus_gallus.GRCg6a.dna.toplevel.fa
```

### If chromosome naming mismatch persists
1. Check VCF chromosome format: `grep -v "^#" chickengtex_input.vcf | head -1`
2. Check FASTA chromosome format: `grep "^>" Gallus_gallus.GRCg6a.dna.toplevel.fa | head -5`
3. Check GTF chromosome format: `head -20 Gallus_gallus.GRCg6a.gtf | grep -v "^#"`
4. All should use Ensembl format: `1`, `2`, `3`... (not `chr1`, `chr2`)

