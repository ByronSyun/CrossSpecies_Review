# Server Execution Commands for Evo2 Rat Evaluation

## Quick Start (Copy-Paste Commands)

### 1. Upload Files to Server

```bash
# From local machine terminal
cd /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/ratGTEx/baselines/Evo2

# Upload prediction script
scp predict_evo2_rat.py user@server:/mnt/userdata4/splicing/Evo2_rat/

# Upload execution script
scp run_evo2_rat_server.sh user@server:/mnt/userdata4/splicing/Evo2_rat/

# Upload input data (if not already on server)
scp /Users/byronsun/Desktop/AS_复现模型/BIB_review/data/processed_data/ratGTEx/ratgtex_silver_benchmark_balanced.tsv \
    user@server:/mnt/userdata4/splicing/Evo2_rat/
```

### 2. SSH to Server and Run

```bash
# SSH to server
ssh user@server

# Navigate to working directory
cd /mnt/userdata4/splicing/Evo2_rat

# Check files are uploaded
ls -lh predict_evo2_rat.py run_evo2_rat_server.sh ratgtex_silver_benchmark_balanced.tsv

# Make script executable
chmod +x run_evo2_rat_server.sh

# Run prediction (this will take 3-6 hours)
bash run_evo2_rat_server.sh
```

**Alternative: Run in background with nohup**

```bash
# Run in background and save log
nohup bash run_evo2_rat_server.sh > evo2_rat.log 2>&1 &

# Get the process ID
echo $!

# Monitor progress
tail -f evo2_rat.log

# Check if still running
ps aux | grep run_evo2_rat_server.sh
```

### 3. Download Results

```bash
# From local machine
cd /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/ratGTEx/baselines/Evo2
mkdir -p results

# Download prediction results
scp user@server:/mnt/userdata4/splicing/Evo2_rat/evo2_rat_predictions.tsv results/

# Verify download
ls -lh results/evo2_rat_predictions.tsv
wc -l results/evo2_rat_predictions.tsv  # Should show ~28,121 lines (28,120 + header)
```

### 4. Run Local Evaluation

```bash
# On local machine
cd /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/ratGTEx/baselines/Evo2

python3 evaluate_evo2_rat.py \
    --predictions results/evo2_rat_predictions.tsv \
    --output_dir results

# Check results
ls -lh results/
cat results/evo2_rat_results.json
open results/evo2_rat_curves.png  # macOS
```

## Monitoring Server Job

### Check Progress

```bash
# SSH to server
ssh user@server

# Check if process is running
ps aux | grep predict_evo2_rat.py

# Monitor GPU usage (should show Evo2 model using ~40GB)
watch -n 2 nvidia-smi

# Check output file size (should be growing)
ls -lh /mnt/userdata4/splicing/Evo2_rat/evo2_rat_predictions.tsv
watch -n 10 'wc -l /mnt/userdata4/splicing/Evo2_rat/evo2_rat_predictions.tsv'

# View recent log output (if using nohup)
tail -50 /mnt/userdata4/splicing/Evo2_rat/evo2_rat.log
```

### Estimated Runtime

- **Total variants**: 28,120
- **Batch size**: 1 (conservative for 8kb sequences)
- **Time per batch**: ~1-2 seconds
- **Total time**: ~8-16 hours (depends on GPU load)

**Progress indicators**:
- Every 100 variants: ~2-3 minutes
- 1,000 variants: ~20-40 minutes
- 10,000 variants: ~3-7 hours
- Full dataset: ~8-16 hours

## Troubleshooting

### If prediction fails

```bash
# Check error log
tail -100 /mnt/userdata4/splicing/Evo2_rat/evo2_rat.log

# Verify input file exists and is readable
head -5 /mnt/userdata4/splicing/Evo2_rat/ratgtex_silver_benchmark_balanced.tsv

# Check container exists
ls -lh /mnt/userdata4/splicing/Evo2Splicevardb/vortex.sif

# Check model checkpoint exists
ls -lh /mnt/userdata4/splicing/Evo2Splicevardb/model_checkpoint/evo2_7b_base.pt

# Check GPU availability
nvidia-smi
```

### If OOM (Out of Memory) occurs

The script is already set to batch_size=1, which should be safe. If OOM still occurs:

```bash
# Check if other processes are using GPU
nvidia-smi

# Kill if necessary (replace PID with actual process ID)
kill -9 <PID>

# Wait for GPU to free up
watch -n 5 nvidia-smi

# Restart prediction
cd /mnt/userdata4/splicing/Evo2_rat
bash run_evo2_rat_server.sh
```

## File Verification

### Check input data format

```bash
# On server
cd /mnt/userdata4/splicing/Evo2_rat

# Count lines (should be 28,120, no header)
wc -l ratgtex_silver_benchmark_balanced.tsv

# Check first few lines
head -3 ratgtex_silver_benchmark_balanced.tsv

# Count labels (should be ~50/50 balanced)
awk '{print $4}' ratgtex_silver_benchmark_balanced.tsv | sort | uniq -c

# Check sequence lengths (should be 8192)
head -1 ratgtex_silver_benchmark_balanced.tsv | awk '{print length($2), length($3)}'
```

### Check output data format

```bash
# After prediction completes
cd /mnt/userdata4/splicing/Evo2_rat

# Check output file
head -10 evo2_rat_predictions.tsv

# Count lines (should be 28,121 = 28,120 + header)
wc -l evo2_rat_predictions.tsv

# Quick statistics
awk 'NR>1 {sum+=$4; sumsq+=$4*$4} END {print "Mean delta_logp:", sum/(NR-1); print "Std:", sqrt(sumsq/(NR-1) - (sum/(NR-1))^2)}' \
    evo2_rat_predictions.tsv
```

## Clean Up (After Successful Completion)

```bash
# On server (optional, to save space)
cd /mnt/userdata4/splicing/Evo2_rat

# Keep only essential files, remove scripts and logs
rm -f predict_evo2_rat.py run_evo2_rat_server.sh evo2_rat.log

# Compress results if needed
gzip evo2_rat_predictions.tsv

# Later, to use compressed file
gunzip evo2_rat_predictions.tsv.gz
```

## Expected Output Example

### Prediction file format (evo2_rat_predictions.tsv)

```
variant_id      logp_ref        logp_alt        delta_logp      label
15:55290237_A/G -2345.678       -2346.123       -0.445          1
7:3089418_A/T   -2890.123       -2889.987       0.136           1
9:25410298_T/A  -3012.456       -3012.501       -0.045          0
...
```

### Evaluation results (evo2_rat_results.json)

```json
{
  "n_variants": 28120,
  "n_positive": 14060,
  "n_negative": 14060,
  "positive_rate": 0.5,
  "auroc": 0.XXXX,
  "auprc": 0.XXXX,
  "pearson_r": 0.XXXX,
  "pearson_p": 0.XXXX
}
```

## Notes

- **Batch size**: Set to 1 for safety with 8kb sequences and 7B model
- **GPU memory**: Expects ~40GB (Tesla A100 or similar)
- **Runtime**: Highly variable depending on server load (3-16 hours)
- **Checkpointing**: Not implemented; if interrupted, must restart from beginning

