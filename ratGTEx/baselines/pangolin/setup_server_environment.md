# Pangolin Server Environment Setup Guide

This guide provides step-by-step instructions for setting up the Pangolin environment on the server for RatGTEx evaluation.

## Server Directory Structure

```
/mnt/userdata4/splicing/
├── conda_envs/
│   └── pangolin-env/              # Pangolin conda environment
├── references/
│   └── rat/
│       ├── ratNor6.fa             # Rat reference genome
│       ├── ratNor6.fa.fai         # FASTA index
│       ├── ratNor6.gtf            # Rat gene annotation
│       └── ratNor6.annotation.db  # Pangolin annotation database
├── ratgtex/
│   └── processed_data/
│       └── ratgtex_silver_benchmark_balanced.tsv
└── pangolin/
    └── Pangolin-main/             # Pangolin source code
```

## Step 1: Create Conda Environment

```bash
# Navigate to splicing directory
cd /mnt/userdata4/splicing/

# Create conda environment for Pangolin (following official README requirements)
conda create -p ./conda_envs/pangolin-env python=3.8 -y

# Activate environment
conda activate ./conda_envs/pangolin-env

# Step 1: Install PyTorch FIRST (as per official README)
# For CPU-only (adjust for GPU if available)
conda install pytorch torchvision torchaudio cpuonly -c pytorch

# For GPU (if available on server)
# conda install pytorch torchvision torchaudio pytorch-cuda=11.8 -c pytorch -c nvidia

# Step 2: Install specific dependencies as per README
conda install -c conda-forge pyvcf
pip install gffutils biopython pandas pyfastx

# Step 3: Install additional packages for our evaluation pipeline
pip install scikit-learn matplotlib seaborn

# Step 4: Clone and install Pangolin from source (as per README)
cd /mnt/userdata4/splicing/
git clone https://github.com/tkzeng/Pangolin.git pangolin/Pangolin-source
cd pangolin/Pangolin-source
conda run -p /mnt/userdata4/splicing/conda_envs/pangolin-env pip install .

# Verify installation
conda run -p /mnt/userdata4/splicing/conda_envs/pangolin-env pangolin --help
```

## Step 2: Use Existing Rat Reference Files on Server

```bash
# Create references directory
mkdir -p /mnt/userdata4/splicing/references/rat/
cd /mnt/userdata4/splicing/references/rat/

# Create symbolic links to existing rat reference files on server
# (No need to download - files already exist)

# Link to existing rat genome FASTA
ln -s /mnt/userdata4/splicing/ratgetx/reference_genome/rn6/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa ratNor6.fa

# Link to existing FASTA index  
ln -s /mnt/userdata4/splicing/ratgetx/reference_genome/rn6/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.fai ratNor6.fa.fai

# Link to existing GTF annotation (keep original name as requested)
ln -s /mnt/userdata4/splicing/ratgetx/reference_genome/Rattus_norvegicus.Rnor_6.0.84.gtf ./Rattus_norvegicus.Rnor_6.0.84.gtf

# Verify files exist and are accessible
ls -lh ratNor6.fa*
ls -lh Rattus_norvegicus.Rnor_6.0.84.gtf
file ratNor6.fa
head -5 ratNor6.fa
head -20 Rattus_norvegicus.Rnor_6.0.84.gtf
```

## Step 3: Verify File Compatibility

```bash
# Check chromosome naming consistency between FASTA and GTF
echo "=== FASTA chromosomes ==="
grep "^>" ratNor6.fa | head -10

echo "=== GTF chromosomes ==="
cut -f1 Rattus_norvegicus.Rnor_6.0.84.gtf | grep -v "^#" | sort | uniq | head -10

# Verify versions match (both should be Rnor_6.0)
echo "=== File versions ==="
echo "FASTA file: $(basename ratNor6.fa)"
echo "GTF file: $(basename Rattus_norvegicus.Rnor_6.0.84.gtf)"
```

## Step 4: Create Pangolin Annotation Database

```bash
# Create annotation database using Pangolin
# 1. 直接进入ratgetx的reference_genome目录
cd /mnt/userdata4/splicing/ratgetx/reference_genome/

# 2. 激活Pangolin环境
conda activate /mnt/userdata4/splicing/conda_envs/pangolin-env

# 3. 检查所有需要的文件都在这里
ls -la Rattus_norvegicus.Rnor_6.0.84.gtf
ls -la rn6/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa*

# 4. 直接在这里创建注释数据库
python -c "
import gffutils
print('Creating annotation database...')
db = gffutils.create_db('Rattus_norvegicus.Rnor_6.0.84.gtf', 
                       'ratNor6.annotation.db', 
                       force=True,
                       keep_order=True)
print('✅ Annotation database created successfully')
"

# 5. 验证创建的数据库
ls -lh ratNor6.annotation.db
```

## Step 5: Verify Pangolin Installation

```bash
# Test Pangolin installation
conda run -p /mnt/userdata4/splicing/conda_envs/pangolin-env python -c "
import pangolin
print('✅ Pangolin module imported successfully')
"

# Test command-line tool
conda run -p /mnt/userdata4/splicing/conda_envs/pangolin-env pangolin --help

# Check PyTorch installation
conda run -p /mnt/userdata4/splicing/conda_envs/pangolin-env python -c "
import torch
print(f'✅ PyTorch version: {torch.__version__}')
print(f'✅ CUDA available: {torch.cuda.is_available()}')
"

# Verify all required dependencies
conda run -p /mnt/userdata4/splicing/conda_envs/pangolin-env python -c "
import gffutils, Bio, pandas, pyfastx, vcf
print('✅ All required dependencies available')
"
```

## Step 6: Verify RatGTEx Data on Server

```bash
# RatGTEx data already exists on server - just verify access
echo "=== Verifying RatGTEx Data ==="
ls -lh /mnt/userdata4/splicing/ratgetx/processed_data/ratgtex_silver_benchmark_balanced.tsv

# Check data format and size
echo "=== Data Statistics ==="
echo "File size: $(du -h /mnt/userdata4/splicing/ratgetx/processed_data/ratgtex_silver_benchmark_balanced.tsv | cut -f1)"
echo "Line count: $(wc -l < /mnt/userdata4/splicing/ratgetx/processed_data/ratgtex_silver_benchmark_balanced.tsv)"
echo "Column count: $(head -1 /mnt/userdata4/splicing/ratgetx/processed_data/ratgtex_silver_benchmark_balanced.tsv | tr '\t' '\n' | wc -l)"

# Show first few lines to verify format
echo "=== Data Format Preview ==="
head -3 /mnt/userdata4/splicing/ratgetx/processed_data/ratgtex_silver_benchmark_balanced.tsv | cut -c1-100
echo "... (sequences truncated for display)"
```

## Step 7: Transfer Pangolin Scripts

```bash
# Create scripts directory
mkdir -p /mnt/userdata4/splicing/scripts/pangolin/

# Transfer all Pangolin scripts from local machine
# (Run this from your local machine)
scp -r /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/ratGTEx/baselines/pangolin/* \
    your_server:/mnt/userdata4/splicing/scripts/pangolin/

# Make scripts executable
chmod +x /mnt/userdata4/splicing/scripts/pangolin/*.sh
chmod +x /mnt/userdata4/splicing/scripts/pangolin/*.py
```

## Step 8: Verify Complete Setup

```bash
# Test environment activation
conda activate /mnt/userdata4/splicing/conda_envs/pangolin-env

# Test Pangolin installation
pangolin --help

# Test file access
ls -la /mnt/userdata4/splicing/references/rat/
ls -la /mnt/userdata4/splicing/ratgtex/processed_data/

# Test Python imports
python -c "
import pandas as pd
import numpy as np
from pangolin import Pangolin
print('✅ All dependencies available')
"
```

## Step 9: Run Test Pipeline

```bash
# Navigate to scripts directory
cd /mnt/userdata4/splicing/scripts/pangolin/

# Run the complete pipeline
./run_pangolin_ratgtex_pipeline.sh

# Monitor progress
tail -f pangolin_pipeline_*.log
```

## Expected File Sizes

After setup, expect these approximate file sizes:

```
ratNor6.fa              ~2.8GB (rat genome)
ratNor6.fa.fai          ~125KB (FASTA index)
ratNor6.gtf             ~150MB (gene annotation)
ratNor6.annotation.db   ~50MB (Pangolin database)
ratgtex_silver_benchmark_balanced.tsv  ~200MB (RatGTEx data)
```

## Troubleshooting

### Issue: Database creation fails
```bash
# Check GTF format
head -50 ratNor6.gtf
grep -c "^#" ratNor6.gtf
grep -c "^[^#]" ratNor6.gtf

# Try with simpler GTF format
# Sometimes Ensembl GTF needs preprocessing
```

### Issue: Memory problems during database creation
```bash
# Monitor memory usage
free -h
htop

# Consider using a compute node with more memory
# or filtering GTF to essential features only
```

### Issue: Pangolin command not found
```bash
# Verify conda environment
which python
which pangolin
conda list | grep pangolin

# Reinstall if needed
pip uninstall pangolin
pip install pangolin
```

### Issue: Reference genome version mismatch
```bash
# Check chromosome naming
head -1 ratNor6.fa
grep "^>" ratNor6.fa | head -10

# Ensure consistency between FA and GTF
awk 'NR==1{print $1}' ratNor6.gtf
```

## Performance Optimization

For large-scale runs:

1. **Use SSD storage** for reference files
2. **Increase memory allocation** for Pangolin
3. **Consider batch processing** for large variant sets
4. **Monitor disk space** during runs

## Maintenance

Regular maintenance tasks:

```bash
# Update conda packages
conda update -p /mnt/userdata4/splicing/conda_envs/pangolin-env --all

# Check disk usage
du -sh /mnt/userdata4/splicing/

# Clean up temporary files
find /mnt/userdata4/splicing/ -name "*.tmp" -delete
find /mnt/userdata4/splicing/ -name "core.*" -delete
```

## Running the Pipeline

Once setup is complete:

```bash
# Navigate to scripts
cd /mnt/userdata4/splicing/scripts/pangolin/

# Run complete pipeline
./run_pangolin_ratgtex_pipeline.sh

# Check results
ls -la pangolin_ratgtex_results.json
ls -la pangolin_performance_metrics.json
```

The pipeline should complete successfully and generate all expected output files for cross-species performance analysis.


human


# 1. 中断当前进程
# 按 Ctrl+C

# 2. 删除未完成的数据库文件
rm -f /mnt/userdata4/splicing/SpliceAI/reference_genome/hg38.annotation.db

# 3. 运行优化版本
cd /mnt/userdata4/splicing/SpliceAI/reference_genome/
python -c "
import gffutils
print('Creating human annotation database with optimizations...')
db = gffutils.create_db(
    '/mnt/userdata4/splicing/SpliceAI/reference_genome/gencode.v38.annotation.gtf',
    'hg38.annotation.db',
    force=True,
    keep_order=True,
    merge_strategy='merge',
    disable_infer_genes=True,
    disable_infer_transcripts=True
)
print('Database created successfully!')
"


pangolin /mnt/userdata4/splicing/SpliceAI/data/spliceai_input_splicevardb.vcf \
         /mnt/userdata4/splicing/SpliceAI/reference_genome/hg38.fa \
         /mnt/userdata4/splicing/SpliceAI/reference_genome/hg38.annotation.db \
         /mnt/userdata4/splicing/pangolin/results/pangolin_splicevardb_full
