# 重新运行修复后的SpliceBERT评估

## 🔧 修复说明

已修复SpliceBERT实现中的关键错误：
- 正确的KL距离计算方法
- 距离加权平均处理
- 遵循原论文的log-sum公式

## 服务器运行命令

### 步骤1：上传修复后的脚本
将以下文件上传到服务器 `/mnt/userdata4/splicing/SpliceBERT/`：
- `splicebert_score_variants.py` (修复后的版本)
- `evaluate_splicebert_scores.py`

### 步骤2：备份旧结果
```bash
cd /mnt/userdata4/splicing/SpliceBERT
mkdir -p results_backup_old
mv results/splicebert_splicevardb_* results_backup_old/ 2>/dev/null || echo "No old results to backup"
```

### 步骤3：重新运行SpliceBERT (人类数据)
```bash
# 激活环境
conda activate /mnt/userdata4/splicing/conda_envs/splicebert-env

# 重新运行人类SpliceVarDB数据集
conda run -p /mnt/userdata4/splicing/conda_envs/splicebert-env python splicebert_score_variants.py \
  --input_tsv "/mnt/userdata4/splicing/Evo2Splicevardb/splicevardb_data/evo2_splicevardb_dataset_dedup.tsv" \
  --model_path "/mnt/userdata4/splicing/SpliceBERT/SpliceBERT-main/models/models/SpliceBERT.1024nt" \
  --output_prefix "results/splicebert_splicevardb_fixed" \
  --flanking_window 100 \
  --batch_size 2 \
  --window_size 1024
```

### 步骤4：评估修复后的结果
```bash
# 评估人类数据结果
conda run -p /mnt/userdata4/splicing/conda_envs/splicebert-env python evaluate_splicebert_scores.py \
  --labels "/mnt/userdata4/splicing/Evo2Splicevardb/splicevardb_data/evo2_splicevardb_dataset_dedup.tsv" \
  --scores "results/splicebert_splicevardb_fixed_scores.csv" \
  --out_dir "results/splicevardb_fixed_eval"
```

### 步骤5：比较修复前后的结果
```bash
# 检查修复前后的分数分布差异
conda run -p /mnt/userdata4/splicing/conda_envs/splicebert-env python -c "
import pandas as pd
import numpy as np

print('=== 修复前后结果对比 ===')

# 加载旧结果 (如果存在)
try:
    old_scores = pd.read_csv('results_backup_old/splicebert_splicevardb_scores.csv')
    print(f'修复前: {len(old_scores)} 样本')
    print(f'  分数范围: [{old_scores[\"kl_context_score\"].min():.4f}, {old_scores[\"kl_context_score\"].max():.4f}]')
    print(f'  分数均值: {old_scores[\"kl_context_score\"].mean():.4f}')
    print(f'  NaN数量: {old_scores[\"kl_context_score\"].isna().sum()}')
except:
    print('修复前: 无旧结果文件')

# 加载新结果
try:
    new_scores = pd.read_csv('results/splicebert_splicevardb_fixed_scores.csv')
    print(f'修复后: {len(new_scores)} 样本')
    print(f'  分数范围: [{new_scores[\"kl_context_score\"].min():.4f}, {new_scores[\"kl_context_score\"].max():.4f}]')
    print(f'  分数均值: {new_scores[\"kl_context_score\"].mean():.4f}')
    print(f'  NaN数量: {new_scores[\"kl_context_score\"].isna().sum()}')
except:
    print('修复后: 新结果文件不存在')
"
```

### 步骤6：运行大鼠数据 (如果需要)
```bash
# 运行大鼠RatGTEx数据集
conda run -p /mnt/userdata4/splicing/conda_envs/splicebert-env python splicebert_score_variants.py \
  --input_tsv "/mnt/userdata4/splicing/ratgetx/processed_data/ratgtex_silver_benchmark_balanced.tsv" \
  --model_path "/mnt/userdata4/splicing/SpliceBERT/SpliceBERT-main/models/models/SpliceBERT.1024nt" \
  --output_prefix "results/splicebert_ratgtex_fixed" \
  --flanking_window 100 \
  --batch_size 2 \
  --window_size 1024

# 评估大鼠数据结果
conda run -p /mnt/userdata4/splicing/conda_envs/splicebert-env python evaluate_splicebert_scores.py \
  --labels "/mnt/userdata4/splicing/ratgetx/processed_data/ratgtex_silver_benchmark_balanced.tsv" \
  --scores "results/splicebert_ratgtex_fixed_scores.csv" \
  --out_dir "results/ratgtex_fixed_eval"
```

## 预期结果

修复后的SpliceBERT应该显示：
1. **更高的AUROC/AUPRC**: 预期从 ~0.52 提升到 >0.7
2. **更少的NaN值**: 更稳定的KL计算
3. **更合理的分数分布**: 符合原论文的表现

## 下载结果到本地

运行完成后，将以下文件下载到本地进行分析：
- `results/splicebert_splicevardb_fixed_scores.csv`
- `results/splicevardb_fixed_eval/splicebert_evaluation_results.json`
- `results/splicevardb_fixed_eval/*.png` (可视化图表)
