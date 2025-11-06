# SpliceBERT 服务器运行指令

## 🔧 重要更新 (2025-09-21)

**关键修复**: 已修复重大的KL-context分数计算错误，现在严格遵循SpliceBERT原论文的实现方法：

- ✅ **正确的KL距离计算**：使用完整序列KL散度而不是窗口片段
- ✅ **距离加权平均**：对变异位点前后对称位置的KL值进行平均 `(before_reversed + after) / 2`
- ✅ **正确的log-sum公式**：`sum(log(clip(kl_distance, 1e-6)))` 而不是 `sum(log(KL_scores))`

**预期改进**: 修复后的实现应该显著提升SpliceBERT在SpliceVarDB上的性能表现。

---

## 模型下载状态检查
```bash
cd /mnt/userdata4/splicing/SpliceBERT
ls -la models.tar.gz
# 检查下载进度，如果还在下载会显示当前大小
```

## 等模型下载完成后执行
```bash
# 解压模型
tar -zxvf models.tar.gz
ls -la models/

# 验证模型结构
find models/ -name "*.json" -o -name "*.bin" -o -name "*.safetensors" | head -10
```

## 上传脚本到服务器
将以下文件上传到服务器的 `/mnt/userdata4/splicing/SpliceBERT/` 目录：
- `splicebert_score_variants.py` - SpliceBERT变异评分脚本
- `evaluate_splicebert_scores.py` - 评估脚本

## 服务器端运行命令

### 步骤1：进入工作目录并设置权限
```bash
cd "/mnt/userdata4/splicing/SpliceBERT"
chmod +x splicebert_score_variants.py
chmod +x evaluate_splicebert_scores.py
mkdir -p results
```

### 步骤2：测试SpliceBERT模型加载
```bash
conda run -p /mnt/userdata4/splicing/conda_envs/splicebert-env python -c "
from transformers import AutoTokenizer, AutoModelForMaskedLM
import torch

model_path = '/mnt/userdata4/splicing/SpliceBERT/models/model_folder'  # 根据实际解压后的路径调整
print('Testing SpliceBERT model loading...')
try:
    tokenizer = AutoTokenizer.from_pretrained(model_path)
    model = AutoModelForMaskedLM.from_pretrained(model_path)
    print('✅ SpliceBERT model loaded successfully!')
    print(f'Vocab size: {tokenizer.vocab_size}')
    print(f'Model parameters: {sum(p.numel() for p in model.parameters())/1e6:.1f}M')
except Exception as e:
    print(f'❌ Error loading model: {e}')
    print('Available model directories:')
    import os
    for item in os.listdir('/mnt/userdata4/splicing/SpliceBERT/models/'):
        print(f'  - {item}')
"
```

### 步骤3：运行人类SpliceVarDB数据集
```bash
# 注意：需要根据实际的模型路径调整 --model_path 参数
conda run -p /mnt/userdata4/splicing/conda_envs/splicebert-env python splicebert_score_variants.py \
  --input_tsv "/mnt/userdata4/splicing/Evo2Splicevardb/splicevardb_data/evo2_splicevardb_dataset_dedup.tsv" \
  --model_path "/mnt/userdata4/splicing/SpliceBERT/models/model_folder" \
  --output_prefix "results/splicebert_splicevardb" \
  --flanking_window 100 \
  --batch_size 2
```

### 步骤4：评估人类数据结果
```bash
conda run -p /mnt/userdata4/splicing/conda_envs/splicebert-env python evaluate_splicebert_scores.py \
  --labels "/mnt/userdata4/splicing/Evo2Splicevardb/splicevardb_data/evo2_splicevardb_dataset_dedup.tsv" \
  --scores "results/splicebert_splicevardb_scores.csv" \
  --out_dir "results/splicevardb_eval"
```

### 步骤5：运行大鼠RatGTEx数据集
```bash
conda run -p /mnt/userdata4/splicing/conda_envs/splicebert-env python splicebert_score_variants.py \
  --input_tsv "/mnt/userdata4/splicing/ratgetx/processed_data/ratgtex_silver_benchmark_balanced.tsv" \
  --model_path "/mnt/userdata4/splicing/SpliceBERT/models/model_folder" \
  --output_prefix "results/splicebert_ratgtex" \
  --flanking_window 100 \
  --batch_size 2
```

### 步骤6：评估大鼠数据结果
```bash
conda run -p /mnt/userdata4/splicing/conda_envs/splicebert-env python evaluate_splicebert_scores.py \
  --labels "/mnt/userdata4/splicing/ratgetx/processed_data/ratgtex_silver_benchmark_balanced.tsv" \
  --scores "results/splicebert_ratgtex_scores.csv" \
  --out_dir "results/ratgtex_eval"
```

## 预期输出文件

### 人类数据集结果：
- `results/splicebert_splicevardb_scores.csv` - SpliceBERT KL-divergence评分结果
- `results/splicebert_splicevardb_scores.jsonl` - 评分结果(JSONL格式)
- `results/splicevardb_eval/splicebert_evaluation_results.csv` - 详细评估指标
- `results/splicevardb_eval/splicebert_evaluation_summary.json` - 评估摘要
- `results/splicevardb_eval/score_distribution.png` - 分数分布图
- `results/splicevardb_eval/roc_curve.png` - ROC曲线
- `results/splicevardb_eval/precision_recall_curve.png` - PR曲线

### 大鼠数据集结果：
- `results/splicebert_ratgtex_scores.csv` - SpliceBERT KL-divergence评分结果
- `results/splicebert_ratgtex_scores.jsonl` - 评分结果(JSONL格式)
- `results/ratgtex_eval/splicebert_evaluation_results.csv` - 详细评估指标
- `results/ratgtex_eval/splicebert_evaluation_summary.json` - 评估摘要
- `results/ratgtex_eval/score_distribution.png` - 分数分布图
- `results/ratgtex_eval/roc_curve.png` - ROC曲线
- `results/ratgtex_eval/precision_recall_curve.png` - PR曲线

## 技术说明

### SpliceBERT方法论：
1. **序列预处理**：DNA序列转RNA格式(T→U)，添加空格分隔符
2. **KL散度评分**：比较参考序列和变异序列的MLM预测分布
3. **Flanking窗口**：计算变异位点上下游±100nt内的KL散度总和
4. **评分公式**：KL-context score = Σ log(KL(alt||ref))

### 关键参数：
- `--flanking_window 100`：根据SpliceBERT论文设置
- `--batch_size 2`：保守设置，避免GPU内存不足
- `--window_size 1024`：SpliceBERT最大序列长度（会从8192bp中心截取1024nt）

## 故障排除

### 如果模型路径错误：
检查解压后的实际模型目录结构：
```bash
find /mnt/userdata4/splicing/SpliceBERT/models/ -name "config.json" -o -name "pytorch_model.bin" -o -name "model.safetensors"
```

### 如果GPU内存不足：
```bash
# 减少batch_size
--batch_size 1

# 或使用CPU（会很慢）
export CUDA_VISIBLE_DEVICES=""
```

### 如果序列长度超限：
SpliceBERT训练时使用64-1024nt长度，如果输入序列过长可能需要截断。
