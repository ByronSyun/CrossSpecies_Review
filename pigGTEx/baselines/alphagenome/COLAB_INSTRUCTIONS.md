# AlphaGenome PigGTEx Prediction - Colab Instructions

## 概述

AlphaGenome 需要在 Google Colab 上运行，因为它依赖 Google Cloud API。本文档详细说明如何准备数据和运行预测。

## 步骤 1：准备 16384bp 序列数据（服务器端）

AlphaGenome 需要 16384bp 的序列长度（而不是我们标准的 8192bp）。

### 1.1 上传脚本到服务器

```bash
# 本地执行
scp /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/pigGTEx/baselines/alphagenome/prepare_16384bp_data.py \
    yinuos@mlerp-monash-node04:/mnt/userdata4/splicing/PigGTEx/
```

### 1.2 在服务器上运行

```bash
# 服务器执行
cd /mnt/userdata4/splicing/PigGTEx
python3 prepare_16384bp_data.py
```

**预期输出：**
- 文件：`/mnt/userdata4/splicing/PigGTEx/processed_data/piggtex_silver_benchmark_balanced_len16384.tsv`
- 格式：5列，无header
  - Column 1: `variant_id` (e.g., `1:12345_A/G`)
  - Column 2: `ref_seq` (16384bp)
  - Column 3: `alt_seq` (16384bp)
  - Column 4: `label` (0 or 1)
  - Column 5: `tissue_id`

### 1.3 下载到本地

```bash
# 本地执行
scp yinuos@mlerp-monash-node04:/mnt/userdata4/splicing/PigGTEx/processed_data/piggtex_silver_benchmark_balanced_len16384.tsv \
    /Users/byronsun/Desktop/AS_复现模型/BIB_review/data/processed_data/pigGTEx/
```

## 步骤 2：上传文件到 Google Drive

### 2.1 在 Google Drive 中创建文件夹

在 Google Drive 中创建以下路径：
```
My Drive/Colab Notebooks/pigGTEx_alphagenome/
```

### 2.2 上传文件

将以下文件上传到该文件夹：

1. **数据文件**（必需）：
   - `piggtex_silver_benchmark_balanced_len16384.tsv` (约 500MB)

2. **Python 脚本**（必需）：
   - `predict_piggtex_alphagenome.py`

文件路径：
```
/Users/byronsun/Desktop/AS_复现模型/BIB_review/data/processed_data/pigGTEx/piggtex_silver_benchmark_balanced_len16384.tsv
/Users/byronsun/Desktop/AS_复现模型/BIB_review/code/pigGTEx/baselines/alphagenome/predict_piggtex_alphagenome.py
```

## 步骤 3：在 Colab 中运行

### 3.1 创建新的 Colab Notebook

访问：https://colab.research.google.com

### 3.2 安装 AlphaGenome

在第一个 cell 中运行：

```python
!pip install alphagenome
```

### 3.3 上传并运行脚本

**方法 A：直接复制代码**

将 `predict_piggtex_alphagenome.py` 的内容复制粘贴到 Colab cell 中。

**方法 B：从 Drive 加载**

```python
# Mount Google Drive
from google.colab import drive
drive.mount('/content/drive')

# Run the script
%run "/content/drive/MyDrive/Colab Notebooks/pigGTEx_alphagenome/predict_piggtex_alphagenome.py"
```

### 3.4 配置 API Key

在运行前，修改脚本中的 API_KEY：

```python
API_KEY = "YOUR_ALPHAGENOME_API_KEY_HERE"
```

如果没有 API key，请访问：https://alphagenome.com/

### 3.5 运行预测

脚本会自动：
1. 挂载 Google Drive
2. 加载 PigGTEx 数据（~26,358 variants）
3. 分3批运行预测（避免网络中断）
4. 保存中间结果到 Drive
5. 合并所有批次结果

**预计运行时间：** 6-12小时（取决于网络和 API 速度）

## 步骤 4：下载结果

### 4.1 结果文件位置

在 Google Drive 中，结果保存在：
```
My Drive/Colab Notebooks/pigGTEx_alphagenome/
```

**输出文件：**
- `piggtex_alphagenome_batch_1_YYYYMMDD_HHMMSS.json` (批次1)
- `piggtex_alphagenome_batch_2_YYYYMMDD_HHMMSS.json` (批次2)
- `piggtex_alphagenome_batch_3_YYYYMMDD_HHMMSS.json` (批次3)
- `piggtex_alphagenome_complete_results_YYYYMMDD_HHMMSS.json` (合并结果)

### 4.2 下载到本地

从 Google Drive 下载 `piggtex_alphagenome_complete_results_*.json` 到：
```
/Users/byronsun/Desktop/AS_复现模型/BIB_review/code/pigGTEx/baselines/alphagenome/results/
```

## 结果格式

完整结果 JSON 文件结构：

```json
{
  "metadata": {
    "total_batches": 3,
    "final_timestamp": "20250110_123456",
    "total_variants": 26358,
    "total_successful": 26300,
    "total_failed": 58
  },
  "results": [
    {
      "variant_id": "1:12345_A/G",
      "label": 1,
      "status": "success",
      "alphagenome_score": 0.123,
      "error": null
    },
    ...
  ]
}
```

## 故障恢复

如果 Colab 中断，脚本会自动检测已完成的批次并跳过它们：

1. 重新运行脚本
2. 脚本会检查 Drive 中已存在的批次文件
3. 只运行缺失的批次
4. 运行完成后，可以手动调用 `merge_all_batches_fixed()` 来合并结果

## 常见问题

### Q1: AlphaGenome API 配额不足

**解决：** 降低 `NUM_BATCHES` 或分多天运行

### Q2: Colab 超时断开连接

**解决：** 脚本支持断点续传，重新运行即可

### Q3: 序列长度不是 16384bp

**解决：** 确认 `prepare_16384bp_data.py` 运行成功，检查输出文件

## 性能预期

基于 RatGTEx 的经验：
- **Coverage**: 100%（AlphaGenome 无需注释）
- **AUROC**: 约 0.65-0.70（跨物种零样本预测）

## 注意事项

1. **API Cost**: AlphaGenome API 可能收费，请确认配额
2. **运行时间**: 约 6-12 小时，建议使用 Colab Pro
3. **Drive 空间**: 确保至少有 1GB 可用空间存储结果
4. **批次安全**: 脚本自动分批并保存，避免数据丢失

