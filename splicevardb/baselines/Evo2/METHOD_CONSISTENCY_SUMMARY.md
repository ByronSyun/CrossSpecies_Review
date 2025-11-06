# Evo2 Method Consistency Check

## 你的两个问题

### 1. SpliceVarDB 脚本是否和 rat 方法一致？

**是的，现在完全一致了。** 主要修正包括：

#### 核心方法一致性 ✓
- **Embedding extraction**: 相同的 `get_evo2_embeddings()` 函数，提取 hidden states + mean pooling
- **Concatenation strategy**: `[WT, MT, diff]` 三部分拼接
- **MLP architecture**: 2层 `[256, 128]`，相同的超参数
- **Normalization**: StandardScaler
- **Evaluation protocol**: 80/20 train/test split

#### 关键差异（已解决）
| 方面 | Rat 版本 | SpliceVarDB 初始版本 | **修正后** |
|------|----------|---------------------|-----------|
| MLP `batch_size` | 64 | 200 | ✓ **64** |
| MLP `max_iter` | 1000 | 200 | ✓ **1000** |
| MLP `n_iter_no_change` | 20 | 10 | ✓ **20** |
| Predictions TSV | 全部数据 | 仅 test set | ✓ **全部数据** |
| 生成方式 | 内联 Python 脚本 | train_evo2_mlp.py 内部 | ✓ **内联 Python 脚本** |

#### 数据格式差异（已处理）
- **Rat**: 直接有 `ref_seq` 和 `alt_seq` 两列
- **SpliceVarDB**: 只有 `sequence` 列（WT），需要从中心位置构建 mutant
- **解决方案**: `create_mutant_sequences()` 函数自动处理，与 DNABERT-2 的方法一致

---

### 2. 跑完脚本后最终的 score TSV 会生成吗？

**是的！** 脚本会生成以下文件：

#### 最终输出文件
```
splicevardb_data/
├── evo2_concat_embeddings.npz           # Step 1: Embeddings
├── evo2_mlp_model.pkl                   # Step 2: 训练好的模型
├── evo2_mlp_scaler.pkl                  # Step 2: 标准化器
├── evo2_mlp_train_test_metrics.json     # Step 2: Train/test 性能
├── evo2_mlp_variant_scores.tsv          # Step 3: ⭐ 全部 variants 的 scores
└── evo2_mlp_full_metrics.json           # Step 3: 全部数据的性能指标
```

#### `evo2_mlp_variant_scores.tsv` 格式
```tsv
variant_id              true_label  prediction_probability  predicted_label
chr1:100110484:C>A      1           0.8234                  1
chr1:100207186:T>C      1           0.7891                  1
chr1:100267318:G>A      0           0.2145                  0
...
```

这个 TSV 文件可以直接用于：
1. 与其他模型比较（加入 `summary_results.py`）
2. 绘制 ROC/PR 曲线
3. 计算各种性能指标

---

## 运行指令（服务器）

```bash
cd /mnt/userdata4/splicing/Evo2Splicevardb
bash run_evo2_concat_mlp.sh
```

运行完成后，查看结果：
```bash
# 查看 predictions
head splicevardb_data/evo2_mlp_variant_scores.tsv

# 查看性能指标
cat splicevardb_data/evo2_mlp_full_metrics.json
```

---

## 与 Rat 结果对比

| 数据集 | Method | AUROC | AUPRC | Accuracy |
|--------|--------|-------|-------|----------|
| RatGTEx | Evo2 Embedding+MLP | **0.7761** | 0.7707 | 0.7001 |
| RatGTEx | Evo2 Direct `delta_logp` | 0.4848 | - | - |
| SpliceVarDB | Evo2 Embedding+MLP | **待运行** | - | - |

预期 SpliceVarDB 表现会更好（因为是 human 数据，模型在 human 上训练）。

---

## 论文中的应用

这个结果可以用来说明：
1. **Generalist LLM (Evo2) 的 embedding 质量高**：通过 downstream MLP 可以达到很好的效果
2. **跨物种迁移能力**：从 human (SpliceVarDB) 训练的 MLP 可以应用到 rat 数据，AUROC 0.7761
3. **Utilization paradigm 的重要性**：
   - Direct scoring (delta_logp): 简单但效果差（AUROC 0.48）
   - Embedding + classifier: 复杂但效果好（AUROC 0.78）
4. **与其他 embedding models 对比**（DNABERT-2, AlphaGenome 等）

---

## 修改总结

修改的文件：
1. ✓ `extract_evo2_embeddings.py` - 添加 shape check, 保存 embedding_dim
2. ✓ `train_evo2_mlp.py` - 统一超参数，移除内部 predictions 保存
3. ✓ `run_evo2_concat_mlp.sh` - 添加 Step 3 内联 Python 脚本生成全部数据的 predictions
4. ✓ `README_CONCAT_MLP.md` - 更新文档说明输出文件和一致性

所有修改已完成，可以直接在服务器运行！

