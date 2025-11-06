# RatGTEx Performance Visualization

本目录包含 RatGTEx 跨物种基准测试的所有性能可视化图表。

## 生成的图表

### 一、Rat 单独性能图表

#### 1. ROC Curves (`rat_roc_curves.png/pdf`)
- **描述**: 所有9个模型在 rat 数据上的 ROC 曲线
- **关键发现**: AlphaGenome 表现最佳 (AUROC 0.600)，多数模型接近随机 (AUROC ≈ 0.50)

#### 2. Precision-Recall Curves (`rat_precision_recall_curves.png/pdf`)
- **描述**: 所有模型的 Precision-Recall 曲线
- **关键发现**: SpliceAI 和 Pangolin 的 AUPRC 相对较高，但覆盖率低

#### 3. Performance Comparison (`rat_performance_comparison.png/pdf`)
- **描述**: AUROC 和 AUPRC 的并排柱状图对比
- **关键发现**: 性能排名：AlphaGenome > SpliceAI > NT > Pangolin > 其他

#### 4. Coverage Analysis (`rat_coverage_analysis.png/pdf`)
- **描述**: 各模型的覆盖率横向柱状图
- **关键发现**: 
  - 100% 覆盖: SpliceBERT, Evo2, NT, DNABERT2_Logistic, SpliceTransformer
  - 低覆盖: Pangolin (38.47%), SpliceAI (38.44%), MMSplice (5.05%)

### 二、跨物种对比图表 (Human vs Rat)

#### 5. Cross-Species AUROC Comparison (`cross_species_auroc_comparison.png/pdf`)
- **描述**: Human (SpliceVarDB) 和 Rat (RatGTEx) 的 AUROC 并排对比
- **关键发现**: 所有模型在 rat 上性能都显著下降
- **颜色**: 蓝色 = Human, 橙色 = Rat

#### 6. Performance Drop Analysis (`cross_species_performance_drop.png/pdf`)
- **描述**: 从 Human 到 Rat 的性能下降百分比（横向柱状图）
- **关键发现**: 
  - 最大下降: DNABERT2_Logistic (100% 下降至随机)
  - 最小下降: AlphaGenome (~37% 下降)
- **颜色**: 每个模型使用其专属颜色（与其他图表一致）

#### 7. Coverage Comparison (`cross_species_coverage_comparison.png/pdf`)
- **描述**: Human 和 Rat 的覆盖率对比
- **关键发现**: 大多数模型在两个数据集上覆盖率一致
- **例外**: Pangolin 和 SpliceAI 在 rat 上覆盖率大幅下降

## 模型颜色方案（与 Human 图表保持一致）

| 模型 | 颜色代码 | 颜色名称 |
|------|---------|---------|
| AlphaGenome | `#1f77b4` | 蓝色 |
| Evo2 | `#ff7f0e` | **橙色** |
| Nucleotide Transformer | `#2ca02c` | 绿色 |
| Pangolin | `#d62728` | 红色 |
| SpliceAI | `#9467bd` | 紫色 |
| SpliceBERT | `#8c564b` | 棕色 |
| SpliceTransformer | `#e377c2` | 粉色 |
| DNABERT-2 (Logistic) | `#17becf` | 青色 |
| MMSplice (pathogenicity) | `#bcbd22` | 橄榄绿 |

## 关键性能指标总结

### Rat (RatGTEx) 最终结果

| 模型 | AUROC | AUPRC | Coverage |
|------|-------|-------|----------|
| AlphaGenome | **0.6005** | 0.5906 | 99.99% |
| SpliceAI | 0.5458 | **0.6496** | 38.44% |
| Nucleotide Transformer | 0.5063 | 0.5067 | 100% |
| Pangolin | 0.5057 | 0.6224 | 38.47% |
| SpliceBERT | 0.5049 | 0.5011 | 100% |
| SpliceTransformer | 0.5006 | 0.5014 | 100% |
| DNABERT2_Logistic | 0.5000 | 0.5000 | 100% |
| Evo2 | 0.4848 | 0.4942 | 100% |
| MMSplice | 0.4831 | 0.7275 | 5.05% |

## 使用说明

### 重新生成图表

```bash
# Rat 单独图表
cd /Users/byronsun/Desktop/AS_复现模型/BIB_review/code/ratGTEx/baselines
python plot_ratgtex_performance.py

# 跨物种对比图表
python plot_cross_species_comparison.py
```

### 依赖包

```bash
pip install pandas numpy matplotlib seaborn scikit-learn
```

## 论文引用

这些图表用于论文的以下章节：
- **4.2.2 Non-human Zero-shot Test (RatGTEx)** - Rat 单独图表
- **4.2.3 Human–Rat: Cross-species Gap Quantification** - 跨物种对比图表

## 脚本说明

1. `plot_ratgtex_performance.py` - 生成 Rat 单独性能图表
2. `plot_cross_species_comparison.py` - 生成跨物种对比图表

两个脚本都使用相同的颜色方案以确保视觉一致性。

