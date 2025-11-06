import pandas as pd
import json
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import (roc_curve, auc, precision_recall_curve, 
                             average_precision_score, classification_report, roc_auc_score)
import numpy as np

# --- 1. 定义最终的数据文件 ---
# 脚本和数据文件都在同一个文件夹下
final_results_file = 'AlphaGenome/Colab/alphagenome_COMPLETE_results_20250630_012057.json'

print(f"Starting final performance evaluation: {final_results_file}")

# --- 2. 加载最终合并的数据 ---
try:
    with open(final_results_file, 'r') as f:
        data = json.load(f)
    all_results = data['results']
except FileNotFoundError:
    print(f"ERROR: Final results JSON not found: {final_results_file}")
    exit()

# --- 新增步骤：加载并合并位置信息 ---
try:
    location_file = 'AlphaGenome/Colab/alphagenome_splicevardb_filtered.tsv'
    location_df = pd.read_csv(location_file, sep='\t', usecols=['hg38', 'location'])
    print(f"Loaded location file: {location_file}")
except FileNotFoundError:
    print(f"WARNING: Location TSV not found: {location_file}. Continuing without location-based analysis.")
    location_df = None

# --- 3. 预处理数据 ---
df = pd.json_normalize(all_results)

# Merge location info if available
if location_df is not None:
    if 'variant_id' in df.columns:
        df = pd.merge(df, location_df, left_on='variant_id', right_on='hg38', how='left')
        print("Location merged via 'variant_id'.")
    else:
        if 'input.hg38' in df.columns:
            df = pd.merge(df, location_df, left_on='input.hg38', right_on='hg38', how='left')
            print("Location merged via 'input.hg38'.")
        else:
            print("WARNING: No suitable ID column ('variant_id' or 'input.hg38') to merge location.")

# Filter and normalise
df = df[df['status'] == 'success'].copy()
df = df.dropna(subset=['classification', 'analysis.differential.max_abs_change']).copy()
df['classification'] = df['classification'].astype(str).str.strip().str.lower()
if 'location' in df.columns:
    df['location'] = df['location'].astype(str).str.strip().str.lower()

df['y_true'] = df['classification'].apply(lambda x: 1 if x in ['splice-altering', 'splicealtering'] else 0)
df['y_score'] = df['analysis.differential.max_abs_change']

y_true = df['y_true']
y_score = df['y_score']

# --- 4. Metrics ---
print("\n" + "="*50)
print("Core performance metrics")
print("="*50)

n_pos = int(y_true.sum())
n_neg = int(len(y_true) - n_pos)
print(f"Class counts -> Positive (Splice-altering): {n_pos}, Negative (Normal): {n_neg}")

both_classes = (n_pos > 0 and n_neg > 0)
if both_classes:
    fpr, tpr, roc_thresholds = roc_curve(y_true, y_score)
    roc_auc = auc(fpr, tpr)
    precision, recall, _ = precision_recall_curve(y_true, y_score)
    avg_precision = average_precision_score(y_true, y_score)
    print(f"Area Under ROC Curve (AUC): {roc_auc:.4f}")
    print(f"Average Precision (AP):   {avg_precision:.4f}")

    # Optimal threshold by Youden's J
    optimal_idx = np.argmax(tpr - fpr)
    optimal_threshold = roc_thresholds[optimal_idx]
    y_pred_optimal = (y_score >= optimal_threshold).astype(int)
    print(f"Optimal Threshold: {optimal_threshold:.4f}")
    print("\n" + "="*50)
    print("Classification report @ optimal threshold")
    print("="*50)
    print(classification_report(y_true, y_pred_optimal, target_names=['Normal', 'Splice-altering']))
else:
    print("WARNING: Only one class present in y_true. ROC/AUPRC are undefined.")
    print(f"y_score summary -> mean: {float(np.mean(y_score)):.4f}, std: {float(np.std(y_score)):.4f}, min: {float(np.min(y_score)):.4f}, max: {float(np.max(y_score)):.4f}")

# --- 6. Location-stratified performance ---
print("\n" + "="*50)
print("Location-stratified analysis")
print("="*50)

if 'location' in df.columns and both_classes:
    exonic_df = df[df['location'] == 'exonic']
    if not exonic_df.empty and exonic_df['y_true'].nunique() > 1:
        exonic_auc = roc_auc_score(exonic_df['y_true'], exonic_df['y_score'])
        print(f"Exonic AUROC: {exonic_auc:.4f} (n={len(exonic_df)})")
    else:
        print("Exonic AUROC: undefined (insufficient class diversity or no exonic records).")

    intronic_df = df[df['location'] == 'intronic']
    if not intronic_df.empty and intronic_df['y_true'].nunique() > 1:
        intronic_auc = roc_auc_score(intronic_df['y_true'], intronic_df['y_score'])
        print(f"Intronic AUROC: {intronic_auc:.4f} (n={len(intronic_df)})")
    else:
        print("Intronic AUROC: undefined (insufficient class diversity or no intronic records).")
else:
    print("Skipping location-based AUROC (missing 'location' or single-class dataset).")


# --- 7. 可视化 ---
plt.style.use('seaborn-v0_8-whitegrid')

if both_classes:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))
    fig.suptitle('AlphaGenome Final Performance Evaluation', fontsize=20)

    # ROC
    ax1.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (AUC = {roc_auc:.4f})')
    ax1.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    ax1.scatter(fpr[optimal_idx], tpr[optimal_idx], marker='o', color='red', s=100, zorder=5,
                label=f'Optimal Threshold ≈ {optimal_threshold:.2f}')
    ax1.set_title('Receiver Operating Characteristic (ROC) Curve', fontsize=16)
    ax1.set_xlabel('False Positive Rate', fontsize=12)
    ax1.set_ylabel('True Positive Rate', fontsize=12)
    ax1.legend(loc="lower right", fontsize=12)
    ax1.grid(True)

    # PR
    ax2.step(recall, precision, color='b', where='post', lw=2,
             label=f'Precision-Recall curve (AP = {avg_precision:.4f})')
    ax2.set_title('Precision-Recall Curve', fontsize=16)
    ax2.set_xlabel('Recall (Sensitivity)', fontsize=12)
    ax2.set_ylabel('Precision', fontsize=12)
    ax2.legend(loc="upper right", fontsize=12)
    ax2.grid(True)
    ax2.set_ylim([0.0, 1.05])
    ax2.set_xlim([0.0, 1.0])

    # Dataset info
    total_variants = len(df)
    n_splice_altering = int(y_true.sum())
    n_normal = int(total_variants - n_splice_altering)
    info_text = (
        f"Total Variants: {total_variants}\n"
        f"Splice-altering: {n_splice_altering}\n"
        f"Normal: {n_normal}"
    )
    fig.text(0.5, 0.02, info_text, ha='center', fontsize=12,
             bbox=dict(boxstyle='round,pad=0.5', fc='aliceblue', alpha=0.8))

    plt.tight_layout(rect=[0, 0.05, 1, 0.96])
    output_filename = 'AlphaGenome/Colab/final_performance_curves.png'
    plt.savefig(output_filename, dpi=300)
    print("\n" + "="*50)
    print(f"Saved performance curves to: {output_filename}")
    print("="*50)

# --- 8. 生成并保存 max_abs_change 分布图 ---
df['classification_label'] = df['y_true'].apply(lambda x: 'Splice-altering' if x == 1 else 'Normal')
plt.figure(figsize=(12, 8))
sns.violinplot(x='classification_label', y='y_score', data=df, order=['Normal', 'Splice-altering'], palette={'Normal': 'skyblue', 'Splice-altering': 'salmon'})
plt.title('Distribution of max_abs_change Score', fontsize=18)
plt.xlabel('Variant Classification', fontsize=14)
plt.ylabel('max_abs_change Score', fontsize=14)
plt.yscale('log')
plt.grid(True, which="both", ls="--")
plt.tight_layout()

distribution_plot_filename = 'AlphaGenome/Colab/max_abs_change_distribution.png'
plt.savefig(distribution_plot_filename, dpi=300)
print(f"Saved distribution figure to: {distribution_plot_filename}")
print("="*50)

plt.show()