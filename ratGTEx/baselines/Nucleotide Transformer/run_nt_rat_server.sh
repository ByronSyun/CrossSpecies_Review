# ============================================================
# Step 1: 准备工作目录和数据
# ============================================================
cd /mnt/userdata4/splicing/NucleotideTransformer
mkdir -p ratgtex/results

# Step 2: 添加header到rat数据（一次性操作）
# ============================================================
cd ratgtex
echo -e "variant_id\tref_seq\talt_seq\tlabel\tchrom" > ratgtex_with_header.tsv
cat /mnt/userdata4/splicing/ratgetx/processed_data/ratgtex_silver_benchmark_balanced.tsv >> ratgtex_with_header.tsv

# 验证header已添加
head -1 ratgtex_with_header.tsv
# 应该输出: variant_id	ref_seq	alt_seq	label	chrom

# Step 3: 复制nt_score_variants.py（如果还没有）
# ============================================================
# 从human目录复制通用脚本
cp ../nt_score_variants.py .

# Step 4: 运行NT inference（无需手动激活环境）
# ============================================================
conda run -p /mnt/userdata4/splicing/conda_envs/nt-repo \
  python nt_score_variants.py \
    --model_id "InstaDeepAI/nucleotide-transformer-v2-500m-multi-species" \
    --input_tsv ratgtex_with_header.tsv \
    --output_prefix results/nt_rat \
    --window 8192 \
    --batch_size 16 \
    --num_workers 2

# 预计运行时间：2-4小时
# 输出：results/nt_rat_scores.csv