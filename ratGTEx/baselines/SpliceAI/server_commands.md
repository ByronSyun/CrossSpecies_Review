# 激活SpliceAI环境
conda activate /mnt/userdata4/splicing/conda_envs/spliceai_env

# 进入工作目录
cd /mnt/userdata4/splicing/SpliceAI

# 运行TSV转VCF转换
python scripts/tsv_to_vcf_rat.py \
    --input_tsv /mnt/userdata4/splicing/ratgetx/processed_data/ratgtex_silver_benchmark_balanced.tsv \
    --output_vcf rat_data/spliceai_input_rat.vcf


python create_rat_annotation.py \
    --input_gtf /mnt/userdata4/splicing/SpliceAI/reference_genome/Rattus_norvegicus.Rnor_6.0.84.final.gtf \
    --output_txt rat_data/rat_rn6_annotation.txt

# 运行SpliceAI预测
spliceai \
    -I rat_data/spliceai_input_rat.vcf \
    -O rat_data/spliceai_predictions_rat.vcf \
    -R /mnt/userdata4/splicing/SpliceAI/reference_genome/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa \
    -A rat_data/rat_rn6_annotation.txt \
    -D 4999


# 运行解析脚本
python parse_spliceai_output.py \
    -I rat_data/spliceai_predictions_rat.vcf \
    -O rat_data/spliceai_parsed_scores_rat.tsv

# 运行评估脚本（使用正确的参数名）
python evaluate_spliceai_rat.py \
    --ground_truth /mnt/userdata4/splicing/ratgetx/processed_data/ratgtex_silver_benchmark_balanced.tsv \
    --predictions rat_data/spliceai_parsed_scores_rat.tsv \
    --output rat_data/results/spliceai_rat_performance.csv