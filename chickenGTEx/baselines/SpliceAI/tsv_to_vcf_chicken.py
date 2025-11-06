import pandas as pd
import argparse
import os

def extract_variant_info(variant_id):
    """
    Extract chromosome, position, ref, alt from variant_id format like '1:15450992_C/T'
    Returns: (chrom, pos_1based, ref, alt)
    """
    if ':' in variant_id and '_' in variant_id and '/' in variant_id:
        # Format: '1:15450992_C/T'
        chrom_pos, ref_alt = variant_id.split('_')
        chrom, pos_str = chrom_pos.split(':')
        ref, alt = ref_alt.split('/')
        pos_1based = int(pos_str)
        return chrom, pos_1based, ref, alt
    else:
        raise ValueError(f"Unexpected variant_id format: {variant_id}")

def tsv_to_vcf_chicken(input_tsv, output_vcf):
    """
    Converts chicken TSV format to VCF format for SpliceAI.
    
    The chicken TSV has format: variant_id, ref_seq, alt_seq, label, tissue_id
    where variant_id is like '1:15450992_C/T'
    """
    
    # Ensure output directory exists
    output_dir = os.path.dirname(output_vcf)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Read TSV - no header, so we need to specify column names
    df = pd.read_csv(input_tsv, sep='\t', header=None, 
                     names=['variant_id', 'ref_seq', 'alt_seq', 'label', 'tissue_id'])
    print(f"Loaded {len(df)} chicken variants from: {input_tsv}")
    
    # Chicken chromosome contigs for VCF header (GRCg6a assembly)
    # Including microchromosomes 29-38 to avoid "contig not defined" errors
    vcf_header = [
        "##fileformat=VCFv4.2",
        "##source=tsv_to_vcf_chicken_converter",
        "##reference=GRCg6a",
        "##contig=<ID=1,length=195276750>",
        "##contig=<ID=2,length=149736717>",
        "##contig=<ID=3,length=106816920>",
        "##contig=<ID=4,length=91407136>",
        "##contig=<ID=5,length=59100713>",
        "##contig=<ID=6,length=33212725>",
        "##contig=<ID=7,length=36332221>",
        "##contig=<ID=8,length=30281471>",
        "##contig=<ID=9,length=25426900>",
        "##contig=<ID=10,length=21075988>",
        "##contig=<ID=11,length=19937869>",
        "##contig=<ID=12,length=21118477>",
        "##contig=<ID=13,length=17944816>",
        "##contig=<ID=14,length=16301643>",
        "##contig=<ID=15,length=12596774>",
        "##contig=<ID=16,length=7865194>",
        "##contig=<ID=17,length=11193993>",
        "##contig=<ID=18,length=10499472>",
        "##contig=<ID=19,length=9809832>",
        "##contig=<ID=20,length=14260733>",
        "##contig=<ID=21,length=6936039>",
        "##contig=<ID=22,length=4382547>",
        "##contig=<ID=23,length=6232768>",
        "##contig=<ID=24,length=6326522>",
        "##contig=<ID=25,length=2321640>",
        "##contig=<ID=26,length=5319930>",
        "##contig=<ID=27,length=4929992>",
        "##contig=<ID=28,length=4803398>",
        "##contig=<ID=29,length=1000000>",
        "##contig=<ID=30,length=1000000>",
        "##contig=<ID=31,length=1000000>",
        "##contig=<ID=32,length=1000000>",
        "##contig=<ID=33,length=1000000>",
        "##contig=<ID=34,length=1000000>",
        "##contig=<ID=35,length=1000000>",
        "##contig=<ID=36,length=1000000>",
        "##contig=<ID=37,length=1000000>",
        "##contig=<ID=38,length=1000000>",
        "##contig=<ID=W,length=6743233>",
        "##contig=<ID=Z,length=81039894>",
        "##contig=<ID=MT,length=16775>",
        '##INFO=<ID=SpliceAI,Number=.,Type=String,Description="SpliceAIv1.3 annotations for SNVs and simple INDELs">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    ]

    with open(output_vcf, 'w') as f:
        # Write header
        for line in vcf_header:
            f.write(line + '\n')

        # Write data rows
        variants_written = 0
        skipped_variants = 0
        
        for index, row in df.iterrows():
            try:
                # Extract variant information
                chrom, pos_1based, ref, alt = extract_variant_info(row['variant_id'])
                
                # VCF columns - using placeholders for non-essential fields
                vcf_id = '.'
                qual = '.'
                filt = 'PASS'
                info = '.'

                vcf_line = f"{chrom}\t{pos_1based}\t{vcf_id}\t{ref}\t{alt}\t{qual}\t{filt}\t{info}\n"
                f.write(vcf_line)
                variants_written += 1
                
            except Exception as e:
                print(f"Warning: Skipping variant {row['variant_id']}: {e}")
                skipped_variants += 1
                continue
            
    print(f"Wrote VCF with {variants_written} variants -> {output_vcf}")
    if skipped_variants > 0:
        print(f"Skipped {skipped_variants} variants due to format issues")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert chicken TSV file to VCF file for SpliceAI.")
    parser.add_argument('--input_tsv', required=True, help='Path to the input chicken TSV file.')
    parser.add_argument('--output_vcf', required=True, help='Path to the output VCF file.')
    
    args = parser.parse_args()
    
    tsv_to_vcf_chicken(args.input_tsv, args.output_vcf)

