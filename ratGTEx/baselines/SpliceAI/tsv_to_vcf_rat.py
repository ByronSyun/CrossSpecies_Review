import pandas as pd
import argparse
import os

def extract_variant_info(variant_id):
    """
    Extract chromosome, position, ref, alt from variant_id format like '15:55290237_A/G'
    Returns: (chrom, pos_1based, ref, alt)
    """
    if ':' in variant_id and '_' in variant_id and '/' in variant_id:
        # Format: '15:55290237_A/G'
        chrom_pos, ref_alt = variant_id.split('_')
        chrom, pos_str = chrom_pos.split(':')
        ref, alt = ref_alt.split('/')
        pos_1based = int(pos_str)
        return chrom, pos_1based, ref, alt
    else:
        raise ValueError(f"Unexpected variant_id format: {variant_id}")

def tsv_to_vcf_rat(input_tsv, output_vcf):
    """
    Converts rat TSV format to VCF format for SpliceAI.
    
    The rat TSV has format: variant_id, ref_seq, alt_seq, label, tissue_id
    where variant_id is like '15:55290237_A/G'
    """
    
    # Ensure output directory exists
    output_dir = os.path.dirname(output_vcf)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Read TSV - no header, so we need to specify column names
    df = pd.read_csv(input_tsv, sep='\t', header=None, 
                     names=['variant_id', 'ref_seq', 'alt_seq', 'label', 'tissue_id'])
    print(f"Loaded {len(df)} rat variants from: {input_tsv}")
    
    # Rat chromosome contigs for VCF header (rn6 assembly)
    vcf_header = [
        "##fileformat=VCFv4.2",
        "##source=tsv_to_vcf_rat_converter",
        "##reference=rn6",
        "##contig=<ID=1,length=282763074>",
        "##contig=<ID=2,length=266435125>", 
        "##contig=<ID=3,length=177699992>",
        "##contig=<ID=4,length=184226339>",
        "##contig=<ID=5,length=173707219>",
        "##contig=<ID=6,length=147991367>",
        "##contig=<ID=7,length=145729302>",
        "##contig=<ID=8,length=133307652>",
        "##contig=<ID=9,length=122095297>",
        "##contig=<ID=10,length=112626471>",
        "##contig=<ID=11,length=90463843>",
        "##contig=<ID=12,length=52716770>",
        "##contig=<ID=13,length=114033958>",
        "##contig=<ID=14,length=115493446>",
        "##contig=<ID=15,length=111246239>",
        "##contig=<ID=16,length=90668790>",
        "##contig=<ID=17,length=90843779>",
        "##contig=<ID=18,length=88201929>",
        "##contig=<ID=19,length=62275575>",
        "##contig=<ID=20,length=56205956>",
        "##contig=<ID=X,length=159970021>",
        "##contig=<ID=Y,length=3310458>",
        "##contig=<ID=MT,length=16313>",
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
    parser = argparse.ArgumentParser(description="Convert rat TSV file to VCF file for SpliceAI.")
    parser.add_argument('--input_tsv', required=True, help='Path to the input rat TSV file.')
    parser.add_argument('--output_vcf', required=True, help='Path to the output VCF file.')
    
    args = parser.parse_args()
    
    tsv_to_vcf_rat(args.input_tsv, args.output_vcf)
