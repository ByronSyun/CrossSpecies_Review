import pandas as pd
import argparse
import os

def extract_variant_info(variant_id):
    """
    Extract chromosome, position, ref, alt from variant_id format like '13:65362642_C/T'
    Returns: (chrom, pos_1based, ref, alt)
    """
    if ':' in variant_id and '_' in variant_id and '/' in variant_id:
        # Format: '13:65362642_C/T'
        chrom_pos, ref_alt = variant_id.split('_')
        chrom, pos_str = chrom_pos.split(':')
        ref, alt = ref_alt.split('/')
        pos_1based = int(pos_str)
        return chrom, pos_1based, ref, alt
    else:
        raise ValueError(f"Unexpected variant_id format: {variant_id}")

def tsv_to_vcf_pig(input_tsv, output_vcf):
    """
    Converts pig TSV format to VCF format for SpliceAI.
    
    The pig TSV has format: variant_id, ref_seq, alt_seq, label, tissue_id
    where variant_id is like '13:65362642_C/T'
    """
    
    # Ensure output directory exists
    output_dir = os.path.dirname(output_vcf)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Read TSV - no header, so we need to specify column names
    df = pd.read_csv(input_tsv, sep='\t', header=None, 
                     names=['variant_id', 'ref_seq', 'alt_seq', 'label', 'tissue_id'])
    print(f"Loaded {len(df)} pig variants from: {input_tsv}")
    
    # Pig chromosome contigs for VCF header (Sscrofa11.1 assembly)
    vcf_header = [
        "##fileformat=VCFv4.2",
        "##source=tsv_to_vcf_pig_converter",
        "##reference=Sscrofa11.1",
        "##contig=<ID=1,length=274330532>",
        "##contig=<ID=2,length=151935994>",
        "##contig=<ID=3,length=132848913>",
        "##contig=<ID=4,length=130910915>",
        "##contig=<ID=5,length=104526007>",
        "##contig=<ID=6,length=170843587>",
        "##contig=<ID=7,length=121844099>",
        "##contig=<ID=8,length=138966237>",
        "##contig=<ID=9,length=139009144>",
        "##contig=<ID=10,length=69359453>",
        "##contig=<ID=11,length=79169978>",
        "##contig=<ID=12,length=61602749>",
        "##contig=<ID=13,length=208334590>",
        "##contig=<ID=14,length=141755446>",
        "##contig=<ID=15,length=140412725>",
        "##contig=<ID=16,length=79944280>",
        "##contig=<ID=17,length=63494081>",
        "##contig=<ID=18,length=55982971>",
        "##contig=<ID=X,length=125939595>",
        "##contig=<ID=Y,length=43547828>",
        "##contig=<ID=MT,length=16613>",
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
    parser = argparse.ArgumentParser(description="Convert pig TSV file to VCF file for SpliceAI.")
    parser.add_argument('--input_tsv', required=True, help='Path to the input pig TSV file.')
    parser.add_argument('--output_vcf', required=True, help='Path to the output VCF file.')
    
    args = parser.parse_args()
    
    tsv_to_vcf_pig(args.input_tsv, args.output_vcf)

