import pandas as pd
import argparse
import os

def tsv_to_vcf(input_tsv, output_vcf):
    """
    Converts a specific TSV format to a VCF format file.

    The TSV is expected to have columns:
    'variant_id', 'chrom', 'pos_0based', 'ref', 'alt', 'label', 'location', 'sequence'

    The VCF format will be:
    #CHROM  POS ID  REF ALT QUAL    FILTER  INFO
    
    Notes:
    - Positions are 0-based in TSV but 1-based in VCF; add +1.
    - Remove 'chr' prefix for SpliceAI compatibility.
    """
    
    # Ensure the output directory exists
    output_dir = os.path.dirname(output_vcf)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        # Directory created

    # Read TSV with proper header
    df = pd.read_csv(input_tsv, sep='\t')
    print(f"Loaded {len(df)} variants from: {input_tsv}")
    
    # Check if required columns exist
    required_cols = ['chrom', 'pos_0based', 'ref', 'alt']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    # Basic VCF header
    vcf_header = [
        "##fileformat=VCFv4.2",
        "##source=tsv_to_vcf_converter",
        "##reference=hg38",
        "##contig=<ID=1,length=248956422>",
        "##contig=<ID=2,length=242193529>",
        "##contig=<ID=3,length=198295559>",
        "##contig=<ID=4,length=190214555>",
        "##contig=<ID=5,length=181538259>",
        "##contig=<ID=6,length=170805979>",
        "##contig=<ID=7,length=159345973>",
        "##contig=<ID=8,length=145138636>",
        "##contig=<ID=9,length=138394717>",
        "##contig=<ID=10,length=133797422>",
        "##contig=<ID=11,length=135086622>",
        "##contig=<ID=12,length=133275309>",
        "##contig=<ID=13,length=114364328>",
        "##contig=<ID=14,length=107043718>",
        "##contig=<ID=15,length=101991189>",
        "##contig=<ID=16,length=90338345>",
        "##contig=<ID=17,length=83257441>",
        "##contig=<ID=18,length=80373285>",
        "##contig=<ID=19,length=58617616>",
        "##contig=<ID=20,length=64444167>",
        "##contig=<ID=21,length=46709983>",
        "##contig=<ID=22,length=50818468>",
        "##contig=<ID=X,length=156040895>",
        "##contig=<ID=Y,length=57227415>",
        "##contig=<ID=MT,length=16569>",
        '##INFO=<ID=SpliceAI,Number=.,Type=String,Description="SpliceAIv1.3 annotations for SNVs and simple INDELs">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    ]

    with open(output_vcf, 'w') as f:
        # Write header
        for line in vcf_header:
            f.write(line + '\n')

        # Write data rows
        variants_written = 0
        for index, row in df.iterrows():
            # Extract data from TSV columns
            chrom = str(row['chrom'])
            pos_0based = int(row['pos_0based'])
            ref = str(row['ref'])
            alt = str(row['alt'])
            
            # Convert 0-based position to 1-based for VCF format
            pos_1based = pos_0based + 1
            
            # Ensure chromosome format consistency: remove 'chr' prefix
            if chrom.startswith('chr'):
                chrom_clean = chrom[3:]
            else:
                chrom_clean = chrom
            
            # VCF columns - using placeholders for non-essential fields
            vcf_id = '.'
            qual = '.'
            filt = 'PASS'
            info = '.'

            vcf_line = f"{chrom_clean}\t{pos_1based}\t{vcf_id}\t{ref}\t{alt}\t{qual}\t{filt}\t{info}\n"
            f.write(vcf_line)
            variants_written += 1
            
    print(f"Wrote VCF with {variants_written} variants -> {output_vcf}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert TSV file to VCF file for SpliceAI.")
    parser.add_argument('--input_tsv', required=True, help='Path to the input TSV file.')
    parser.add_argument('--output_vcf', required=True, help='Path to the output VCF file.')
    
    args = parser.parse_args()
    
    tsv_to_vcf(args.input_tsv, args.output_vcf) 