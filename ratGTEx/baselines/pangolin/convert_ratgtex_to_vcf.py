#!/usr/bin/env python3
"""
Convert RatGTEx TSV data to VCF format for Pangolin input.

This script converts the RatGTEx silver-standard benchmark TSV file 
to VCF format that Pangolin can process directly.

Input format (TSV):
- variant_id (e.g., 15:55290237_A/G)
- ref_sequence (8192bp)
- alt_sequence (8192bp) 
- label (0/1)
- tissue_id

Output format (VCF):
- Standard VCF with CHROM, POS, REF, ALT
- Saves mapping file with labels for evaluation
"""

import pandas as pd
import os
import argparse
import logging
import re
from datetime import datetime

logging.basicConfig(level=logging.INFO, format='%(asctime)s - [%(levelname)s] - %(message)s')

def parse_variant_id(variant_id):
    """
    Parse variant_id to extract CHROM, POS, REF, ALT.
    
    Args:
        variant_id (str): Format like "15:55290237_A/G" or "chr15:55290237_A/G"
        
    Returns:
        dict: {chrom, pos, ref, alt} or None if parsing fails
    """
    # Handle both chr15:pos_ref/alt and 15:pos_ref/alt formats
    pattern = r'^(chr)?(\w+):(\d+)_([ACGT]+)/([ACGT]+)$'
    match = re.match(pattern, variant_id)
    
    if not match:
        logging.warning(f"Could not parse variant_id: {variant_id}")
        return None
    
    chr_prefix, chrom, pos_str, ref, alt = match.groups()
    
    # For rat genome (rn6), remove 'chr' prefix to match FASTA/GTF naming (1, 2, 3, etc.)
    if chr_prefix:
        chrom = chrom  # Remove chr prefix, just use the number/letter part
    # chrom is already without prefix if chr_prefix is None
        
    pos = int(pos_str)
    
    return {
        'CHROM': chrom,
        'POS': pos,
        'REF': ref,
        'ALT': alt
    }

def convert_ratgtex_to_vcf(input_tsv, output_vcf, sequence_length=8192):
    """
    Convert RatGTEx TSV to VCF format for Pangolin.
    
    Args:
        input_tsv (str): Path to input TSV file
        output_vcf (str): Path to output VCF file
        sequence_length (int): Expected sequence length (should be 8192)
    """
    logging.info("=== RatGTEx to Pangolin VCF Conversion ===")
    logging.info(f"Input: {input_tsv}")
    logging.info(f"Output: {output_vcf}")
    logging.info(f"Expected sequence length: {sequence_length}bp")
    
    # Check input file
    if not os.path.exists(input_tsv):
        raise FileNotFoundError(f"Input file not found: {input_tsv}")
    
    # Read TSV file
    logging.info("Loading RatGTEx data...")
    df = pd.read_csv(input_tsv, sep='\t', header=None)
    df.columns = ['variant_id', 'ref_sequence', 'alt_sequence', 'label', 'tissue_id']
    
    total_variants = len(df)
    logging.info(f"Loaded {total_variants:,} variants")
    
    # Validate sequence lengths
    ref_lengths = df['ref_sequence'].str.len()
    alt_lengths = df['alt_sequence'].str.len()
    
    if not all(ref_lengths == sequence_length):
        logging.warning(f"Some reference sequences are not {sequence_length}bp")
    if not all(alt_lengths == sequence_length):
        logging.warning(f"Some alternative sequences are not {sequence_length}bp")
    
    # Parse variant IDs
    logging.info("Parsing variant coordinates...")
    parsed_variants = []
    failed_count = 0
    
    for idx, row in df.iterrows():
        parsed = parse_variant_id(row['variant_id'])
        if parsed:
            parsed['variant_id'] = row['variant_id']
            parsed['label'] = int(row['label'])
            parsed['tissue_id'] = row['tissue_id']
            parsed_variants.append(parsed)
        else:
            failed_count += 1
            if failed_count <= 5:  # Show first 5 failures
                logging.warning(f"Failed to parse: {row['variant_id']}")
    
    if failed_count > 0:
        logging.warning(f"Failed to parse {failed_count} variant IDs")
    
    logging.info(f"Successfully parsed {len(parsed_variants):,} variants")
    
    # Convert to DataFrame
    vcf_df = pd.DataFrame(parsed_variants)
    
    # Create VCF content
    logging.info("Generating VCF format...")
    
    vcf_lines = []
    
    # VCF header
    vcf_lines.append("##fileformat=VCFv4.2")
    vcf_lines.append(f"##fileDate={datetime.now().strftime('%Y%m%d')}")
    vcf_lines.append("##reference=ratNor6")
    vcf_lines.append("##source=RatGTEx_Silver_Standard")
    vcf_lines.append("##INFO=<ID=LABEL,Number=1,Type=Integer,Description=\"Silver standard label (0=negative, 1=positive)\">")
    vcf_lines.append("##INFO=<ID=TISSUE,Number=1,Type=String,Description=\"Tissue ID from RatGTEx\">")
    
    # Add contig information
    unique_chroms = sorted(vcf_df['CHROM'].unique())
    for chrom in unique_chroms:
        vcf_lines.append(f"##contig=<ID={chrom}>")
    
    # Column header
    vcf_lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    
    # Data lines
    for _, row in vcf_df.iterrows():
        info_field = f"LABEL={row['label']};TISSUE={row['tissue_id']}"
        vcf_line = f"{row['CHROM']}\t{row['POS']}\t{row['variant_id']}\t{row['REF']}\t{row['ALT']}\t.\tPASS\t{info_field}"
        vcf_lines.append(vcf_line)
    
    # Write VCF file
    os.makedirs(os.path.dirname(output_vcf), exist_ok=True)
    with open(output_vcf, 'w') as f:
        f.write('\n'.join(vcf_lines))
    
    logging.info(f"VCF written to: {output_vcf}")
    
    # Also save a mapping file for evaluation
    mapping_file = output_vcf.replace('.vcf', '_mapping.csv')
    vcf_df.to_csv(mapping_file, index=False)
    logging.info(f"Mapping file saved to: {mapping_file}")
    
    # Summary statistics
    logging.info("\n=== Conversion Summary ===")
    logging.info(f"Total variants processed: {len(vcf_df):,}")
    logging.info(f"Positive samples (label=1): {vcf_df['label'].sum():,}")
    logging.info(f"Negative samples (label=0): {(len(vcf_df) - vcf_df['label'].sum()):,}")
    logging.info(f"Unique chromosomes: {len(unique_chroms)}")
    logging.info(f"Chromosome range: {unique_chroms}")
    logging.info("============================\n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Convert RatGTEx TSV to VCF format for Pangolin input"
    )
    
    # Default paths (server only)
    default_input = "/mnt/userdata4/splicing/ratgetx/processed_data/ratgtex_silver_benchmark_balanced.tsv"
    default_output = "/mnt/userdata4/splicing/pangolin/data/ratgtex_for_pangolin.vcf"
    
    parser.add_argument(
        '--input',
        type=str,
        default=default_input,
        help=f'Path to input TSV file (default: {default_input})'
    )
    
    parser.add_argument(
        '--output',
        type=str,
        default=default_output,
        help=f'Path to output VCF file (default: {default_output})'
    )
    
    parser.add_argument(
        '--sequence_length',
        type=int,
        default=8192,
        help='Expected sequence length (default: 8192)'
    )
    
    args = parser.parse_args()
    
    try:
        convert_ratgtex_to_vcf(args.input, args.output, args.sequence_length)
        logging.info("✅ Conversion completed successfully!")
    except Exception as e:
        logging.error(f"❌ Conversion failed: {e}")
        raise
