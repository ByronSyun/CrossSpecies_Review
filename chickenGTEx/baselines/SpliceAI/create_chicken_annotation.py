#!/usr/bin/env python3
"""
Create SpliceAI-compatible annotation file from chicken GTF.
Converts GTF format to the tab-separated format required by SpliceAI:
#NAME    CHROM    STRAND    TX_START    TX_END    EXON_START    EXON_END
"""

import argparse
import gffutils
import pandas as pd
from collections import defaultdict
import os

def gtf_to_spliceai_annotation(gtf_path, output_path):
    """
    Convert chicken GTF to SpliceAI annotation format.
    """
    print(f"Reading GTF file: {gtf_path}")
    
    # Create a temporary database from GTF
    db_path = gtf_path + ".db"
    if os.path.exists(db_path):
        print(f"Using existing database: {db_path}")
        db = gffutils.FeatureDB(db_path)
    else:
        print(f"Creating database from GTF (this may take a while)...")
        db = gffutils.create_db(gtf_path, db_path, force=True, keep_order=True,
                               merge_strategy='merge', sort_attribute_values=True)
    
    print("Processing transcripts...")
    
    # Collect transcript information
    transcript_data = []
    processed_count = 0
    
    for transcript in db.features_of_type('transcript'):
        processed_count += 1
        if processed_count % 1000 == 0:
            print(f"Processed {processed_count} transcripts...")
        
        # Get transcript info
        chrom = transcript.chrom
        if chrom.startswith('chr'):
            chrom = chrom[3:]  # Remove 'chr' prefix for SpliceAI compatibility
        
        strand = transcript.strand
        tx_start = transcript.start - 1  # Convert to 0-based for SpliceAI
        tx_end = transcript.end - 1      # Convert to 0-based for SpliceAI
        
        # Get gene name/ID
        gene_name = transcript.attributes.get('gene_name', transcript.attributes.get('gene_id', [transcript.id]))[0]
        
        # Get exons for this transcript
        exons = list(db.children(transcript, featuretype='exon', order_by='start'))
        
        if not exons:
            continue  # Skip transcripts without exons
        
        # Collect exon coordinates (0-based for SpliceAI)
        exon_starts = [str(exon.start - 1) for exon in exons]
        exon_ends = [str(exon.end - 1) for exon in exons]
        
        transcript_data.append({
            '#NAME': gene_name,
            'CHROM': chrom,
            'STRAND': strand,
            'TX_START': tx_start,
            'TX_END': tx_end,
            'EXON_START': ','.join(exon_starts) + ',',
            'EXON_END': ','.join(exon_ends) + ','
        })
    
    print(f"Processed {processed_count} transcripts total")
    print(f"Creating annotation file with {len(transcript_data)} entries...")
    
    # Create DataFrame and save
    df = pd.DataFrame(transcript_data)
    
    # Sort by chromosome and position
    def sort_key(chrom):
        if chrom.isdigit():
            return (0, int(chrom))
        elif chrom in ['W', 'Z', 'MT']:
            return (1, {'W': 1, 'Z': 2, 'MT': 3}[chrom])
        else:
            return (2, chrom)
    
    df['sort_key'] = df['CHROM'].apply(sort_key)
    df = df.sort_values(['sort_key', 'TX_START']).drop('sort_key', axis=1)
    
    # Save to file
    df.to_csv(output_path, sep='\t', index=False)
    print(f"Saved SpliceAI annotation file: {output_path}")
    print(f"Total entries: {len(df)}")
    
    # Show sample entries
    print("\nSample entries:")
    print(df.head())
    
    return output_path

def main():
    parser = argparse.ArgumentParser(description="Convert chicken GTF to SpliceAI annotation format")
    parser.add_argument('--input_gtf', required=True, help='Path to chicken GTF file')
    parser.add_argument('--output_txt', required=True, help='Path to output SpliceAI annotation file')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input_gtf):
        print(f"Error: GTF file not found: {args.input_gtf}")
        return 1
    
    try:
        gtf_to_spliceai_annotation(args.input_gtf, args.output_txt)
        print("Conversion completed successfully!")
        return 0
    except Exception as e:
        print(f"Error during conversion: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    exit(main())

