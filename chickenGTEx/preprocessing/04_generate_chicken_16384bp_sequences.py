#!/usr/bin/env python3
"""
Prepare 16384bp sequence data for AlphaGenome from ChickenGTEx 8192bp data
===========================================================================
AlphaGenome requires longer context (16384bp) for optimal performance.
This script extends the existing 8192bp sequences to 16384bp by fetching
additional flanking sequence from the chicken reference genome.

Input: chickengtex_silver_benchmark_balanced.tsv (8192bp, no header)
Output: chickengtex_silver_benchmark_balanced_len16384.tsv (16384bp, no header)
"""

import logging
import re
from typing import Tuple, Optional
import pandas as pd
from pyfaidx import Fasta

logging.basicConfig(level=logging.INFO, format='%(message)s')

# =============================================================================
# Configuration
# =============================================================================
INPUT_TSV = "/mnt/userdata4/splicing/ChickenGTEx/processed_data/chickengtex_silver_benchmark_balanced.tsv"
OUTPUT_TSV = "/mnt/userdata4/splicing/ChickenGTEx/processed_data/chickengtex_silver_benchmark_balanced_len16384.tsv"
REFERENCE_FASTA = "/mnt/userdata4/splicing/ChickenGTEx/reference_genome/Gallus_gallus.GRCg6a.dna.toplevel.fa"

ORIGINAL_LENGTH = 8192
TARGET_LENGTH = 16384
HALF_EXTENSION = (TARGET_LENGTH - ORIGINAL_LENGTH) // 2  # 4096bp on each side

# =============================================================================
# Helper Functions
# =============================================================================
def normalise_chrom(chrom: str) -> Tuple[str, str]:
    """Return (primary, fallback) chrom names after normalisation.
    - Strip leading 'chr'
    - Map M/chrM -> MT
    - Provide a fallback with 'chr' prefix for UCSC-style FASTA if needed
    """
    c = chrom
    if c.lower().startswith('chr'):
        c = c[3:]
    if c in ('M', 'Mt', 'm', 'MT'):
        c = 'MT'
    primary = c
    fallback = 'chr' + c
    return primary, fallback

def parse_variant_id(variant_id: str) -> Optional[Tuple[str, int, str, str]]:
    """Parse variant_id to extract chrom, pos, ref, alt.
    Supports multiple formats:
    - chrom:pos_ref/alt (e.g., 1:15450992_C/T)
    - chrom_pos_ref_alt (e.g., 1_59848360_T_A)
    """
    # Try format 1: chrom:pos_ref/alt
    match = re.match(r'^(.+):(\d+)_([ACGTacgt])/([ACGTacgt])$', variant_id)
    if match:
        return match.group(1), int(match.group(2)), match.group(3).upper(), match.group(4).upper()
    
    # Try format 2: chrom_pos_ref_alt
    match = re.match(r'^(.+?)_(\d+)_([ACGTacgt])_([ACGTacgt])$', variant_id)
    if match:
        return match.group(1), int(match.group(2)), match.group(3).upper(), match.group(4).upper()
    
    return None

def fetch_slice(fasta: Fasta, chrom: str, start: int, end: int) -> Optional[str]:
    """Fetch sequence from FASTA with error handling."""
    try:
        return str(fasta[str(chrom)][start - 1:end]).upper()
    except (KeyError, IndexError):
        return None

def extract_sequence_robust(fasta: Fasta, chrom: str, pos: int, ref: str, alt: str, seq_len: int) -> Tuple[Optional[str], Optional[str]]:
    """Extract reference and alternative sequences with chromosome name normalization.
    
    Args:
        fasta: Pyfaidx Fasta object
        chrom: Chromosome name
        pos: 1-based position
        ref: Reference allele
        alt: Alternative allele
        seq_len: Total sequence length (e.g., 16384)
    
    Returns:
        (ref_seq, alt_seq) or (None, None) if failed
    """
    start = pos - (seq_len // 2) + 1
    end = pos + (seq_len // 2)
    
    if start < 1:
        return None, None
    
    # Try primary chromosome name
    primary, fallback = normalise_chrom(chrom)
    ref_seq = fetch_slice(fasta, primary, start, end)
    
    # If failed, try fallback name
    if ref_seq is None:
        ref_seq = fetch_slice(fasta, fallback, start, end)
    
    if ref_seq is None:
        return None, None
    
    # Verify sequence length
    if len(ref_seq) != seq_len:
        return None, None
    
    # Verify reference allele at center position
    center_idx = (seq_len // 2) - 1
    if ref_seq[center_idx] != ref:
        return None, None
    
    # Create alternative sequence
    alt_seq = ref_seq[:center_idx] + alt + ref_seq[center_idx + 1:]
    
    return ref_seq, alt_seq

# =============================================================================
# Main Processing
# =============================================================================
def extend_sequences():
    """Extend 8192bp sequences to 16384bp by fetching from reference genome."""
    
    logging.info("="*80)
    logging.info("AlphaGenome Data Preparation: 8192bp → 16384bp Extension")
    logging.info("="*80)
    
    # Load input data
    logging.info(f"Loading input data: {INPUT_TSV}")
    df = pd.read_csv(INPUT_TSV, sep='\t', header=None, 
                     names=['variant_id', 'ref_seq', 'alt_seq', 'label', 'tissue_id'])
    logging.info(f"  Loaded {len(df):,} variants")
    
    # Load reference genome
    logging.info(f"Loading reference genome: {REFERENCE_FASTA}")
    fasta = Fasta(REFERENCE_FASTA, key_function=lambda k: k.split(' ')[0])
    logging.info(f"  Reference genome loaded with {len(fasta.keys())} chromosomes")
    
    # Parse variant coordinates and extract sequences
    logging.info(f"Parsing variants and extracting {TARGET_LENGTH}bp sequences...")
    
    rebuilt_rows = []
    skipped_parse = 0
    skipped_extraction = 0
    
    for idx, row in df.iterrows():
        variant_id = str(row['variant_id'])
        parsed = parse_variant_id(variant_id)
        
        if not parsed:
            skipped_parse += 1
            continue
        
        chrom, pos, ref, alt = parsed
        
        # Extract sequences with robust chromosome name handling
        ref_seq, alt_seq = extract_sequence_robust(fasta, chrom, pos, ref, alt, TARGET_LENGTH)
        
        if ref_seq is None:
            skipped_extraction += 1
            continue
        
        rebuilt_rows.append([
            variant_id,
            ref_seq,
            alt_seq,
            int(row['label']),
            int(row['tissue_id']) if pd.notnull(row['tissue_id']) else 15
        ])
    
    if not rebuilt_rows:
        logging.error("No variants were successfully extended!")
        logging.error(f"  Skipped (parse failed): {skipped_parse:,}")
        logging.error(f"  Skipped (extraction failed): {skipped_extraction:,}")
        return
    
    # Create output dataframe
    output_df = pd.DataFrame(rebuilt_rows, columns=['variant_id', 'ref_sequence', 'alt_sequence', 'label', 'tissue_id'])
    
    logging.info("="*80)
    logging.info("✅ Data preparation completed!")
    logging.info(f"  Input variants: {len(df):,}")
    logging.info(f"  Successfully extended: {len(output_df):,}")
    logging.info(f"  Skipped (parse failed): {skipped_parse:,}")
    logging.info(f"  Skipped (extraction failed): {skipped_extraction:,}")
    logging.info(f"  Success rate: {len(output_df)/len(df)*100:.1f}%")
    logging.info("="*80)
    
    # Save output (same format as input: no header, 5 columns)
    logging.info(f"Saving to: {OUTPUT_TSV}")
    output_df.to_csv(OUTPUT_TSV, sep='\t', header=False, index=False)
    logging.info(f"Output format: variant_id | ref_sequence | alt_sequence | label | tissue_id")
    logging.info(f"Sequence length: {TARGET_LENGTH}bp")
    logging.info("="*80)

if __name__ == "__main__":
    extend_sequences()

