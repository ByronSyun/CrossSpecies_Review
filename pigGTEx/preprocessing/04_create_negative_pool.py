import pandas as pd
import logging
import os
import gzip
from pyfaidx import Fasta
from tqdm import tqdm

logging.basicConfig(level=logging.INFO, format='%(asctime)s - [%(levelname)s] - %(message)s')

SEQ_LENGTH = 8192
UNKNOWN_TISSUE_ID = 15

def load_chromosome_mapping():
    """Load chromosome mapping from simple names to NCBI RefSeq IDs."""
    mapping_file = './processed_data/chromosome_mapping.tsv'
    if not os.path.exists(mapping_file):
        logging.error(f"Chromosome mapping file not found: {mapping_file}")
        return None
    
    mapping_df = pd.read_csv(mapping_file, sep='\t')
    return dict(zip(mapping_df['simple_chrom'].astype(str), mapping_df['ncbi_chrom']))

def get_sequence_context(chrom, pos, fasta_handle, chrom_mapping):
    """Extract 8192bp sequence context around a given position."""
    simple_chrom = str(chrom)
    if simple_chrom not in chrom_mapping:
        return None
    
    ncbi_chrom = chrom_mapping[simple_chrom]
    start = pos - (SEQ_LENGTH // 2)
    end = pos + (SEQ_LENGTH // 2)
    
    try:
        if start < 1 or end > len(fasta_handle[ncbi_chrom]):
            return None
            
        return fasta_handle[ncbi_chrom][start - 1:end].seq.upper()
    except (KeyError, IndexError):
        return None

def build_negative_sequences(ref_allele, alt_allele, context_sequence):
    """Build ref/alt sequences for a negative sample."""
    variant_offset = SEQ_LENGTH // 2
    
    ref_sequence = context_sequence
    alt_sequence = (
        context_sequence[:variant_offset - 1] +
        alt_allele +
        context_sequence[variant_offset - 1 + len(ref_allele):]
    )
    
    final_ref_seq = (ref_sequence + 'N' * SEQ_LENGTH)[:SEQ_LENGTH]
    final_alt_seq = (alt_sequence + 'N' * SEQ_LENGTH)[:SEQ_LENGTH]

    return final_ref_seq, final_alt_seq

def create_negative_pool():
    """
    Create high-quality negative pool by excluding all positive sQTL locations.
    """
    positive_samples_file = './processed_data/pig_positive_samples_with_tissue.tsv'
    vcf_file = './pig_raw_data/sus_scrofa.vcf.gz'
    fasta_file = './reference_genome/Sscrofa11.1_genomic.fna'
    output_file = './processed_data/pig_negative_pool_HIGH_QUALITY.tsv'
    
    MAX_NEGATIVE_SAMPLES = 225000
    
    logging.info("Creating High-Quality Negative Pool")
    logging.info(f"Target: {MAX_NEGATIVE_SAMPLES:,} high-quality negative samples")

    logging.info("Loading chromosome mapping")
    chrom_mapping = load_chromosome_mapping()
    if not chrom_mapping:
        return
    logging.info(f"Loaded mapping for {len(chrom_mapping)} chromosomes")

    logging.info("Loading positive sample locations for exclusion")
    
    try:
        positive_df = pd.read_csv(positive_samples_file, sep='\t')
        positive_locations = set()
        
        for _, row in positive_df.iterrows():
            try:
                chrom = str(row['CHR'])
                pos = int(row['POS'])
                positive_locations.add((chrom, pos))
            except (ValueError, IndexError, KeyError):
                continue
                
        logging.info(f"Loaded {len(positive_locations):,} unique positive locations to exclude")
    except FileNotFoundError:
        logging.error(f"Positive samples file not found: {positive_samples_file}")
        logging.error("Please run '01_merge_positive_sqtls.py' first")
        return

    logging.info("Loading reference genome")
    try:
        fasta_handle = Fasta(fasta_file, key_function=lambda key: key.split(' ')[0], sequence_always_upper=True)
    except FileNotFoundError:
        logging.error(f"FASTA file not found: {fasta_file}")
        return

    logging.info("Processing VCF to create negative pool")
    
    negative_samples = []
    vcf_line_count = 0
    snp_count = 0
    excluded_positive_count = 0
    ref_mismatch_count = 0
    boundary_error_count = 0
    sequence_error_count = 0
    
    with gzip.open(vcf_file, 'rt') as f:
        for line in tqdm(f, desc="Processing VCF"):
            if line.startswith('#'):
                continue
                
            vcf_line_count += 1
            parts = line.strip().split('\t')
            
            if len(parts) < 5:
                continue
                
            chrom, pos_str, _, ref, alt_field = parts[0], parts[1], parts[2], parts[3], parts[4]
            
            try:
                pos = int(pos_str)
            except ValueError:
                continue
                
            alt = alt_field.split(',')[0]
            
            if not (len(ref) == 1 and len(alt) == 1):
                continue
            snp_count += 1
            
            if (str(chrom), pos) in positive_locations:
                excluded_positive_count += 1
                continue
            
            context = get_sequence_context(chrom, pos, fasta_handle, chrom_mapping)
            if not context:
                boundary_error_count += 1
                continue
            
            variant_offset = SEQ_LENGTH // 2
            if context[variant_offset - 1] != ref.upper():
                ref_mismatch_count += 1
                continue
            
            try:
                ref_seq, alt_seq = build_negative_sequences(ref.upper(), alt.upper(), context)
                
                negative_samples.append({
                    'variant_id': f"{chrom}:{pos}_{ref}/{alt}",
                    'ref_sequence': ref_seq,
                    'alt_sequence': alt_seq,
                    'label': 0,
                    'tissue_id': UNKNOWN_TISSUE_ID
                })
                
                if len(negative_samples) % 50000 == 0:
                    logging.info(f"Progress: {len(negative_samples):,} negative samples collected")
                
                if len(negative_samples) >= MAX_NEGATIVE_SAMPLES:
                    logging.info(f"Reached target of {MAX_NEGATIVE_SAMPLES:,} samples, stopping")
                    break
                    
            except Exception:
                sequence_error_count += 1
                continue

    logging.info("VCF Processing Statistics")
    logging.info(f"  Total VCF lines processed: {vcf_line_count:,}")
    if vcf_line_count > 0:
        logging.info(f"  Simple SNPs found: {snp_count:,} ({snp_count/vcf_line_count:.2%})")
    if snp_count > 0:
        logging.info(f"  Excluded (positive locations): {excluded_positive_count:,} ({excluded_positive_count/snp_count:.2%})")
        logging.info(f"  Excluded (reference mismatch): {ref_mismatch_count:,} ({ref_mismatch_count/snp_count:.2%})")
        logging.info(f"  Excluded (boundary errors): {boundary_error_count:,} ({boundary_error_count/snp_count:.2%})")
    logging.info(f"  Excluded (sequence errors): {sequence_error_count:,}")
    logging.info(f"  Final high-quality negatives: {len(negative_samples):,}")
    
    if not negative_samples:
        logging.error("No negative samples found after filtering")
        return
    
    final_df = pd.DataFrame(negative_samples)
    cols_to_save = ['variant_id', 'ref_sequence', 'alt_sequence', 'label', 'tissue_id']
    final_df[cols_to_save].to_csv(output_file, sep='\t', index=False, header=False)
    
    logging.info(f"✅ Output saved to: {output_file}")
    logging.info(f"Dataset contains {len(final_df):,} samples with 100% reference concordance")

if __name__ == '__main__':
    create_negative_pool()

