#!/usr/bin/env python3
"""
SpliceBERT Variant Effect Scoring Script
Based on SpliceBERT paper methodology: KL-divergence scoring for variant effect prediction
"""

import os
import json
import argparse
import time
import logging
import numpy as np
import pandas as pd
from tqdm import tqdm
from transformers import AutoTokenizer, AutoModelForMaskedLM
import torch
import torch.nn.functional as F
from torch.utils.data import DataLoader, Dataset

logging.basicConfig(level=logging.INFO, format='%(asctime)s - [%(levelname)s] - %(message)s')

class SpliceBERTVariantDataset(Dataset):
    def __init__(self, df, tokenizer, window_size=1024):
        self.df = df
        self.tokenizer = tokenizer
        self.window_size = window_size

    def __len__(self):
        return len(self.df)

    def __getitem__(self, idx):
        row = self.df.iloc[idx]
        variant_id = row['variant_id']
        full_ref_seq = row['sequence']  # Full reference sequence (8192bp)
        
        # Truncate to SpliceBERT's maximum length (1024nt)
        # Center the variant position within the window
        center_idx = len(full_ref_seq) // 2
        half_window = self.window_size // 2
        
        start_pos = max(0, center_idx - half_window)
        end_pos = min(len(full_ref_seq), start_pos + self.window_size)
        
        # Adjust start_pos if we're at the end of sequence
        if end_pos - start_pos < self.window_size:
            start_pos = max(0, end_pos - self.window_size)
        
        ref_seq = full_ref_seq[start_pos:end_pos]
        
        # Generate alternative sequence by replacing center position
        variant_pos_in_window = center_idx - start_pos
        alt_seq = list(ref_seq)
        alt_seq[variant_pos_in_window] = row['alt']
        alt_seq = "".join(alt_seq)
        
        return variant_id, ref_seq, alt_seq, variant_pos_in_window

def prepare_sequence_for_splicebert(seq):
    """
    Prepare DNA sequence for SpliceBERT input
    - Convert T to U (RNA format)
    - Add spaces between nucleotides
    - Ensure uppercase
    """
    seq = seq.upper().replace("T", "U")
    seq = ' '.join(list(seq))
    return seq

def kl_divergence(ref_logits, alt_logits):
    """
    Compute KL divergence between reference and alternative logits
    KL(alt || ref) = sum(alt * log(alt / ref))
    """
    ref_probs = F.softmax(ref_logits, dim=-1)
    alt_probs = F.softmax(alt_logits, dim=-1)
    
    # Add small epsilon to avoid log(0)
    eps = 1e-8
    ref_probs = ref_probs + eps
    alt_probs = alt_probs + eps
    
    kl_div = (alt_probs * torch.log(alt_probs / ref_probs)).sum(dim=-1)
    return kl_div

@torch.no_grad()
def score_variants(args):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    logging.info(f"Device: {device}")

    # Load SpliceBERT model and tokenizer
    logging.info(f"Loading SpliceBERT model from: {args.model_path}")
    tokenizer = AutoTokenizer.from_pretrained(args.model_path)
    model = AutoModelForMaskedLM.from_pretrained(args.model_path)
    model.to(device)
    model.eval()
    
    logging.info(f"✅ Successfully loaded SpliceBERT model")
    logging.info(f"Model vocab size: {tokenizer.vocab_size}")

    # Load data
    logging.info(f"Loading variant data from: {args.input_tsv}")
    df = pd.read_csv(args.input_tsv, sep='\t')
    required = ['variant_id', 'sequence', 'ref', 'alt']
    for c in required:
        if c not in df.columns:
            raise ValueError(f"Missing required column: {c}")

    logging.info(f"Loaded {len(df)} variants")
    
    # Ensure sequence length is consistent
    seq_lengths = df['sequence'].str.len()
    logging.info(f"Sequence lengths - mean={seq_lengths.mean():.0f}, max={seq_lengths.max():.0f}, min={seq_lengths.min():.0f}")

    # Create dataset and dataloader
    dataset = SpliceBERTVariantDataset(df, tokenizer, args.window_size)
    loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=False, num_workers=args.num_workers)

    # Prepare output files
    out_jsonl = f"{args.output_prefix}_scores.jsonl"
    out_csv = f"{args.output_prefix}_scores.csv"
    os.makedirs(os.path.dirname(args.output_prefix), exist_ok=True)

    # Process variants
    rows = []
    t0 = time.time()
    n = 0
    batch_errors = 0
    flanking_window = args.flanking_window
    
    logging.info(f"Starting variant scoring with flanking window: ±{flanking_window} nucleotides...")
    
    for batch_data in tqdm(loader, total=len(loader), desc="Scoring variants"):
        try:
            vids, ref_list, alt_list, variant_positions = batch_data
            batch_kl_scores = []
            
            for ref_seq, alt_seq, var_pos in zip(ref_list, alt_list, variant_positions):
                # Prepare sequences for SpliceBERT
                ref_seq_formatted = prepare_sequence_for_splicebert(ref_seq)
                alt_seq_formatted = prepare_sequence_for_splicebert(alt_seq)
                
                # Tokenize sequences
                ref_tokens = tokenizer.encode(ref_seq_formatted, return_tensors='pt').to(device)
                alt_tokens = tokenizer.encode(alt_seq_formatted, return_tensors='pt').to(device)
                
                # Get logits from SpliceBERT
                ref_logits = model(ref_tokens).logits[0]  # Remove batch dimension
                alt_logits = model(alt_tokens).logits[0]  # Remove batch dimension
                
                # Find variant position in tokenized sequence
                # Account for [CLS] token at start and space-separated tokenization
                variant_token_pos = var_pos + 1  # +1 for [CLS] token
                
                # Calculate KL divergence in flanking window
                start_pos = max(1, variant_token_pos - flanking_window)  # Skip [CLS] token
                end_pos = min(len(ref_tokens[0]) - 1, variant_token_pos + flanking_window + 1)  # Skip [SEP] token
                
                # Compute KL divergence for entire sequence (following original implementation)
                full_kl = kl_divergence(ref_logits, alt_logits)  # Shape: [seq_len]
                
                # Extract flanking window around variant (±100nt as in original)
                variant_pos_in_seq = variant_token_pos - 1  # Adjust for [CLS] token
                flanking_start = max(0, variant_pos_in_seq - flanking_window)
                flanking_end = min(len(full_kl), variant_pos_in_seq + flanking_window + 1)
                
                # Get neighboring KL values
                neighboring = full_kl[flanking_start:flanking_end]
                
                # Create distance-weighted KL (following original implementation)
                # Split into before and after variant, then average symmetrically
                if len(neighboring) == 2 * flanking_window + 1:  # Full window available
                    before_variant = neighboring[:flanking_window].flip(0)  # Reverse order using flip
                    after_variant = neighboring[flanking_window + 1:]
                    # Average corresponding positions
                    min_len = min(len(before_variant), len(after_variant))
                    kl_distance = (before_variant[:min_len] + after_variant[:min_len]) / 2
                else:
                    # Partial window, use as-is
                    kl_distance = neighboring
                
                # Apply original formula: sum(log(clip(kl_distance)))
                kl_context_score = torch.sum(torch.log(torch.clamp(kl_distance, min=1e-6))).item()
                
                batch_kl_scores.append(kl_context_score)
            
            # Store results
            for vid, kl_score in zip(vids, batch_kl_scores):
                rec = {
                    'variant_id': vid,
                    'kl_context_score': float(kl_score),
                    'flanking_window': flanking_window
                }
                rows.append(rec)
                n += 1
                
        except Exception as e:
            logging.error(f"Error processing batch for VIDs {vids[0]}...{vids[-1]}: {e}")
            batch_errors += 1
            continue

        # Periodic logging
        if n % 100 == 0 and n > 0:
            logging.info(f"Processed {n} variants, {batch_errors} batch errors...")

    # Write outputs
    logging.info("Writing outputs...")
    with open(out_jsonl, 'w') as f:
        for r in rows:
            f.write(json.dumps(r) + '\n')

    pd.DataFrame(rows).to_csv(out_csv, index=False)
    
    elapsed = time.time() - t0
    logging.info(f"✅ Completed SpliceBERT variant scoring!")
    logging.info(f"Processed {n} variants in {elapsed:.1f}s ({n/elapsed:.1f} variants/sec)")
    logging.info(f"Batch errors: {batch_errors}")
    logging.info(f"Outputs: {out_jsonl} and {out_csv}")
    
    # Print summary stats
    if rows:
        scores = [r['kl_context_score'] for r in rows]
        non_zero_scores = sum(1 for s in scores if abs(s) > 1e-9)
        logging.info(f"Non-zero scores: {non_zero_scores}/{len(rows)} ({non_zero_scores/len(rows)*100:.1f}%)")
        logging.info(f"Score range: [{np.min(scores):.4f}, {np.max(scores):.4f}]")
        logging.info(f"Score summary - Mean: {np.mean(scores):.4f}, Std: {np.std(scores):.4f}, Median: {np.median(scores):.4f}")
    else:
        logging.warning("No scores generated.")

def main():
    parser = argparse.ArgumentParser(
        description="Score variants using SpliceBERT KL-divergence method"
    )
    parser.add_argument('--input_tsv', required=True, 
                       help='Input TSV file with variant_id, sequence, ref, alt columns')
    parser.add_argument('--model_path', required=True,
                       help='Path to SpliceBERT model directory')
    parser.add_argument('--output_prefix', required=True, 
                       help='Output file prefix (without extension)')
    parser.add_argument('--window_size', type=int, default=1024,
                       help='Maximum sequence length for SpliceBERT (default: 1024)')
    parser.add_argument('--flanking_window', type=int, default=100,
                       help='Flanking window size for KL-divergence calculation (default: 100)')
    parser.add_argument('--batch_size', type=int, default=4,
                       help='Batch size for processing (default: 4)')
    parser.add_argument('--num_workers', type=int, default=1,
                       help='Number of data loader workers (default: 1)')
    
    args = parser.parse_args()
    
    logging.info("=== SpliceBERT Variant Effect Scoring ===")
    logging.info(f"Input: {args.input_tsv}")
    logging.info(f"Model: {args.model_path}")
    logging.info(f"Output prefix: {args.output_prefix}")
    logging.info(f"Flanking window: ±{args.flanking_window} nucleotides")
    logging.info(f"Batch size: {args.batch_size}")
    
    score_variants(args)

if __name__ == '__main__':
    main()
