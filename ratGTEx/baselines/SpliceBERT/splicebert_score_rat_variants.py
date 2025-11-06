#!/usr/bin/env python3
"""
SpliceBERT Variant Effect Scoring Script for Rat Data
Based on SpliceBERT paper methodology: KL-divergence scoring for variant effect prediction
Adapted for rat ratGTEx dataset
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

class SpliceBERTRatDataset(Dataset):
    def __init__(self, df, tokenizer, window_size=1024):
        self.df = df
        self.tokenizer = tokenizer
        self.window_size = window_size

    def __len__(self):
        return len(self.df)

    def __getitem__(self, idx):
        row = self.df.iloc[idx]
        variant_id = row['variant_id']
        full_ref_seq = row['ref_seq']  # Reference sequence from rat dataset
        ref_allele = row['ref']
        alt_allele = row['alt']
        
        # Truncate to SpliceBERT's maximum length (1024nt)
        # For rat data, the sequence is already 8192bp centered on the variant
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
        
        # Verify reference allele matches
        if variant_pos_in_window < len(ref_seq) and ref_seq[variant_pos_in_window].upper() == ref_allele.upper():
            alt_seq = list(ref_seq)
            alt_seq[variant_pos_in_window] = alt_allele
            alt_seq = "".join(alt_seq)
        else:
            # If reference doesn't match, skip this variant
            logging.warning(f"Reference allele mismatch for {variant_id}: expected {ref_allele}, got {ref_seq[variant_pos_in_window] if variant_pos_in_window < len(ref_seq) else 'N/A'}")
            return None
        
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
def score_rat_variants(args):
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

    # Load rat dataset
    logging.info(f"Loading rat dataset from: {args.input_tsv}")
    df = pd.read_csv(args.input_tsv, sep='\t', header=None, 
                     names=['variant_id', 'ref_seq', 'alt_seq', 'label', 'tissue_id'])
    
    # Parse variant information from variant_id (format: chr:pos_ref/alt)
    # Handle potential parsing failures gracefully
    variant_parts = df['variant_id'].str.extract(r'(.+):(\d+)_(.+)/(.+)')
    
    # Check for parsing failures
    failed_parsing = variant_parts.isnull().any(axis=1)
    if failed_parsing.sum() > 0:
        logging.warning(f"Failed to parse {failed_parsing.sum()} variant IDs. Dropping these variants.")
        logging.info(f"Example failed variant IDs: {df.loc[failed_parsing, 'variant_id'].head().tolist()}")
        
        # Drop rows with failed parsing
        valid_mask = ~failed_parsing
        df = df[valid_mask].reset_index(drop=True)
        variant_parts = variant_parts[valid_mask].reset_index(drop=True)
    
    df['chrom'] = variant_parts[0]
    df['pos'] = variant_parts[1].astype(int)
    df['ref'] = variant_parts[2]
    df['alt'] = variant_parts[3]
    
    logging.info(f"Loaded {len(df)} rat variants")

    # Create dataset
    dataset = SpliceBERTRatDataset(df, tokenizer, window_size=args.window_size)
    dataloader = DataLoader(dataset, batch_size=args.batch_size, shuffle=False, 
                           collate_fn=lambda x: [item for item in x if item is not None])

    # Prepare output files
    out_jsonl = f"{args.output_prefix}_scores.jsonl"
    out_csv = f"{args.output_prefix}_scores.csv"
    
    logging.info(f"Output files: {out_jsonl}, {out_csv}")

    rows = []
    n = 0
    batch_errors = 0
    flanking_window = args.flanking_window

    for batch_data in tqdm(dataloader, desc="Processing batches"):
        if not batch_data:  # Skip empty batches
            continue
            
        try:
            vids, ref_seqs, alt_seqs, var_positions = zip(*batch_data)
            
            # Prepare sequences for SpliceBERT
            ref_inputs = [prepare_sequence_for_splicebert(seq) for seq in ref_seqs]
            alt_inputs = [prepare_sequence_for_splicebert(seq) for seq in alt_seqs]
            
            # Tokenize sequences
            ref_tokens = tokenizer(ref_inputs, return_tensors='pt', padding=True, truncation=True, max_length=1024)
            alt_tokens = tokenizer(alt_inputs, return_tensors='pt', padding=True, truncation=True, max_length=1024)
            
            # Move to device
            ref_tokens = {k: v.to(device) for k, v in ref_tokens.items()}
            alt_tokens = {k: v.to(device) for k, v in alt_tokens.items()}
            
            # Get model predictions
            ref_outputs = model(**ref_tokens)
            alt_outputs = model(**alt_tokens)
            
            ref_logits = ref_outputs.logits  # [batch_size, seq_len, vocab_size]
            alt_logits = alt_outputs.logits  # [batch_size, seq_len, vocab_size]
            
            batch_kl_scores = []
            for i, (vid, var_pos) in enumerate(zip(vids, var_positions)):
                # Find variant position in tokenized sequence
                # Account for [CLS] token at start and space-separated tokenization
                variant_token_pos = var_pos + 1  # +1 for [CLS] token
                
                # Calculate KL divergence in flanking window
                start_pos = max(1, variant_token_pos - flanking_window)  # Skip [CLS] token
                end_pos = min(len(ref_tokens['input_ids'][i]) - 1, variant_token_pos + flanking_window + 1)  # Skip [SEP] token
                
                # Compute KL divergence for entire sequence (following original implementation)
                full_kl = kl_divergence(ref_logits[i], alt_logits[i])  # Shape: [seq_len]
                
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
    
    df_out = pd.DataFrame(rows)
    df_out.to_csv(out_csv, index=False)
    
    logging.info(f"✅ Completed! Processed {n} variants with {batch_errors} batch errors")
    logging.info(f"Results saved to: {out_jsonl}, {out_csv}")
    
    return out_csv, out_jsonl

def main():
    parser = argparse.ArgumentParser(description="Score rat variants using SpliceBERT")
    parser.add_argument('--input_tsv', required=True, help='Input TSV file with rat variants')
    parser.add_argument('--model_path', required=True, help='Path to SpliceBERT model directory')
    parser.add_argument('--output_prefix', required=True, help='Output file prefix')
    parser.add_argument('--flanking_window', type=int, default=100, help='Flanking window size (default: 100)')
    parser.add_argument('--window_size', type=int, default=1024, help='Sequence window size (default: 1024)')
    parser.add_argument('--batch_size', type=int, default=4, help='Batch size (default: 4)')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.input_tsv):
        raise FileNotFoundError(f"Input file not found: {args.input_tsv}")
    
    if not os.path.exists(args.model_path):
        raise FileNotFoundError(f"Model path not found: {args.model_path}")
    
    # Create output directory
    output_dir = os.path.dirname(args.output_prefix)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    score_rat_variants(args)

if __name__ == "__main__":
    main()
