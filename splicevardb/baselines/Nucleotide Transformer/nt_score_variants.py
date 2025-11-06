#!/usr/bin/env python3
import os
import json
import argparse
import time
import logging
import numpy as np
import pandas as pd
from tqdm import tqdm
from transformers import AutoTokenizer, AutoModel
import torch
from torch.utils.data import DataLoader, Dataset

logging.basicConfig(level=logging.INFO, format='%(asctime)s - [%(levelname)s] - %(message)s')

class VariantSeqDataset(Dataset):
    def __init__(self, df, window, tokenizer):
        self.df = df
        self.window = window
        self.tokenizer = tokenizer

    def __len__(self):
        return len(self.df)

    def __getitem__(self, idx):
        row = self.df.iloc[idx]
        vid = row['variant_id']
        ref_seq = row['ref_seq']
        alt_seq = row['alt_seq']
        return vid, ref_seq, alt_seq

@torch.no_grad()
def embed_batch(model, tokenizer, sequences, device):
    """Embed a batch of sequences using the model"""
    tokens = tokenizer(sequences, return_tensors='pt', padding=True, truncation=True)
    tokens = {k: v.to(device) for k, v in tokens.items()}
    outputs = model(**tokens)
    # Mean pooling over sequence length
    last_hidden = outputs.last_hidden_state
    emb = last_hidden.mean(dim=1)  # [B, D]
    return emb

@torch.no_grad()
def score_variants(args):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    logging.info(f"Device: {device}")

    # Load model and tokenizer
    logging.info(f"Loading model: {args.model_id}")
    tokenizer = AutoTokenizer.from_pretrained(args.model_id, trust_remote_code=True)
    model = AutoModel.from_pretrained(args.model_id, trust_remote_code=True)
    model.to(device)
    model.eval()

    # Load data
    logging.info(f"Loading data from: {args.input_tsv}")
    df = pd.read_csv(args.input_tsv, sep='\t')
    required = ['variant_id', 'ref_seq', 'alt_seq']
    for c in required:
        if c not in df.columns:
            raise ValueError(f"Missing required column: {c}")

    logging.info(f"Loaded {len(df)} variants")

    # Create dataset and dataloader
    dataset = VariantSeqDataset(df, args.window, tokenizer)
    loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=False, num_workers=args.num_workers)

    # Prepare output files
    out_jsonl = f"{args.output_prefix}_scores.jsonl"
    out_csv = f"{args.output_prefix}_scores.csv"
    os.makedirs(os.path.dirname(args.output_prefix), exist_ok=True)

    # Process variants
    rows = []
    t0 = time.time()
    n = 0
    
    logging.info("Starting variant scoring...")
    for vids, ref_list, alt_list in tqdm(loader, total=len(loader), desc="Scoring variants"):
        # Get embeddings for reference and alternative sequences
        ref_emb = embed_batch(model, tokenizer, list(ref_list), device)
        alt_emb = embed_batch(model, tokenizer, list(alt_list), device)
        
        # Compute similarity scores
        cos = torch.nn.functional.cosine_similarity(ref_emb, alt_emb, dim=1)
        score_cosine = (1.0 - cos).cpu().numpy()  # Cosine distance: 1 - cosine_similarity
        score_l2 = torch.linalg.vector_norm(ref_emb - alt_emb, ord=2, dim=1).cpu().numpy()

        # Store results
        for vid, sc, sl2 in zip(vids, score_cosine, score_l2):
            rec = {
                'variant_id': vid,
                'score_cosine': float(sc),
                'score_l2': float(sl2)
            }
            rows.append(rec)
            n += 1

        # Periodic logging
        if n % 1000 == 0:
            logging.info(f"Processed {n} variants...")

    # Write outputs
    logging.info("Writing outputs...")
    with open(out_jsonl, 'w') as f:
        for r in rows:
            f.write(json.dumps(r) + '\n')

    pd.DataFrame(rows).to_csv(out_csv, index=False)
    
    elapsed = time.time() - t0
    logging.info(f"✅ Completed! Processed {n} variants in {elapsed:.1f}s")
    logging.info(f"Outputs: {out_jsonl} and {out_csv}")
    
    # Print summary stats
    scores = [r['score_cosine'] for r in rows]
    logging.info(f"Score summary - Mean: {np.mean(scores):.4f}, Std: {np.std(scores):.4f}, Range: [{np.min(scores):.4f}, {np.max(scores):.4f}]")

def main():
    parser = argparse.ArgumentParser(
        description="Score variants using Nucleotide Transformer embeddings"
    )
    parser.add_argument('--model_id', required=True, 
                       help='HuggingFace model ID (e.g., InstaDeepAI/nucleotide-transformer-v2-500m-multi-species)')
    parser.add_argument('--input_tsv', required=True, 
                       help='Input TSV file with variant_id, ref_seq, alt_seq columns')
    parser.add_argument('--output_prefix', required=True, 
                       help='Output file prefix (without extension)')
    parser.add_argument('--window', type=int, default=8192,
                       help='Sequence window size (default: 8192)')
    parser.add_argument('--batch_size', type=int, default=4,
                       help='Batch size for processing (default: 4)')
    parser.add_argument('--num_workers', type=int, default=2,
                       help='Number of data loader workers (default: 2)')
    
    args = parser.parse_args()
    
    logging.info("=== Nucleotide Transformer Variant Scoring ===")
    logging.info(f"Model: {args.model_id}")
    logging.info(f"Input: {args.input_tsv}")
    logging.info(f"Output prefix: {args.output_prefix}")
    logging.info(f"Batch size: {args.batch_size}")
    
    score_variants(args)

if __name__ == '__main__':
    main()