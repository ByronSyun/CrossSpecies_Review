#!/usr/bin/env python3
import os
import json
import argparse
import time
import logging
import numpy as np
import pandas as pd
from tqdm import tqdm
import torch
from torch.utils.data import DataLoader, Dataset

logging.basicConfig(level=logging.INFO, format='%(asctime)s - [%(levelname)s] - %(message)s')

class SpliceVarDBDataset(Dataset):
    def __init__(self, df):
        self.df = df

    def __len__(self):
        return len(self.df)

    def __getitem__(self, idx):
        row = self.df.iloc[idx]
        variant_id = row['variant_id']
        ref_seq = row['sequence']  # Reference sequence
        
        # Generate alternative sequence by replacing center position
        center_idx = len(ref_seq) // 2
        alt_seq = list(ref_seq)
        alt_seq[center_idx] = row['alt']
        alt_seq = "".join(alt_seq)
        
        return variant_id, ref_seq, alt_seq

def try_load_model(model_ids):
    """Try loading different model versions until one works"""
    for model_id in model_ids:
        try:
            logging.info(f"Trying model: {model_id}")
            from transformers import AutoTokenizer, AutoModel
            
            tokenizer = AutoTokenizer.from_pretrained(model_id, trust_remote_code=True)
            
            # Try loading with different configurations
            try:
                model = AutoModel.from_pretrained(model_id, trust_remote_code=True)
            except Exception as e:
                logging.warning(f"Standard loading failed for {model_id}: {e}")
                # Try with torch_dtype specification
                try:
                    model = AutoModel.from_pretrained(
                        model_id, 
                        trust_remote_code=True,
                        torch_dtype=torch.float32,
                        low_cpu_mem_usage=True
                    )
                except Exception as e2:
                    logging.warning(f"Alternative loading failed for {model_id}: {e2}")
                    continue
            
            logging.info(f"✅ Successfully loaded model: {model_id}")
            return tokenizer, model, model_id
            
        except Exception as e:
            logging.warning(f"Failed to load {model_id}: {e}")
            continue
    
    raise RuntimeError("Failed to load any Nucleotide Transformer model")

@torch.no_grad()
def embed_batch(model, tokenizer, sequences, device, max_length=1000):
    """Embed sequences with error handling"""
    try:
        # Truncate sequences if too long
        truncated_seqs = [seq[:max_length] for seq in sequences]
        
        tokens = tokenizer(
            truncated_seqs, 
            return_tensors='pt', 
            padding=True, 
            truncation=True,
            max_length=max_length
        )
        tokens = {k: v.to(device) for k, v in tokens.items()}
        
        outputs = model(**tokens)
        last_hidden = outputs.last_hidden_state
        emb = last_hidden.mean(dim=1)
        return emb
        
    except Exception as e:
        logging.warning(f"Batch embedding failed: {e}")
        # Fallback: process sequences one by one
        embeddings = []
        for seq in sequences:
            try:
                truncated_seq = seq[:max_length]
                tokens = tokenizer(
                    truncated_seq, 
                    return_tensors='pt', 
                    truncation=True,
                    max_length=max_length
                )
                tokens = {k: v.to(device) for k, v in tokens.items()}
                outputs = model(**tokens)
                emb = outputs.last_hidden_state.mean(dim=1)
                embeddings.append(emb)
            except:
                # If individual sequence fails, use zero embedding
                embeddings.append(torch.zeros(1, model.config.hidden_size, device=device))
        
        return torch.cat(embeddings, dim=0)

def score_variants(args):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    logging.info(f"Device: {device}")

    # Try multiple model versions in order of preference
    model_candidates = [
        "InstaDeepAI/nucleotide-transformer-500m-1000g",  # v1, most stable
        "InstaDeepAI/nucleotide-transformer-2.5b-1000g",  # v1, larger
        "InstaDeepAI/nucleotide-transformer-500m-human-ref",  # v1, human-focused
        "InstaDeepAI/nucleotide-transformer-v2-500m-multi-species",  # v2, potentially problematic
    ]
    
    tokenizer, model, used_model_id = try_load_model(model_candidates)
    model.to(device)
    model.eval()

    # Load SpliceVarDB data
    logging.info(f"Loading SpliceVarDB data from: {args.input_tsv}")
    df = pd.read_csv(args.input_tsv, sep='\t')
    
    # Verify required columns for SpliceVarDB format
    required = ['variant_id', 'sequence', 'ref', 'alt', 'label']
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns for SpliceVarDB format: {missing}")

    logging.info(f"Loaded {len(df)} variants")
    
    # Log sequence length statistics
    seq_lens = df['sequence'].str.len()
    logging.info(f"Sequence lengths - mean={seq_lens.mean():.0f}, max={seq_lens.max()}, min={seq_lens.min()}")

    dataset = SpliceVarDBDataset(df)
    loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=False, num_workers=0)

    # Output files
    out_jsonl = f"{args.output_prefix}_scores.jsonl"
    out_csv = f"{args.output_prefix}_scores.csv"
    os.makedirs(os.path.dirname(args.output_prefix), exist_ok=True)

    rows = []
    t0 = time.time()
    n = 0
    errors = 0
    
    logging.info("Starting variant scoring...")
    for vids, ref_list, alt_list in tqdm(loader, total=len(loader)):
        try:
            # Use conservative max_length to avoid memory issues
            max_len = min(1000, max(len(s) for s in ref_list + alt_list))
            
            ref_emb = embed_batch(model, tokenizer, list(ref_list), device, max_len)
            alt_emb = embed_batch(model, tokenizer, list(alt_list), device, max_len)
            
            # Compute distances
            cos = torch.nn.functional.cosine_similarity(ref_emb, alt_emb, dim=1)
            score_cosine = (1.0 - cos).cpu().numpy()
            score_l2 = torch.linalg.vector_norm(ref_emb - alt_emb, ord=2, dim=1).cpu().numpy()

            for vid, sc, sl2 in zip(vids, score_cosine, score_l2):
                rows.append({
                    'variant_id': vid,
                    'score_cosine': float(sc),
                    'score_l2': float(sl2)
                })
                n += 1

        except Exception as e:
            logging.error(f"Error processing batch: {e}")
            errors += 1
            # Add zero scores for failed batch
            for vid in vids:
                rows.append({
                    'variant_id': vid,
                    'score_cosine': 0.0,
                    'score_l2': 0.0
                })
                n += 1

        if n % 1000 == 0:
            logging.info(f"Processed {n} variants, {errors} batch errors...")

    # Write outputs
    logging.info("Writing outputs...")
    with open(out_jsonl, 'w') as f:
        for r in rows:
            f.write(json.dumps(r) + '\n')

    pd.DataFrame(rows).to_csv(out_csv, index=False)
    
    elapsed = time.time() - t0
    
    # Summary statistics
    scores = [r['score_cosine'] for r in rows]
    non_zero_scores = sum(1 for s in scores if s > 0)
    
    logging.info(f"✅ Completed! Used model: {used_model_id}")
    logging.info(f"Processed {n} variants in {elapsed:.1f}s ({n/elapsed:.1f} variants/sec)")
    logging.info(f"Batch errors: {errors}")
    logging.info(f"Non-zero scores: {non_zero_scores}/{len(scores)} ({non_zero_scores/len(scores)*100:.1f}%)")
    logging.info(f"Score range: [{min(scores):.4f}, {max(scores):.4f}]")
    logging.info(f"Outputs: {out_jsonl} and {out_csv}")

def main():
    parser = argparse.ArgumentParser(
        description="Nucleotide Transformer variant scoring for SpliceVarDB format"
    )
    parser.add_argument('--input_tsv', required=True, 
                       help='SpliceVarDB TSV file with variant_id, sequence, ref, alt, label columns')
    parser.add_argument('--output_prefix', required=True, 
                       help='Output file prefix (without extension)')
    parser.add_argument('--batch_size', type=int, default=2,
                       help='Batch size for processing (default: 2)')
    
    args = parser.parse_args()
    
    logging.info("=== Nucleotide Transformer SpliceVarDB Scoring ===")
    logging.info(f"Input: {args.input_tsv}")
    logging.info(f"Output prefix: {args.output_prefix}")
    logging.info(f"Batch size: {args.batch_size}")
    
    score_variants(args)

if __name__ == '__main__':
    main()
