#!/usr/bin/env python3
"""Generate DNABERT-2 differential embeddings (ALT - REF) for chicken variants."""
import pandas as pd
import torch
from transformers import AutoTokenizer, BertModel
import numpy as np
from tqdm import tqdm
import argparse

def get_dnabert2_embeddings(sequences, model, tokenizer, device, batch_size=32):
    """Generate embeddings for DNA sequences using DNABERT-2."""
    all_embeddings = []
    
    for i in tqdm(range(0, len(sequences), batch_size), desc="Generating Embeddings"):
        batch_sequences = sequences[i:i+batch_size]
        inputs = tokenizer(batch_sequences, padding="longest", return_tensors='pt', truncation=True, max_length=512)
        inputs = {key: val.to(device) for key, val in inputs.items()}
        
        with torch.no_grad():
            hidden_states = model(**inputs)[0]
        
        batch_embeddings = torch.mean(hidden_states, dim=1)
        all_embeddings.append(batch_embeddings.cpu().numpy())
        
    return np.vstack(all_embeddings)

def load_chicken_data(input_file):
    """Load chicken data from TSV file (no header: variant_id, ref_seq, alt_seq, label)."""
    df = pd.read_csv(input_file, sep='\t', header=None,
                     names=['variant_id', 'ref_seq', 'alt_seq', 'label'],
                     dtype={'variant_id': str})
    return df

def main(args):
    """Generate DIFFERENTIAL embeddings (ALT - REF) for chicken data."""
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    print("Loading DNABERT-2 model...")
    tokenizer = AutoTokenizer.from_pretrained("zhihan1996/DNABERT-2-117M", trust_remote_code=True)
    model = BertModel.from_pretrained("zhihan1996/DNABERT-2-117M", trust_remote_code=True)
    model.to(device)
    model.eval()
    print("Model loaded")

    print(f"Loading data from {args.input_file}...")
    df = load_chicken_data(args.input_file)
    print(f"Loaded {len(df)} variants")
    
    wt_sequences = df['ref_seq'].tolist()
    mt_sequences = df['alt_seq'].tolist()
    
    print("Generating embeddings for REF sequences...")
    wt_embeddings = get_dnabert2_embeddings(wt_sequences, model, tokenizer, device, args.batch_size)
    
    print("Generating embeddings for ALT sequences...")
    mt_embeddings = get_dnabert2_embeddings(mt_sequences, model, tokenizer, device, args.batch_size)
    
    diff_embeddings = mt_embeddings - wt_embeddings
    
    labels = df['label'].values.astype(int)
    variant_ids = df['variant_id'].values

    print(f"Saving to {args.output_file}...")
    np.savez_compressed(
        args.output_file, 
        variant_ids=variant_ids,
        embeddings=diff_embeddings,
        labels=labels
    )
    print(f"Saved {len(variant_ids)} variants, embedding shape: {diff_embeddings.shape}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate DNABERT-2 differential embeddings for chicken variants")
    parser.add_argument("--input_file", type=str, required=True, help="Input TSV file")
    parser.add_argument("--output_file", type=str, default="dnabert2_chicken_logistic_embeddings.npz", help="Output NPZ file")
    parser.add_argument("--batch_size", type=int, default=32, help="Batch size")
    args = parser.parse_args()
    main(args)

