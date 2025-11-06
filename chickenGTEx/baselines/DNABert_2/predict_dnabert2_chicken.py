#!/usr/bin/env python3
"""DNABERT-2 Concatenated Embedding Generation for Chicken Data [WT, MT, diff] (2304D)"""
import pandas as pd
import torch
from transformers import AutoTokenizer, BertModel
import numpy as np
from tqdm import tqdm
import argparse

def get_dnabert2_embeddings(sequences, model, tokenizer, device, batch_size=32):
    """Generate DNABERT-2 mean pooling embeddings."""
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

def main(args):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")
    
    print("Loading DNABERT-2 model...")
    tokenizer = AutoTokenizer.from_pretrained("zhihan1996/DNABERT-2-117M", trust_remote_code=True)
    model = BertModel.from_pretrained("zhihan1996/DNABERT-2-117M", trust_remote_code=True)
    model.to(device)
    model.eval()
    
    print(f"Loading data: {args.input_file}")
    df = pd.read_csv(args.input_file, sep='\t', header=None,
                     names=['variant_id', 'ref_seq', 'alt_seq', 'label'])
    print(f"Variants: {len(df)}")
    
    wt_sequences = df['ref_seq'].tolist()
    mt_sequences = df['alt_seq'].tolist()
    
    wt_embeddings = get_dnabert2_embeddings(wt_sequences, model, tokenizer, device, args.batch_size)
    mt_embeddings = get_dnabert2_embeddings(mt_sequences, model, tokenizer, device, args.batch_size)
    
    diff_embeddings = mt_embeddings - wt_embeddings
    concatenated_embeddings = np.concatenate((wt_embeddings, mt_embeddings, diff_embeddings), axis=1)
    
    labels = df['label'].values.astype(int)
    variant_ids = df['variant_id'].values
    
    np.savez_compressed(
        args.output_file, 
        variant_ids=variant_ids,
        embeddings=concatenated_embeddings,
        labels=labels
    )
    print(f"Saved: {args.output_file} ({len(variant_ids)} variants, {concatenated_embeddings.shape[1]}D)")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate DNABERT-2 concatenated embeddings for chicken variants")
    parser.add_argument("--input_file", type=str, required=True, help="Input TSV file")
    parser.add_argument("--output_file", type=str, required=True, help="Output NPZ file")
    parser.add_argument("--batch_size", type=int, default=16, help="Batch size")
    args = parser.parse_args()
    main(args)

