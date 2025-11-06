#!/usr/bin/env python3
"""
DNABERT-2 Embedding Generation for Rat Data
Strategy: Concatenate [WT, MT, diff] embeddings (2304D)
"""
import pandas as pd
import torch
from transformers import AutoTokenizer, BertModel
import numpy as np
from tqdm import tqdm
import argparse

def parse_variant_id(variant_id):
    """Parse variant_id to extract ref and alt alleles.
    Format: chr:pos_ref/alt (e.g., 15:55290237_A/G)
    """
    parts = variant_id.split(':')
    if len(parts) == 2:
        chrom = parts[0]
        pos_alleles = parts[1]
        if '_' in pos_alleles:
            pos, alleles = pos_alleles.split('_')
            if '/' in alleles:
                ref, alt = alleles.split('/')
                return ref, alt
    raise ValueError(f"Cannot parse variant_id: {variant_id}")

def get_dnabert2_embeddings(sequences, model, tokenizer, device, batch_size=32):
    """Generate DNABERT-2 mean pooling embeddings."""
    all_embeddings = []
    for i in tqdm(range(0, len(sequences), batch_size), desc="Generating Embeddings"):
        batch_sequences = sequences[i:i+batch_size]
        inputs = tokenizer(batch_sequences, padding="longest", return_tensors='pt', 
                         truncation=True, max_length=512)
        inputs = {key: val.to(device) for key, val in inputs.items()}
        
        with torch.no_grad():
            hidden_states = model(**inputs)[0]
        
        batch_embeddings = torch.mean(hidden_states, dim=1)
        all_embeddings.append(batch_embeddings.cpu().numpy())
    
    return np.vstack(all_embeddings)

def create_mutant_sequences(df):
    """Create mutant sequences by replacing the center position."""
    wt_sequences = []
    mt_sequences = []
    
    for index, row in tqdm(df.iterrows(), total=df.shape[0], desc="Creating Mutant Sequences"):
        wt_seq = row['ref_seq'] if 'ref_seq' in row else row['sequence']
        variant_id = row['variant_id']
        
        try:
            ref, alt = parse_variant_id(variant_id)
        except ValueError as e:
            print(f"Warning: {e}, skipping variant")
            continue
        
        center_index = len(wt_seq) // 2
        
        if wt_seq[center_index].upper() != ref.upper():
            print(f"Warning: Mismatch at variant {variant_id}: "
                  f"expected {ref}, got {wt_seq[center_index]}")
        
        wt_sequences.append(wt_seq)
        
        mut_seq_list = list(wt_seq)
        mut_seq_list[center_index] = alt
        mt_sequences.append("".join(mut_seq_list))
    
    return wt_sequences, mt_sequences

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
                     names=['variant_id', 'ref_seq', 'alt_seq', 'label', 'chrom'])
    print(f"Variants: {len(df)}")
    
    wt_sequences, mt_sequences = create_mutant_sequences(df)
    
    wt_embeddings = get_dnabert2_embeddings(wt_sequences, model, tokenizer, device, args.batch_size)
    mt_embeddings = get_dnabert2_embeddings(mt_sequences, model, tokenizer, device, args.batch_size)
    
    diff_embeddings = mt_embeddings - wt_embeddings
    concatenated_embeddings = np.concatenate((wt_embeddings, mt_embeddings, diff_embeddings), axis=1)
    
    labels = df['label'].apply(lambda x: 1 if x == 1 or x == '1' else 0).values[:len(wt_sequences)]
    variant_ids = df['variant_id'].values[:len(wt_sequences)]
    
    np.savez_compressed(
        args.output_file, 
        variant_ids=variant_ids,
        embeddings=concatenated_embeddings,
        labels=labels
    )
    print(f"Saved: {args.output_file} ({len(variant_ids)} variants, {concatenated_embeddings.shape[1]}D)")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate DNABERT-2 concatenated embeddings for rat variants")
    parser.add_argument("--input_file", type=str, required=True, help="Input TSV file")
    parser.add_argument("--output_file", type=str, required=True, help="Output NPZ file")
    parser.add_argument("--batch_size", type=int, default=16, help="Batch size for inference")
    args = parser.parse_args()
    main(args)

