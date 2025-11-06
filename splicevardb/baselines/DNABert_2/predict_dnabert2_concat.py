#!/usr/bin/env python3
import pandas as pd
import torch
from transformers import AutoTokenizer, BertModel
import numpy as np
from tqdm import tqdm
import argparse

def get_dnabert2_embeddings(sequences, model, tokenizer, device, batch_size=32):
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

def create_mutant_sequences(df):
    wt_sequences = df['sequence'].tolist()
    mut_sequences = []
    for index, row in tqdm(df.iterrows(), total=df.shape[0], desc="Creating Mutant Sequences"):
        wt_seq = row['sequence']
        ref = row['ref']
        alt = row['alt']
        center_index = len(wt_seq) // 2
        if wt_seq[center_index].upper() != ref.upper():
            print(f"Warning: Mismatch at variant {row['variant_id']}")
        mut_seq_list = list(wt_seq)
        mut_seq_list[center_index] = alt
        mut_sequences.append("".join(mut_seq_list))
    return wt_sequences, mut_sequences

def main(args):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")
    
    print("Loading DNABERT-2 model...")
    tokenizer = AutoTokenizer.from_pretrained("zhihan1996/DNABERT-2-117M", trust_remote_code=True)
    model = BertModel.from_pretrained("zhihan1996/DNABERT-2-117M", trust_remote_code=True)
    model.to(device)
    model.eval()
    
    print(f"Loading data from {args.input_file}...")
    df = pd.read_csv(args.input_file, sep='\t')
    print(f"Loaded {len(df)} variants")
    
    wt_sequences, mt_sequences = create_mutant_sequences(df)
    
    print("Generating WT embeddings...")
    wt_embeddings = get_dnabert2_embeddings(wt_sequences, model, tokenizer, device, args.batch_size)
    
    print("Generating MT embeddings...")
    mt_embeddings = get_dnabert2_embeddings(mt_sequences, model, tokenizer, device, args.batch_size)
    
    diff = mt_embeddings - wt_embeddings
    features = np.concatenate([wt_embeddings, mt_embeddings, diff], axis=1)
    
    print(f"Feature shape: {features.shape} (WT + MT + diff)")
    
    labels = df['label'].apply(lambda x: 1 if x == 'splice-altering' else 0).values
    variant_ids = df['variant_id'].values
    
    print(f"Saving to {args.output_file}...")
    np.savez_compressed(
        args.output_file, 
        variant_ids=variant_ids,
        embeddings=features,
        labels=labels
    )
    print("Done.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", type=str, required=True)
    parser.add_argument("--output_file", type=str, required=True)
    parser.add_argument("--batch_size", type=int, default=32)
    args = parser.parse_args()
    main(args)

