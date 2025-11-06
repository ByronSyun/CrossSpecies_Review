#!/usr/bin/env python3
import argparse
import pandas as pd
import torch
import yaml
import numpy as np
from tqdm import tqdm
import torch.nn.functional as F

from vortex.model.model import StripedHyena
from vortex.model.tokenizer import CharLevelTokenizer
from vortex.model.utils import dotdict, load_checkpoint


def get_evo2_embeddings(seqs, model, tokenizer, device):
    """
    Extract hidden state embeddings from Evo2 model.
    Returns embeddings of shape (batch_size, hidden_dim).
    """
    if not seqs:
        return np.array([])
    
    tokenized_seqs = [tokenizer.tokenize(s) for s in seqs]
    max_len = max(len(s) for s in tokenized_seqs) if tokenized_seqs else 0
    
    if max_len == 0:
        return np.zeros((len(seqs), 1))
    
    padded_tokens = []
    for s in tokenized_seqs:
        pad_len = max_len - len(s)
        padded_tokens.append(
            F.pad(torch.tensor(s, dtype=torch.long), (0, pad_len), 'constant', 0)
        )
    input_tensor = torch.stack(padded_tokens).to(device)
    
    with torch.no_grad():
        try:
            outputs = model(input_tensor, output_hidden_states=True)
            if isinstance(outputs, tuple) and len(outputs) > 1:
                hidden_states = outputs[1]
            else:
                hidden_states = outputs[0] if isinstance(outputs, tuple) else outputs
        except:
            logits, _ = model(input_tensor)
            hidden_states = logits
    
    # Mean pooling over sequence length if 3D, otherwise use as-is
    if len(hidden_states.shape) == 3:
        embeddings = torch.mean(hidden_states, dim=1)
    else:
        embeddings = hidden_states
    
    return embeddings.cpu().float().numpy()


def create_mutant_sequences(df):
    """Create mutant sequences by substituting alt allele at center position."""
    wt_sequences = df['sequence'].tolist()
    mut_sequences = []
    
    for _, row in tqdm(df.iterrows(), total=len(df), desc="Creating Mutant Sequences"):
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
    print(f"Loading Evo2 model from {args.checkpoint_path}...")
    config = dotdict(yaml.load(open(args.config_path), Loader=yaml.FullLoader))
    tokenizer = CharLevelTokenizer(config.vocab_size)
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")
    
    model = StripedHyena(config)
    load_checkpoint(model, args.checkpoint_path)
    model.to(device)
    model.eval()
    
    print(f"Loading data from {args.input_file}...")
    df = pd.read_csv(args.input_file, sep='\t')
    print(f"Loaded {len(df)} variants")
    
    wt_sequences, mt_sequences = create_mutant_sequences(df)
    
    print("Generating WT embeddings...")
    wt_embeddings_list = []
    for i in tqdm(range(0, len(wt_sequences), args.batch_size), desc="WT Embeddings"):
        batch_seqs = wt_sequences[i:i+args.batch_size]
        batch_emb = get_evo2_embeddings(batch_seqs, model, tokenizer, device)
        wt_embeddings_list.append(batch_emb)
    wt_embeddings = np.vstack(wt_embeddings_list)
    
    print("Generating MT embeddings...")
    mt_embeddings_list = []
    for i in tqdm(range(0, len(mt_sequences), args.batch_size), desc="MT Embeddings"):
        batch_seqs = mt_sequences[i:i+args.batch_size]
        batch_emb = get_evo2_embeddings(batch_seqs, model, tokenizer, device)
        mt_embeddings_list.append(batch_emb)
    mt_embeddings = np.vstack(mt_embeddings_list)
    
    # Concatenate [WT, MT, diff]
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
        labels=labels,
        embedding_dim=wt_embeddings.shape[1]
    )
    print(f"Saved: {args.output_file}")
    print(f"Embedding dimension: {wt_embeddings.shape[1]}, Total features: {features.shape[1]} (WT + MT + diff)")
    print("Done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract Evo2 concatenated embeddings for variant prediction.")
    parser.add_argument("--config_path", required=True, help="Path to Evo2 config YAML file.")
    parser.add_argument("--checkpoint_path", required=True, help="Path to Evo2 checkpoint (.pt file).")
    parser.add_argument("--input_file", required=True, help="Path to input TSV file with variants and sequences.")
    parser.add_argument("--output_file", required=True, help="Path to save output NPZ file with embeddings.")
    parser.add_argument("--batch_size", type=int, default=1, help="Batch size for embedding extraction.")
    args = parser.parse_args()
    main(args)

