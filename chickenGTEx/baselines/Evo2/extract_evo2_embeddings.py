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
    
    if len(hidden_states.shape) == 3:
        embeddings = torch.mean(hidden_states, dim=1)
    else:
        embeddings = hidden_states
    
    return embeddings.cpu().float().numpy()


def main(args):
    config = dotdict(yaml.load(open(args.config_path), Loader=yaml.FullLoader))
    tokenizer = CharLevelTokenizer(config.vocab_size)
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    
    print(f"Device: {device}")
    print(f"Loading model from: {args.checkpoint_path}")
    
    model = StripedHyena(config)
    load_checkpoint(model, args.checkpoint_path)
    model.to(device)
    model.eval()

    df = pd.read_csv(args.input_file, sep=r'\s+', header=None,
                     names=['variant_id', 'ref_seq', 'alt_seq', 'label', 'tissue_id'])
    df.dropna(subset=['ref_seq', 'alt_seq'], inplace=True)
    df = df[df['ref_seq'].str.len() > 0]
    df = df[df['alt_seq'].str.len() > 0]
    
    print(f"Processing {len(df):,} variants")
    
    all_ref_embeddings = []
    all_alt_embeddings = []
    
    for i in tqdm(range(0, len(df), args.batch_size), desc="Extracting"):
        batch_df = df.iloc[i:i+args.batch_size]
        ref_seqs = batch_df['ref_seq'].tolist()
        alt_seqs = batch_df['alt_seq'].tolist()
        
        ref_embeddings = get_evo2_embeddings(ref_seqs, model, tokenizer, device)
        all_ref_embeddings.append(ref_embeddings)
        
        alt_embeddings = get_evo2_embeddings(alt_seqs, model, tokenizer, device)
        all_alt_embeddings.append(alt_embeddings)
    
    ref_embeddings_array = np.vstack(all_ref_embeddings)
    alt_embeddings_array = np.vstack(all_alt_embeddings)
    diff_embeddings = alt_embeddings_array - ref_embeddings_array
    concatenated_features = np.concatenate(
        [ref_embeddings_array, alt_embeddings_array, diff_embeddings], axis=1
    )
    
    print(f"Embedding shape: {concatenated_features.shape}")
    
    np.savez_compressed(
        args.output_file,
        variant_ids=df['variant_id'].values,
        embeddings=concatenated_features,
        labels=df['label'].values,
        embedding_dim=ref_embeddings_array.shape[1]
    )
    
    print(f"Saved: {args.output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract Evo2 embeddings for variant effect prediction"
    )
    parser.add_argument("--config_path", required=True,
                        help="Path to Evo2 config file")
    parser.add_argument("--checkpoint_path", required=True,
                        help="Path to Evo2 checkpoint")
    parser.add_argument("--input_file", required=True,
                        help="Input TSV file with variants")
    parser.add_argument("--output_file", required=True,
                        help="Output NPZ file for embeddings")
    parser.add_argument("--batch_size", type=int, default=1,
                        help="Batch size for processing")
    
    args = parser.parse_args()
    main(args)

