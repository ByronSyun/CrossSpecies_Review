#!/usr/bin/env python3
import argparse
import pandas as pd
import torch
import yaml
from tqdm import tqdm
import torch.nn.functional as F
from vortex.model.model import StripedHyena
from vortex.model.tokenizer import CharLevelTokenizer
from vortex.model.utils import dotdict, load_checkpoint

def get_log_probs(seqs, model, tokenizer, device):
    if not seqs:
        return []
        
    tokenized_seqs = [tokenizer.tokenize(s) for s in seqs]
    max_len = max(len(s) for s in tokenized_seqs) if tokenized_seqs else 0
    if max_len == 0:
        return [0.0] * len(seqs)

    padded_tokens = []
    for s in tokenized_seqs:
        pad_len = max_len - len(s)
        padded_tokens.append(F.pad(torch.tensor(s, dtype=torch.long), (0, pad_len), 'constant', 0))
    
    input_tensor = torch.stack(padded_tokens).to(device)

    with torch.no_grad():
        logits, _ = model(input_tensor)

    log_probs = F.log_softmax(logits, dim=-1)
    
    results = []
    for i in range(len(tokenized_seqs)):
        input_ids = torch.tensor(tokenized_seqs[i], dtype=torch.long, device=device)[1:]
        if input_ids.nelement() == 0:
            results.append(0.0)
            continue
            
        sequence_log_probs = log_probs[i, :len(input_ids)]
        total_log_prob = torch.gather(sequence_log_probs, 1, input_ids.unsqueeze(-1)).sum()
        results.append(total_log_prob.item())
        
    return results

def main(args):
    print(f"Loading config from: {args.config_path}")
    config = dotdict(yaml.load(open(args.config_path), Loader=yaml.FullLoader))
    
    print("Initializing tokenizer...")
    tokenizer = CharLevelTokenizer(config.vocab_size)

    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    print("Initializing model...")
    model = StripedHyena(config)
    
    print(f"Loading checkpoint from: {args.checkpoint_path}")
    load_checkpoint(model, args.checkpoint_path)
    model.to(device)
    model.eval()

    print(f"Loading dataset from: {args.input_file}")
    df = pd.read_csv(args.input_file, sep=r'\s+', header=None,
                     names=['variant_id', 'ref_seq', 'alt_seq', 'label', 'tissue_id'])
    df.dropna(subset=['ref_seq', 'alt_seq'], inplace=True)
    df = df[df['ref_seq'].str.len() > 0]
    df = df[df['alt_seq'].str.len() > 0]
    print(f"Found {len(df)} valid variants to process.")

    all_results = []
    for i in tqdm(range(0, len(df), args.batch_size), desc="Predicting"):
        batch_df = df.iloc[i:i+args.batch_size]
        
        ref_seqs = batch_df['ref_seq'].tolist()
        alt_seqs = batch_df['alt_seq'].tolist()
        all_seqs = ref_seqs + alt_seqs
        
        batch_log_probs = get_log_probs(all_seqs, model, tokenizer, device)
        
        ref_log_probs = batch_log_probs[:len(ref_seqs)]
        alt_log_probs = batch_log_probs[len(ref_seqs):]
        
        for j in range(len(ref_log_probs)):
            delta_log_prob = alt_log_probs[j] - ref_log_probs[j]
            all_results.append({
                'variant_id': batch_df['variant_id'].iloc[j],
                'logp_ref': ref_log_probs[j],
                'logp_alt': alt_log_probs[j],
                'delta_logp': delta_log_prob,
                'label': batch_df['label'].iloc[j]
            })

    results_df = pd.DataFrame(all_results)
    
    print(f"Saving predictions to: {args.output_file}")
    results_df.to_csv(args.output_file, sep='\t', index=False)
    print("Prediction complete!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config_path", required=True)
    parser.add_argument("--checkpoint_path", required=True)
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--output_file", required=True)
    parser.add_argument("--batch_size", type=int, default=1)
    
    args = parser.parse_args()
    main(args)

