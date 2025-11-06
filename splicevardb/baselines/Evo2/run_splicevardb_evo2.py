import argparse
import os
import pandas as pd
import torch
import yaml
from tqdm import tqdm
import torch.nn.functional as F

from vortex.model.model import StripedHyena
from vortex.model.tokenizer import CharLevelTokenizer
from vortex.model.utils import dotdict, load_checkpoint


def get_log_probs(seqs, model, tokenizer, device):
    """Return total log-probability for each sequence in seqs."""
    tokenized_seqs = [tokenizer.tokenize(s) for s in seqs]

    max_len = max(len(s) for s in tokenized_seqs)
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
        sequence_log_probs = log_probs[i, :len(input_ids)]
        total_log_prob = torch.gather(sequence_log_probs, 1, input_ids.unsqueeze(-1)).sum()
        results.append(total_log_prob.item())
    return results


def main(args):
    # 1) Load model and tokenizer
    config = dotdict(yaml.load(open(args.config_path), Loader=yaml.FullLoader))
    tokenizer = CharLevelTokenizer(config.vocab_size)
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    model = StripedHyena(config)
    load_checkpoint(model, args.checkpoint_path)
    model.to(device)
    model.eval()

    # 2) Load data
    df = pd.read_csv(args.input_file, sep='\t')

    # 3) Evaluate in batches
    results = []
    for i in tqdm(range(0, len(df), args.batch_size), desc="Evaluating variants"):
        batch_df = df.iloc[i:i+args.batch_size]
        ref_seqs = batch_df['sequence'].tolist()

        alt_seqs = []
        for _, row in batch_df.iterrows():
            seq = row['sequence']
            center_idx = len(seq) // 2
            alt_seq = list(seq)
            alt_seq[center_idx] = row['alt']
            alt_seqs.append("".join(alt_seq))

        all_seqs = ref_seqs + alt_seqs
        batch_log_probs = get_log_probs(all_seqs, model, tokenizer, device)
        logp_refs = batch_log_probs[:len(ref_seqs)]
        logp_alts = batch_log_probs[len(ref_seqs):]

        for j in range(len(batch_df)):
            delta_logp = logp_alts[j] - logp_refs[j]
            res = batch_df.iloc[j].to_dict()
            res['logp_ref'] = logp_refs[j]
            res['logp_alt'] = logp_alts[j]
            res['delta_logp'] = delta_logp
            results.append(res)

    # 4) Save results
    output_df = pd.DataFrame(results)
    cols = ['variant_id', 'chrom', 'pos_0based', 'ref', 'alt', 'label', 'logp_ref', 'logp_alt', 'delta_logp', 'sequence']
    output_df = output_df[cols]
    output_df.to_csv(args.output_file, sep='\t', index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Splice Variant Effect Prediction using Evo2 model.")
    parser.add_argument("--config_path", required=True, help="Path to model configuration file (e.g., evo2-7b-8k.yml).")
    parser.add_argument("--checkpoint_path", required=True, help="Path to the model checkpoint (.pt file).")
    parser.add_argument("--input_file", required=True, help="Path to the input TSV file prepared by `prepare_dataset_for_evo2.py`.")
    parser.add_argument("--output_file", required=True, help="Path to save the output TSV file with predictions.")
    parser.add_argument("--batch_size", type=int, default=4, help="Number of sequences to process in a batch.")
    args = parser.parse_args()
    main(args) 

