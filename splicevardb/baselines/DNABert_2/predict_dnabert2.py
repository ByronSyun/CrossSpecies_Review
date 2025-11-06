import pandas as pd
import torch
from transformers import AutoTokenizer, BertModel
import numpy as np
from tqdm import tqdm
import argparse
import os

def get_dnabert2_embeddings(sequences, model, tokenizer, device, batch_size=32):
    """
    Generates embeddings for a list of DNA sequences using DNABERT-2.
    
    Args:
        sequences (list): A list of DNA sequences (strings).
        model: The pre-trained DNABERT-2 model.
        tokenizer: The pre-trained DNABERT-2 tokenizer.
        device: The torch device (e.g., 'cuda' or 'cpu').
        batch_size (int): The number of sequences to process in each batch.
        
    Returns:
        np.array: A numpy array of embeddings.
    """
    all_embeddings = []
    
    # Process sequences in batches
    for i in tqdm(range(0, len(sequences), batch_size), desc="Generating Embeddings"):
        batch_sequences = sequences[i:i+batch_size]
        
        inputs = tokenizer(batch_sequences, padding="longest", return_tensors='pt', truncation=True, max_length=512)
        inputs = {key: val.to(device) for key, val in inputs.items()}
        
        with torch.no_grad():
            hidden_states = model(**inputs)[0]
        
        # Mean pooling of the last hidden states
        batch_embeddings = torch.mean(hidden_states, dim=1)
        all_embeddings.append(batch_embeddings.cpu().numpy())
        
    return np.vstack(all_embeddings)

def create_mutant_sequences(df):
    """
    Creates mutant sequences from a dataframe containing WT sequences and mutation info.
    It assumes the mutation is at the center of the 'sequence'.
    """
    wt_sequences = df['sequence'].tolist()
    mut_sequences = []
    for index, row in tqdm(df.iterrows(), total=df.shape[0], desc="Creating Mutant Sequences"):
        wt_seq = row['sequence']
        ref = row['ref']
        alt = row['alt']
        
        # Assume sequence length is odd and mutation is in the center
        center_index = len(wt_seq) // 2
        
        # Sanity check
        if wt_seq[center_index].upper() != ref.upper():
            print(f"Warning: Mismatch at index {index}. WT sequence center '{wt_seq[center_index]}' != ref '{ref}'. Using ref for replacement.")
        
        # Create mutant sequence
        mut_seq_list = list(wt_seq)
        mut_seq_list[center_index] = alt
        mut_sequences.append("".join(mut_seq_list))
        
    return wt_sequences, mut_sequences

def main(args):
    """
    Main function to generate differential embeddings from SpliceVarDB data.
    """
    # Setup device
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    # Load DNABERT-2 model and tokenizer
    print("Loading DNABERT-2 model and tokenizer...")
    tokenizer = AutoTokenizer.from_pretrained("zhihan1996/DNABERT-2-117M", trust_remote_code=True)
    model = BertModel.from_pretrained("zhihan1996/DNABERT-2-117M", trust_remote_code=True)
    model.to(device)
    model.eval()
    print("Model and tokenizer loaded.")

    # Load data
    print(f"Loading data from {args.input_file}...")
    df = pd.read_csv(args.input_file, sep='\t')
    
    # Create WT and MT sequences
    wt_sequences, mt_sequences = create_mutant_sequences(df)
    
    # Generate embeddings for WT and MT sequences
    print("Generating embeddings for Wild-Type sequences...")
    wt_embeddings = get_dnabert2_embeddings(wt_sequences, model, tokenizer, device, args.batch_size)
    
    print("Generating embeddings for Mutant-Type sequences...")
    mt_embeddings = get_dnabert2_embeddings(mt_sequences, model, tokenizer, device, args.batch_size)
    
    # Calculate differential embeddings
    diff_embeddings = mt_embeddings - wt_embeddings
    
    # Prepare labels
    labels = df['label'].apply(lambda x: 1 if x == 'splice-altering' else 0).values
    variant_ids = df['variant_id'].values

    # Save results
    print(f"Saving results to {args.output_file}...")
    np.savez_compressed(
        args.output_file, 
        variant_ids=variant_ids,
        embeddings=diff_embeddings,
        labels=labels
    )
    print("Done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate DNABERT-2 differential embeddings for splice variants.")
    parser.add_argument(
        "--input_file", 
        type=str, 
        default="../Evo2Splicevardb/splicevardb_data/evo2_splicevardb_dataset_dedup.tsv",
        help="Path to the input TSV file from SpliceVarDB."
    )
    parser.add_argument(
        "--output_file", 
        type=str, 
        default="dnabert2_predictions.npz",
        help="Path to save the output NPZ file."
    )
    parser.add_argument(
        "--batch_size", 
        type=int, 
        default=32,
        help="Batch size for generating embeddings."
    )
    
    args = parser.parse_args()
    main(args) 