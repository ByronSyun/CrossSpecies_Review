
import pandas as pd
import argparse
import os
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - [%(levelname)s] - %(message)s')

def create_balanced_dataset(positive_file, negative_file, output_file):
    """
    Creates a class-balanced and shuffled dataset for the RatGTEx benchmark.
    
    It reads all positive samples, and randomly samples an equal number of negative
    samples from the negative pool. The two are then combined and shuffled.
    """
    logging.info("--- Starting Balanced RatGTEx Benchmark Creation ---")

    # --- 1. Load Data ---
    try:
        logging.info(f"Loading positive samples from: {positive_file}")
        pos_df = pd.read_csv(positive_file, sep='\t', header=None)
        num_positives = len(pos_df)
        logging.info(f"Found {num_positives:,} positive samples.")

        logging.info(f"Loading negative samples from: {negative_file}")
        neg_df = pd.read_csv(negative_file, sep='\t', header=None)
        logging.info(f"Loaded a pool of {len(neg_df):,} negative samples.")

    except FileNotFoundError as e:
        logging.error(f"Error: Input file not found. {e}")
        logging.error("Please ensure you have run the 'ratgetx/run.sh' script to generate the source files.")
        return

    # --- 2. Sample Negative Data ---
    if len(neg_df) < num_positives:
        logging.error(f"Negative sample pool ({len(neg_df)}) is smaller than the positive set ({num_positives}). Cannot create a balanced set.")
        return
        
    logging.info(f"Randomly sampling {num_positives:,} negative samples to match positives...")
    neg_df_sampled = neg_df.sample(n=num_positives, random_state=42)
    logging.info("Sampling complete.")

    # --- 3. Combine and Shuffle ---
    logging.info("Combining positive and sampled negative dataframes...")
    combined_df = pd.concat([pos_df, neg_df_sampled], ignore_index=True)
    
    logging.info("Shuffling the combined dataset to ensure random order...")
    final_df = combined_df.sample(frac=1, random_state=42).reset_index(drop=True)
    logging.info("Shuffling complete.")
    
    # --- 4. Save Final Dataset ---
    # Ensure column names are not written to the file
    final_df.columns = [0, 1, 2, 3, 4] 
    
    # Ensure parent directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    final_df.to_csv(output_file, sep='\t', index=False, header=False)
    
    logging.info("-" * 50)
    logging.info("✅ Success! Balanced benchmark created.")
    logging.info(f"   - Total Samples: {len(final_df):,}")
    logging.info(f"   - Positive Samples: {len(final_df[final_df[3] == 1]):,}")
    logging.info(f"   - Negative Samples: {len(final_df[final_df[3] == 0]):,}")
    logging.info(f"   - Output file saved to: {output_file}")
    logging.info("----------------------------------------------------")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Create a class-balanced version of the RatGTEx silver-standard benchmark.")
    
    BASE_DIR = '/Users/byronsun/Desktop/AS_复现模型/'
    
    parser.add_argument('--positive_file', type=str, 
                        default=os.path.join(BASE_DIR, 'ratgetx/processed_data/rat_positive_samples_with_tissue.tsv'),
                        help='Path to the TSV file containing positive samples.')
                        
    parser.add_argument('--negative_file', type=str, 
                        default=os.path.join(BASE_DIR, 'ratgetx/processed_data/rat_negative_subset_for_training.tsv'),
                        help='Path to the TSV file containing the pool of negative samples.')
                        
    parser.add_argument('--output_file', type=str, 
                        default=os.path.join(BASE_DIR, 'BIB_review/data/processed_data/ratGTEx/ratgtex_silver_benchmark_balanced.tsv'),
                        help='Path to save the final balanced and shuffled benchmark file.')
                        
    args = parser.parse_args()
    
    create_balanced_dataset(args.positive_file, args.negative_file, args.output_file)
