"""
AlphaGenome Predictor for RatGTEx Silver-Standard Benchmark (Binary Classification)
====================================================================================
This script is designed for the BIB review project to run zero-shot AlphaGenome 
predictions on the RatGTEx silver-standard binary classification dataset.

Core Logic:
  - Reads a TSV file containing reference and alternate sequences.
  - As AlphaGenome is not pre-trained on the Rattus norvegicus genome, this script
    utilises the full sequence context for predictions.
  - It performs batch predictions via the AlphaGenome API and calculates the max
    absolute difference between ref/alt scores as the prediction score for this
    binary classification task.
  - The script's structure and style are aligned with the AlphaGenome Colab script
    used for regression-based prediction.
"""

# =============================================================================
# --- Step 1: Library Imports ---
# =============================================================================
import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import json
import time
import warnings
import logging
from datetime import datetime
import glob

warnings.filterwarnings('ignore')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - [%(levelname)s] - %(message)s')

try:
    from alphagenome.models import dna_client
    logging.info("AlphaGenome library successfully imported.")
except ImportError:
    logging.error("AlphaGenome library not found. Please run: pip install alphagenome")
    exit()

# =============================================================================
# --- Step 2: Environment Setup & Path Configuration (Server-Optimized) ---
# =============================================================================
# This script is optimized for a server/Colab environment where the data
# and the script reside in the same directory.

try:
    from google.colab import drive
    logging.info("Mounting Google Drive for Colab environment...")
    drive.mount('/content/drive', force_remount=True)
    # --- Server/Colab Environment Paths ---
    BASE_DIR = '/content/drive/MyDrive/Colab Notebooks/ratGetx_classifier/'
    logging.info(f"Google Drive mounted. Base directory: {BASE_DIR}")
except (ImportError, ModuleNotFoundError):
    logging.warning("Not in a Colab environment. Assuming a flat directory structure.")
    # In a non-Colab server, you might run this from the directory containing the data.
    BASE_DIR = './' 

INPUT_FILE = os.path.join(BASE_DIR, 'ratgtex_silver_benchmark_balanced_len16384.tsv')
BATCH_SAVE_PATH = BASE_DIR

# =============================================================================
# --- Step 3: Global Configuration ---
# =============================================================================
# --- A. Authentication ---
# IMPORTANT: Please insert your API key here.
API_KEY = "AIzaSyBuH5MNAdx5M3wqR87b3mYhwigpGjefYro"

# --- B. Batch Processing Parameters ---
NUM_BATCHES = 3  # Split the dataset into 3 batches for safer processing
DELAY_SECONDS = 0.5  # Delay between individual predictions
BATCH_SIZE = 20  # Number of sequence pairs to process in a single API call

# --- C. Test Mode ---
TEST_MODE = False
TEST_SAMPLE_SIZE = 100 if TEST_MODE else None

logging.info("\n--- Prediction Task Configuration ---")
logging.info(f"   - Input File: {INPUT_FILE}")
logging.info(f"   - Batch Save Path: {BATCH_SAVE_PATH}")
logging.info(f"   - Number of Batches: {NUM_BATCHES}")
logging.info(f"   - Batch Size: {BATCH_SIZE}")
logging.info(f"   - Delay Between Predictions: {DELAY_SECONDS}s")
logging.info(f"   - Test Mode: {TEST_MODE}")
logging.info("-------------------------------------\n")

# =============================================================================
# --- Step 4: Batch Predictor Class Definition ---
# =============================================================================
class RatGTExBatchPredictor:
    """Handles batch prediction for the RatGTEx binary classification dataset with network interruption safety."""

    def __init__(self, api_key, batch_save_path):
        if not api_key or "YOUR_API_KEY" in api_key:
            raise ValueError("API_KEY is not set. Please provide a valid AlphaGenome API key.")
        self.api_key = api_key
        self.client = None
        self.batch_save_path = batch_save_path
        self._setup_client()

    def _setup_client(self):
        try:
            self.client = dna_client.create(api_key=self.api_key)
            logging.info("AlphaGenome API client initialized successfully.")
        except Exception as e:
            logging.error(f"Failed to initialize AlphaGenome client: {e}")
            raise

    def load_data(self, data_path, test_sample_size=None):
        """Loads and validates the silver-standard binary classification data."""
        logging.info(f"Loading data from {os.path.basename(data_path)}...")
        if not os.path.exists(data_path):
            raise FileNotFoundError(f"Input file not found: {data_path}")
        
        df = pd.read_csv(data_path, sep='\t', header=None)
        df.columns = ['variant_id', 'ref_sequence', 'alt_sequence', 'label', 'tissue_id']
        logging.info(f"   - Records: {len(df):,}")

        if test_sample_size is not None:
            logging.info(f"   - Test mode enabled. Using a random subset of {test_sample_size} samples.")
            df = df.sample(n=min(test_sample_size, len(df)), random_state=42).reset_index(drop=True)
        
        return df

    def create_batches(self, df_variants, num_batches):
        """Split DataFrame into specified number of batches."""
        batch_size = len(df_variants) // num_batches
        batches = []
        for i in range(num_batches):
            start_idx = i * batch_size
            if i == num_batches - 1:
                end_idx = len(df_variants)
            else:
                end_idx = (i + 1) * batch_size
            batch_df = df_variants.iloc[start_idx:end_idx].copy()
            batch_df.reset_index(drop=True, inplace=True)
            batches.append(batch_df)
        return batches

    def check_existing_batches(self):
        """Check which batches have already been processed."""
        existing_batches = []
        for i in range(NUM_BATCHES):
            batch_num = i + 1
            search_pattern = f"{self.batch_save_path}ratgtex_alphagenome_batch_{batch_num}_*.json"
            found_files = glob.glob(search_pattern)
            if found_files:
                existing_batches.append(batch_num)
        return existing_batches

    def analyze_batch(self, batch_df, batch_num, delay=0.5):
        """Process a single batch of variants."""
        results = []
        for i, row in tqdm(batch_df.iterrows(), total=len(batch_df), desc=f"Batch {batch_num}"):
            result = self.analyze_variant(row['ref_sequence'], row['alt_sequence'])
            result.update({
                'variant_id': row['variant_id'], 
                'label': row['label']
            })
            results.append(result)
            if i < len(batch_df) - 1:
                time.sleep(delay)
        return results

    def analyze_variant(self, ref_sequence, alt_sequence):
        """Analyze a single variant using AlphaGenome API."""
        result = {'status': 'error', 'alphagenome_score': -1.0, 'error': None}
        try:
            sequences_to_predict = [ref_sequence, alt_sequence]
            api_outputs = self.client.predict_sequences(
                sequences_to_predict,
                organism=dna_client.Organism.HOMO_SAPIENS,
                ontology_terms=["CL:0002551"],
                requested_outputs=[dna_client.OutputType.SPLICE_SITES]
            )
            
            ref_output = api_outputs[0]
            alt_output = api_outputs[1]
            
            score = 0.0
            if hasattr(ref_output, 'splice_sites') and hasattr(alt_output, 'splice_sites'):
                ref_scores = np.array(ref_output.splice_sites.values)
                alt_scores = np.array(alt_output.splice_sites.values)

                if ref_scores.size > 0 and alt_scores.size > 0:
                    diff = alt_scores - ref_scores
                    score = float(np.max(np.abs(diff))) if diff.size > 0 else 0.0
            
            result.update({'status': 'success', 'alphagenome_score': score})
        except Exception as e:
            result['error'] = str(e)
        return result

    def save_batch_results(self, results, batch_num):
        """Save batch results to JSON file."""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        batch_filename = f"ratgtex_alphagenome_batch_{batch_num}_{timestamp}.json"
        batch_filepath = os.path.join(self.batch_save_path, batch_filename)
        
        successful = [r for r in results if r['status'] == 'success']
        failed = [r for r in results if r['status'] == 'error']
        
        batch_data = {
            'metadata': {
                'batch_number': batch_num,
                'timestamp': timestamp,
                'total': len(results),
                'successful': len(successful),
                'failed': len(failed)
            },
            'results': results
        }
        
        with open(batch_filepath, 'w', encoding='utf-8') as f:
            json.dump(batch_data, f, indent=2)
        return batch_filepath

    def merge_all_batches(self):
        """Merge all completed batches into a single result file."""
        timestamp = datetime.now().strftime("%Y%m%d")
        all_results, total_successful, total_failed = [], 0, 0
        
        for i in range(NUM_BATCHES):
            batch_pattern = f"ratgtex_alphagenome_batch_{i+1}_{timestamp}"
            batch_files = [f for f in os.listdir(self.batch_save_path) if f.startswith(batch_pattern)]
            if batch_files:
                latest_batch_file = sorted(batch_files)[-1]
                batch_filepath = os.path.join(self.batch_save_path, latest_batch_file)
                with open(batch_filepath, 'r', encoding='utf-8') as f:
                    batch_data = json.load(f)
                all_results.extend(batch_data['results'])
                total_successful += batch_data['metadata']['successful']
                total_failed += batch_data['metadata']['failed']
        
        if all_results:
            final_timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            final_filename = f"ratgtex_alphagenome_complete_results_{final_timestamp}.json"
            final_filepath = os.path.join(self.batch_save_path, final_filename)
            
            final_data = {
                'metadata': {
                    'total_batches': NUM_BATCHES,
                    'final_timestamp': final_timestamp,
                    'total_variants': len(all_results),
                    'total_successful': total_successful,
                    'total_failed': total_failed
                },
                'results': all_results
            }
            
            with open(final_filepath, 'w', encoding='utf-8') as f:
                json.dump(final_data, f, indent=2)
            return final_filepath
        else:
            return None

# =============================================================================
# --- Step 5: Main Execution Block ---
# =============================================================================
def run_batch_analysis():
    """Main function to coordinate the entire prediction workflow with batch safety."""
    try:
        logging.info("[AlphaGenome] Initialising client and loading variants...")
        analyzer = RatGTExBatchPredictor(api_key=API_KEY, batch_save_path=BATCH_SAVE_PATH)
        df_variants = analyzer.load_data(INPUT_FILE, TEST_SAMPLE_SIZE)
        logging.info(f"[AlphaGenome] Loaded {len(df_variants)} variants from {INPUT_FILE}")

        existing_batches = analyzer.check_existing_batches()
        batches = analyzer.create_batches(df_variants, NUM_BATCHES)
        logging.info(f"[AlphaGenome] Starting analysis over {len(batches)} batches")

        for i, batch_df in enumerate(batches):
            batch_num = i + 1
            if batch_num in existing_batches:
                logging.info(f"[AlphaGenome] Skipping existing batch {batch_num}")
                continue
            logging.info(f"[AlphaGenome] Analyzing batch {batch_num} ({len(batch_df)} variants)...")
            batch_results = analyzer.analyze_batch(batch_df, batch_num, DELAY_SECONDS)
            analyzer.save_batch_results(batch_results, batch_num)
            logging.info(f"[AlphaGenome] Batch {batch_num} saved")

        merged = analyzer.merge_all_batches()
        if merged:
            logging.info(f"[AlphaGenome] Merged results: {os.path.basename(merged)}")
        else:
            logging.info("[AlphaGenome] No batches to merge")
    except Exception as e:
        logging.error(f"Pipeline error: {e}")

def merge_all_batches_fixed():
    """Re-merge utility function (if needed)."""
    all_results, total_successful, total_failed = [], 0, 0
    for i in range(NUM_BATCHES):
        batch_num = i + 1
        search_pattern = f"{BATCH_SAVE_PATH}ratgtex_alphagenome_batch_{batch_num}_*.json"
        batch_files = glob.glob(search_pattern)
        if batch_files:
            latest_batch_file = sorted(batch_files)[-1]
            with open(latest_batch_file, 'r', encoding='utf-8') as f:
                batch_data = json.load(f)
            all_results.extend(batch_data['results'])
            total_successful += batch_data['metadata']['successful']
            total_failed += batch_data['metadata']['failed']
    if all_results:
        final_timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        final_filename = f"ratgtex_alphagenome_COMPLETE_results_{final_timestamp}.json"
        final_filepath = os.path.join(BATCH_SAVE_PATH, final_filename)
        final_data = {
            'metadata': {
                'total_batches': NUM_BATCHES,
                'final_timestamp': final_timestamp,
                'total_variants': len(all_results),
                'total_successful': total_successful,
                'total_failed': total_failed
            },
            'results': all_results
        }
        with open(final_filepath, 'w', encoding='utf-8') as f:
            json.dump(final_data, f, indent=2)
        logging.info(f"Manual merge completed: {final_filename}")

if __name__ == '__main__':
    run_batch_analysis()
    
    # Uncomment the line below if you want to force a manual merge
    # merge_all_batches_fixed()
