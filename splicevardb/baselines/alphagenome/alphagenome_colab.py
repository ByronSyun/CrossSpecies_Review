# -*- coding: utf-8 -*-
"""
AlphaGenome SpliceVarDB batch processing script (network interruption safe)
This script supports batch processing and avoids data loss due to network issues.
"""

import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import json
import time
from datetime import datetime
import warnings
import glob
warnings.filterwarnings('ignore')

from alphagenome.data import genome
from alphagenome.models import dna_client
from alphagenome.visualization import plot_components

from google.colab import drive
try:
    drive.mount('/content/drive', force_remount=True)
except Exception as e:
    raise

# IMPORTANT: Replace with your API key
API_KEY = "AIzaSyBuH5MNAdx5M3wqR87b3mYhwigpGjefYro"

# Use benchmark variant list produced by preprocessing
DATA_PATH = "/content/drive/MyDrive/Colab Notebooks/splicevardb_filter_benchmark.tsv"

NUM_BATCHES = 3
DELAY_SECONDS = 0.5
BATCH_SAVE_PATH = "/content/drive/MyDrive/Colab Notebooks/"

TEST_MODE = False
TEST_SAMPLE_SIZE = 30 if TEST_MODE else None

class BatchSpliceAnalyzer:
    """Batch splice analysis with interruption safety"""

    def __init__(self, api_key, batch_save_path):
        self.api_key = api_key
        self.client = None
        self.batch_save_path = batch_save_path
        self._setup_client()

    def _setup_client(self):
        try:
            self.client = dna_client.create(api_key=self.api_key)
        except Exception as e:
            raise

    def load_data(self, data_path, test_sample_size=None):
        if not os.path.exists(data_path):
            raise FileNotFoundError(f"File not found: {data_path}")

        df = pd.read_csv(data_path, sep='\t', low_memory=False)
        
        valid_variants = []
        for _, row in df.iterrows():
            coords = self._parse_coords(row['hg38'])
            if coords:
                chrom, pos, ref, alt = coords
                if len(ref) == 1 and len(alt) == 1:
                    valid_variants.append({
                        'chrom': chrom, 'pos': pos, 'ref': ref, 'alt': alt,
                        'classification': row['classification'],
                        'location': row['location'],
                        'variant_id': f"{chrom}:{pos}:{ref}>{alt}"
                    })
        
        df_clean = pd.DataFrame(valid_variants)
        
        if test_sample_size and len(df_clean) > test_sample_size:
            df_clean = df_clean.sample(n=test_sample_size, random_state=42)
        
        return df_clean

    def _parse_coords(self, hg38_string):
        try:
            if pd.isna(hg38_string): return None
            parts = str(hg38_string).strip('"').split('-')
            if len(parts) != 4: return None
            chrom, pos_str, ref, alt = parts
            if not chrom.startswith('chr'): chrom = f'chr{chrom}'
            pos = int(pos_str)
            if ref in 'ATGC' and alt in 'ATGC': return (chrom, pos, ref, alt)
            return None
        except: return None

    def analyze_variant(self, chrom, pos, ref, alt):
        try:
            variant = genome.Variant(chromosome=chrom, position=pos, reference_bases=ref, alternate_bases=alt)
            interval_length = 1048576
            start_pos = max(1, pos - interval_length // 2)
            end_pos = start_pos + interval_length
            interval = genome.Interval(chromosome=chrom, start=start_pos, end=end_pos)
            outputs = self.client.predict_variant(
                variant=variant, interval=interval,
                ontology_terms=["CL:0002551"],
                requested_outputs=[
                    dna_client.OutputType.SPLICE_SITES,
                    dna_client.OutputType.SPLICE_JUNCTIONS,
                ]
            )
            return {'status': 'success', 'analysis': self._analyze_outputs(outputs)}
        except Exception as e:
            return {'status': 'error', 'error': str(e)}

    def _analyze_outputs(self, outputs):
        analysis = {'reference': {}, 'alternate': {}, 'differential': {}}
        try:
            if hasattr(outputs.reference, 'splice_junctions'):
                analysis['reference']['n_junctions'] = len(outputs.reference.splice_junctions.junctions)
                analysis['alternate']['n_junctions'] = len(outputs.alternate.splice_junctions.junctions)
            if hasattr(outputs.reference, 'splice_sites'):
                ref_vals = np.array(outputs.reference.splice_sites.values)
                alt_vals = np.array(outputs.alternate.splice_sites.values)
                if ref_vals.size > 0:
                    diff = alt_vals - ref_vals
                    analysis['differential']['max_abs_change'] = float(np.max(np.abs(diff)))
                    analysis['differential']['n_significant'] = int(np.sum(np.abs(diff) > 0.1))
        except Exception as e: analysis['error'] = str(e)
        return analysis

    def create_batches(self, df_variants, num_batches):
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
        existing_batches = []
        for i in range(NUM_BATCHES):
            batch_num = i + 1
            search_pattern = f"{self.batch_save_path}alphagenome_batch_{batch_num}_*.json"
            found_files = glob.glob(search_pattern)
            if found_files:
                existing_batches.append(batch_num)
        return existing_batches

    def analyze_batch(self, batch_df, batch_num, delay=0.5):
        results = []
        for i, row in tqdm(batch_df.iterrows(), total=len(batch_df), desc=f"Batch {batch_num}"):
            result = self.analyze_variant(chrom=row['chrom'], pos=row['pos'], ref=row['ref'], alt=row['alt'])
            result.update({'variant_id': row['variant_id'], 'classification': row['classification']})
            results.append(result)
            if i < len(batch_df) - 1:
                time.sleep(delay)
        return results

    def save_batch_results(self, results, batch_num):
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        batch_filename = f"alphagenome_batch_{batch_num}_{timestamp}.json"
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
        timestamp = datetime.now().strftime("%Y%m%d")
        all_results, total_successful, total_failed = [], 0, 0
        for i in range(NUM_BATCHES):
            batch_pattern = f"alphagenome_batch_{i+1}_{timestamp}"
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
            final_filename = f"alphagenome_complete_results_{final_timestamp}.json"
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

def run_batch_analysis():
    try:
        print("[AlphaGenome] Initialising client and loading variants...")
        analyzer = BatchSpliceAnalyzer(api_key=API_KEY, batch_save_path=BATCH_SAVE_PATH)
        df_variants = analyzer.load_data(DATA_PATH, TEST_SAMPLE_SIZE)
        print(f"[AlphaGenome] Loaded {len(df_variants)} variants from {DATA_PATH}")

        existing_batches = analyzer.check_existing_batches()
        batches = analyzer.create_batches(df_variants, NUM_BATCHES)
        print(f"[AlphaGenome] Starting analysis over {len(batches)} batches")

        for i, batch_df in enumerate(batches):
            batch_num = i + 1
            if batch_num in existing_batches:
                print(f"[AlphaGenome] Skipping existing batch {batch_num}")
                continue
            print(f"[AlphaGenome] Analyzing batch {batch_num} ({len(batch_df)} variants)...")
            batch_results = analyzer.analyze_batch(batch_df, batch_num, DELAY_SECONDS)
            analyzer.save_batch_results(batch_results, batch_num)
            print(f"[AlphaGenome] Batch {batch_num} saved")

        merged = analyzer.merge_all_batches()
        if merged:
            print(f"[AlphaGenome] Merged results: {os.path.basename(merged)}")
        else:
            print("[AlphaGenome] No batches to merge")
    except Exception as e:
        print(f"Pipeline error: {e}")

run_batch_analysis()

# Re-merge utility (if needed)
def merge_all_batches_fixed():
    all_results, total_successful, total_failed = [], 0, 0
    for i in range(NUM_BATCHES):
        batch_num = i + 1
        search_pattern = f"{BATCH_SAVE_PATH}alphagenome_batch_{batch_num}_*.json"
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
        final_filename = f"alphagenome_COMPLETE_results_{final_timestamp}.json"
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

merge_all_batches_fixed()
