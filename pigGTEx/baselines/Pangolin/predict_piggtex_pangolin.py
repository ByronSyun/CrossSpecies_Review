#!/usr/bin/env python3
"""
Run Pangolin predictions on PigGTEx data.

This script uses Pangolin to predict variant effects on the PigGTEx 
silver-standard dataset and saves results for evaluation.
"""

import os
import subprocess
import json
import pandas as pd
import logging
import argparse
from datetime import datetime
import tempfile

logging.basicConfig(level=logging.INFO, format='%(message)s')

PANGOLIN_ENV = "/mnt/userdata4/splicing/conda_envs/pangolin-env"
PIG_REFERENCE = "/mnt/userdata4/splicing/PigGTEx/reference_genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa"
PIG_ANNOTATION = "/mnt/userdata4/splicing/PigGTEx/reference_genome/Sscrofa11.1.annotation.db"

def check_environment():
    """Check server environment and validate dependencies."""
    if not os.path.exists("/mnt/userdata4/splicing/"):
        logging.error("This script only runs on the server environment!")
        return False, None, None, None
    
    logging.info("Server environment detected")
    
    missing_deps = []
    
    if not os.path.exists(PANGOLIN_ENV):
        missing_deps.append(f"Pangolin environment: {PANGOLIN_ENV}")
    if not os.path.exists(PIG_REFERENCE):
        missing_deps.append(f"Pig reference genome: {PIG_REFERENCE}")
    if not os.path.exists(PIG_ANNOTATION):
        missing_deps.append(f"Pig annotation database: {PIG_ANNOTATION}")
    
    if missing_deps:
        logging.error("Missing dependencies:")
        for dep in missing_deps:
            logging.error(f"  - {dep}")
        return False, None, None, None
    
    return True, PANGOLIN_ENV, PIG_REFERENCE, PIG_ANNOTATION

def run_pangolin_prediction(vcf_file, reference_genome, annotation_db, output_prefix, env_path=None):
    """Run Pangolin prediction on VCF file."""
    logging.info("=== Running Pangolin Prediction ===")
    logging.info(f"Input VCF: {vcf_file}")
    logging.info(f"Reference: {reference_genome}")
    logging.info(f"Annotation: {annotation_db}")
    logging.info(f"Output prefix: {output_prefix}")
    
    cmd = ["pangolin", vcf_file, reference_genome, annotation_db, output_prefix]
    
    if env_path:
        cmd = ["conda", "run", "-p", env_path] + cmd
    
    logging.info(f"Command: {' '.join(cmd)}")
    
    try:
        # Pangolin on ~26k variants can exceed 1 hour; extend timeout to 6 hours
        result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=21600)
        
        logging.info("Pangolin stdout:")
        logging.info(result.stdout)
        
        if result.stderr:
            logging.warning("Pangolin stderr:")
            logging.warning(result.stderr)
        
        output_file = output_prefix + ".vcf"
        if os.path.exists(output_file):
            logging.info(f"✅ Pangolin prediction completed: {output_file}")
            return output_file
        else:
            logging.error(f"Expected output file not found: {output_file}")
            return None
            
    except subprocess.TimeoutExpired:
        logging.error("Pangolin prediction timed out (1 hour)")
        return None
    except subprocess.CalledProcessError as e:
        logging.error(f"Pangolin prediction failed with exit code {e.returncode}")
        logging.error(f"Stdout: {e.stdout}")
        logging.error(f"Stderr: {e.stderr}")
        return None
    except Exception as e:
        logging.error(f"Unexpected error running Pangolin: {e}")
        return None

def parse_pangolin_output(pangolin_vcf, mapping_csv, output_json):
    """Parse Pangolin VCF output and merge with original labels."""
    logging.info("=== Processing Pangolin Output ===")
    
    try:
        logging.info(f"Reading Pangolin output: {pangolin_vcf}")
        
        pangolin_results = []
        with open(pangolin_vcf, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) >= 8:
                    chrom, pos, variant_id, ref, alt = fields[0:5]
                    info = fields[7] if len(fields) > 7 else ""
                    
                    pangolin_score = 0.0
                    
                    if '|' in info:
                        parts = info.split('|')
                        if len(parts) >= 3:
                            try:
                                increase_part = parts[1].split(':')[-1] if ':' in parts[1] else parts[1]
                                decrease_part = parts[2].split(':')[-1] if ':' in parts[2] else parts[2]
                                
                                increase_score = float(increase_part) if increase_part.replace('.','').replace('-','').isdigit() else 0.0
                                decrease_score = float(decrease_part) if decrease_part.replace('.','').replace('-','').isdigit() else 0.0
                                
                                pangolin_score = max(abs(increase_score), abs(decrease_score))
                            except (ValueError, IndexError):
                                pangolin_score = 0.0
                    
                    pangolin_results.append({
                        'CHROM': chrom,
                        'POS': int(pos),
                        'variant_id': variant_id,
                        'REF': ref,
                        'ALT': alt,
                        'pangolin_score': pangolin_score,
                        'pangolin_raw_info': info
                    })
        
        logging.info(f"Parsed {len(pangolin_results)} Pangolin predictions")
        
        logging.info(f"Reading original labels: {mapping_csv}")
        labels_df = pd.read_csv(mapping_csv)
        
        # Ensure data types match for merge
        labels_df['POS'] = labels_df['POS'].astype(int)
        labels_df['CHROM'] = labels_df['CHROM'].astype(str)
        labels_df['REF'] = labels_df['REF'].astype(str)
        labels_df['ALT'] = labels_df['ALT'].astype(str)
        
        pangolin_df = pd.DataFrame(pangolin_results)
        pangolin_df['CHROM'] = pangolin_df['CHROM'].astype(str)
        pangolin_df['REF'] = pangolin_df['REF'].astype(str)
        pangolin_df['ALT'] = pangolin_df['ALT'].astype(str)
        
        merged_df = pd.merge(
            pangolin_df,
            labels_df[['CHROM', 'POS', 'REF', 'ALT', 'label', 'variant_id']],
            on=['CHROM', 'POS', 'REF', 'ALT'],
            how='inner',
            suffixes=('_pred', '_orig')
        )
        
        if merged_df.empty:
            logging.error("No variants matched between Pangolin output and original labels")
            return False
        
        logging.info(f"Successfully merged {len(merged_df)} variants")
        
        results = []
        for _, row in merged_df.iterrows():
            results.append({
                'variant_id': row['variant_id_orig'],
                'chrom': row['CHROM'],
                'pos': int(row['POS']),
                'ref': row['REF'],
                'alt': row['ALT'],
                'label': int(row['label']),
                'pangolin_score': float(row['pangolin_score']),
                'pangolin_raw_info': row['pangolin_raw_info']
            })
        
        output_data = {
            'metadata': {
                'timestamp': datetime.now().strftime("%Y%m%d_%H%M%S"),
                'model': 'Pangolin',
                'dataset': 'PigGTEx_Silver_Standard',
                'total_variants': len(results),
                'positive_variants': sum(1 for r in results if r['label'] == 1),
                'negative_variants': sum(1 for r in results if r['label'] == 0),
                'input_vcf': pangolin_vcf,
                'mapping_csv': mapping_csv
            },
            'results': results
        }
        
        os.makedirs(os.path.dirname(output_json), exist_ok=True)
        with open(output_json, 'w') as f:
            json.dump(output_data, f, indent=2)
        
        logging.info(f"✅ Results saved to: {output_json}")
        return True
        
    except Exception as e:
        logging.error(f"Error processing Pangolin output: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description="Run Pangolin predictions on PigGTEx data")
    
    parser.add_argument('--input_vcf', type=str, required=True, help='Input VCF file')
    parser.add_argument('--mapping_csv', type=str, required=True, help='Mapping CSV file with labels')
    parser.add_argument('--output_json', type=str, required=True, help='Output JSON file')
    parser.add_argument('--temp_dir', type=str, help='Temporary directory for intermediate files')
    
    args = parser.parse_args()
    
    deps_ok, env_path, ref_genome, annotation_db = check_environment()
    
    if not deps_ok:
        logging.error("❌ Environment check failed. Please set up dependencies.")
        return 1
    
    if not os.path.exists(args.input_vcf):
        logging.error(f"Input VCF file not found: {args.input_vcf}")
        return 1
    
    if not os.path.exists(args.mapping_csv):
        logging.error(f"Mapping CSV file not found: {args.mapping_csv}")
        return 1
    
    if args.temp_dir:
        os.makedirs(args.temp_dir, exist_ok=True)
        temp_prefix = os.path.join(args.temp_dir, "pangolin_piggtex")
    else:
        temp_dir = tempfile.mkdtemp(prefix="pangolin_piggtex_")
        temp_prefix = os.path.join(temp_dir, "pangolin_output")
    
    try:
        logging.info("Starting Pangolin pipeline...")
        
        pangolin_output = run_pangolin_prediction(args.input_vcf, ref_genome, annotation_db, temp_prefix, env_path)
        
        if not pangolin_output:
            logging.error("❌ Pangolin prediction failed")
            return 1
        
        success = parse_pangolin_output(pangolin_output, args.mapping_csv, args.output_json)
        
        if success:
            logging.info("🎉 Pangolin pipeline completed successfully!")
            logging.info(f"Results saved to: {args.output_json}")
            return 0
        else:
            logging.error("❌ Results processing failed")
            return 1
            
    except Exception as e:
        logging.error(f"❌ Pipeline failed: {e}")
        return 1

if __name__ == '__main__':
    exit(main())
