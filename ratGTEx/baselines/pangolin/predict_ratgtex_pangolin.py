#!/usr/bin/env python3
"""
Run Pangolin predictions on RatGTEx data.

This script uses Pangolin to predict variant effects on the RatGTEx 
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

logging.basicConfig(level=logging.INFO, format='%(asctime)s - [%(levelname)s] - %(message)s')

# Server paths (script only runs on server)
PANGOLIN_ENV = "/mnt/userdata4/splicing/conda_envs/pangolin-env"
RAT_REFERENCE = "/mnt/userdata4/splicing/ratgetx/reference_genome/rn6/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa"
RAT_ANNOTATION = "/mnt/userdata4/splicing/ratgetx/reference_genome/ratNor6.annotation.db"
RATGTEX_DATA = "/mnt/userdata4/splicing/ratgetx/processed_data/ratgtex_silver_benchmark_balanced.tsv"

def check_environment():
    """
    Check server environment and validate dependencies.
    
    Note: This script only runs on server with official Pangolin installation:
    - PyTorch installed first
    - Dependencies: pyvcf, gffutils, biopython, pandas, pyfastx
    - Pangolin installed from source (git clone + pip install .)
    """
    # Verify we're on the server
    if not os.path.exists("/mnt/userdata4/splicing/"):
        logging.error("This script only runs on the server environment!")
        return False, None, None, None
    
    logging.info("Server environment detected")
    
    # Check dependencies
    missing_deps = []
    
    if not os.path.exists(PANGOLIN_ENV):
        missing_deps.append(f"Pangolin environment: {PANGOLIN_ENV}")
    if not os.path.exists(RAT_REFERENCE):
        missing_deps.append(f"Rat reference genome: {RAT_REFERENCE}")
    if not os.path.exists(RAT_ANNOTATION):
        missing_deps.append(f"Rat annotation database: {RAT_ANNOTATION}")
    
    if missing_deps:
        logging.error("Missing dependencies:")
        for dep in missing_deps:
            logging.error(f"  - {dep}")
        return False, None, None, None
    
    return True, PANGOLIN_ENV, RAT_REFERENCE, RAT_ANNOTATION

def run_pangolin_prediction(vcf_file, reference_genome, annotation_db, output_prefix, env_path=None):
    """
    Run Pangolin prediction on VCF file.
    
    Args:
        vcf_file (str): Input VCF file
        reference_genome (str): Path to reference genome FASTA
        annotation_db (str): Path to annotation database
        output_prefix (str): Output file prefix
        env_path (str): Path to conda environment (if on server)
    
    Returns:
        str: Path to output file
    """
    logging.info("=== Running Pangolin Prediction ===")
    logging.info(f"Input VCF: {vcf_file}")
    logging.info(f"Reference: {reference_genome}")
    logging.info(f"Annotation: {annotation_db}")
    logging.info(f"Output prefix: {output_prefix}")
    
    # Pangolin command
    cmd = [
        "pangolin",
        vcf_file,
        reference_genome,
        annotation_db,
        output_prefix
    ]
    
    # If on server, activate conda environment
    if env_path:
        # Use conda run to execute in the specific environment
        cmd = ["conda", "run", "-p", env_path] + cmd
    
    logging.info(f"Command: {' '.join(cmd)}")
    
    # Run Pangolin
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
            timeout=3600  # 1 hour timeout
        )
        
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
    """
    Parse Pangolin VCF output and merge with original labels.
    
    Args:
        pangolin_vcf (str): Pangolin output VCF file
        mapping_csv (str): CSV file with original labels
        output_json (str): Output JSON file path
    
    Returns:
        bool: Success status
    """
    logging.info("=== Processing Pangolin Output ===")
    
    try:
        # Load Pangolin predictions
        logging.info(f"Reading Pangolin output: {pangolin_vcf}")
        
        # Parse VCF file (simplified parser for Pangolin output)
        pangolin_results = []
        with open(pangolin_vcf, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) >= 8:
                    chrom, pos, variant_id, ref, alt = fields[0:5]
                    info = fields[7] if len(fields) > 7 else ""
                    
                    # Parse Pangolin scores from INFO field
                    # Pangolin format: gene|pos:largest_increase|pos:largest_decrease|
                    pangolin_score = 0.0  # Default
                    
                    # Simple parsing - extract numeric values from info field
                    # This is a simplified parser, may need adjustment based on actual Pangolin output
                    if '|' in info:
                        parts = info.split('|')
                        if len(parts) >= 3:
                            try:
                                # Try to extract the maximum absolute score
                                increase_part = parts[1].split(':')[-1] if ':' in parts[1] else parts[1]
                                decrease_part = parts[2].split(':')[-1] if ':' in parts[2] else parts[2]
                                
                                increase_score = float(increase_part) if increase_part.replace('.','').replace('-','').isdigit() else 0.0
                                decrease_score = float(decrease_part) if decrease_part.replace('.','').replace('-','').isdigit() else 0.0
                                
                                # Take maximum absolute value as the variant effect score
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
        
        # Load original labels
        logging.info(f"Reading original labels: {mapping_csv}")
        labels_df = pd.read_csv(mapping_csv)
        
        # Merge predictions with labels
        pangolin_df = pd.DataFrame(pangolin_results)
        merged_df = pd.merge(
            pangolin_df,
            labels_df[['CHROM', 'POS', 'REF', 'ALT', 'label', 'variant_id', 'tissue_id']],
            on=['CHROM', 'POS', 'REF', 'ALT'],
            how='inner',
            suffixes=('_pred', '_orig')
        )
        
        if merged_df.empty:
            logging.error("No variants matched between Pangolin output and original labels")
            return False
        
        logging.info(f"Successfully merged {len(merged_df)} variants")
        
        # Convert to final format
        results = []
        for _, row in merged_df.iterrows():
            results.append({
                'variant_id': row['variant_id_orig'],
                'chrom': row['CHROM'],
                'pos': int(row['POS']),
                'ref': row['REF'],
                'alt': row['ALT'],
                'label': int(row['label']),
                'tissue_id': row['tissue_id'],
                'pangolin_score': float(row['pangolin_score']),
                'pangolin_raw_info': row['pangolin_raw_info']
            })
        
        # Save results
        output_data = {
            'metadata': {
                'timestamp': datetime.now().strftime("%Y%m%d_%H%M%S"),
                'model': 'Pangolin',
                'dataset': 'RatGTEx_Silver_Standard',
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
    parser = argparse.ArgumentParser(
        description="Run Pangolin predictions on RatGTEx data"
    )
    
    # Default paths (server only)
    default_vcf = "/mnt/userdata4/splicing/pangolin/data/ratgtex_for_pangolin.vcf"
    default_mapping = "/mnt/userdata4/splicing/pangolin/data/ratgtex_for_pangolin_mapping.csv"
    default_output = "/mnt/userdata4/splicing/pangolin/results/pangolin_ratgtex_results.json"
    
    parser.add_argument(
        '--input_vcf',
        type=str,
        default=default_vcf,
        help=f'Input VCF file (default: {default_vcf})'
    )
    
    parser.add_argument(
        '--mapping_csv',
        type=str,
        default=default_mapping,
        help=f'Mapping CSV file with labels (default: {default_mapping})'
    )
    
    parser.add_argument(
        '--output_json',
        type=str,
        default=default_output,
        help=f'Output JSON file (default: {default_output})'
    )
    
    parser.add_argument(
        '--temp_dir',
        type=str,
        help='Temporary directory for intermediate files'
    )
    
    args = parser.parse_args()
    
    # Check environment and dependencies
    deps_ok, env_path, ref_genome, annotation_db = check_environment()
    
    if not deps_ok:
        logging.error("❌ Environment check failed. Please set up dependencies.")
        return 1
    
    # Validate input files
    if not os.path.exists(args.input_vcf):
        logging.error(f"Input VCF file not found: {args.input_vcf}")
        return 1
    
    if not os.path.exists(args.mapping_csv):
        logging.error(f"Mapping CSV file not found: {args.mapping_csv}")
        return 1
    
    # Set up temporary directory
    if args.temp_dir:
        os.makedirs(args.temp_dir, exist_ok=True)
        temp_prefix = os.path.join(args.temp_dir, "pangolin_ratgtex")
    else:
        temp_dir = tempfile.mkdtemp(prefix="pangolin_ratgtex_")
        temp_prefix = os.path.join(temp_dir, "pangolin_output")
    
    try:
        # Run Pangolin prediction
        logging.info("Starting Pangolin pipeline...")
        
        pangolin_output = run_pangolin_prediction(
            args.input_vcf,
            ref_genome,
            annotation_db,
            temp_prefix,
            env_path
        )
        
        if not pangolin_output:
            logging.error("❌ Pangolin prediction failed")
            return 1
        
        # Process results
        success = parse_pangolin_output(
            pangolin_output,
            args.mapping_csv,
            args.output_json
        )
        
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
