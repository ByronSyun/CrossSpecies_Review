
import os
import re
import sys
import argparse
import logging
from typing import Tuple, Optional

import pandas as pd
from tqdm import tqdm
from pyfaidx import Fasta

logging.basicConfig(level=logging.INFO, format='%(asctime)s - [%(levelname)s] - %(message)s')

VARIANT_ID_PATTERNS = [
    re.compile(r"^(?P<chrom>[^:]+):(?P<pos>\d+)[_:](?P<ref>[ACGTacgt])[/|>](?P<alt>[ACGTacgt])$"),
]


def normalise_chrom(chrom: str) -> Tuple[str, str]:
    """Return (primary, fallback) chrom names after normalisation.
    - Strip leading 'chr'
    - Map M/chrM -> MT
    - Provide a fallback with 'chr' prefix for UCSC-style FASTA if needed
    """
    c = chrom
    if c.lower().startswith('chr'):
        c = c[3:]
    if c in ('M', 'Mt', 'm', 'MT'):
        c = 'MT'
    primary = c
    fallback = 'chr' + c  # in case FASTA uses UCSC names
    return primary, fallback


def parse_variant_id(variant_id: str) -> Optional[Tuple[str, int, str, str]]:
    for pat in VARIANT_ID_PATTERNS:
        m = pat.match(variant_id)
        if m:
            chrom = m.group('chrom')
            pos = int(m.group('pos'))
            ref = m.group('ref').upper()
            alt = m.group('alt').upper()
            return chrom, pos, ref, alt
    # Fallback: try very simple split (chrom:pos_ref/alt)
    try:
        chrom, rest = variant_id.split(':', 1)
        pos_str, refalt = rest.split('_', 1)
        pos = int(pos_str)
        ref, alt = re.split(r"[/|>]", refalt)
        return chrom, pos, ref.upper(), alt.upper()
    except Exception:
        return None


def fetch_slice(fasta: Fasta, chrom: str, start: int, end: int) -> Optional[str]:
    try:
        return fasta[str(chrom)][start - 1:end].seq.upper()
    except (KeyError, IndexError):
        return None


def extract_context_sequence(fasta: Fasta, chrom: str, pos: int, seq_len: int, strip_chr: bool) -> Optional[str]:
    """
    Build a window of exact length seq_len, placing the variant at index (seq_len//2 - 1)
    (0-based) in the returned sequence.
    Coordinates are 1-based. We set:
      start = pos - (seq_len//2) + 1
      end   = pos + (seq_len//2)
    which yields length = end - start + 1 = seq_len.
    """
    start = pos - (seq_len // 2) + 1
    end = pos + (seq_len // 2)
    if start < 1:
        return None
    if not strip_chr:
        return fetch_slice(fasta, chrom, start, end)
    # Try normalised primary then fallback with 'chr'
    primary, fallback = normalise_chrom(chrom)
    seq = fetch_slice(fasta, primary, start, end)
    if seq is None:
        seq = fetch_slice(fasta, fallback, start, end)
    return seq


def rebuild_sequences(input_tsv: str, fasta_path: str, output_tsv: str, seq_len: int, strip_chr: bool) -> None:
    logging.info("--- Rebuilding sequences for AlphaGenome ---")
    logging.info(f"Input TSV: {input_tsv}")
    logging.info(f"FASTA: {fasta_path}")
    logging.info(f"Output TSV: {output_tsv}")
    logging.info(f"Sequence length: {seq_len}")

    if not os.path.exists(input_tsv):
        logging.error(f"Input TSV not found: {input_tsv}")
        sys.exit(1)
    if not os.path.exists(fasta_path):
        logging.error(f"FASTA not found: {fasta_path}")
        sys.exit(1)

    df = pd.read_csv(input_tsv, sep='\t', header=None)
    # Expected columns in balanced TSV: variant_id, ref_sequence, alt_sequence, label, tissue_id
    if df.shape[1] < 5:
        logging.error("Input TSV must have at least 5 columns: variant_id, ref_sequence, alt_sequence, label, tissue_id")
        sys.exit(1)
    df.columns = ['variant_id', 'ref_sequence_old', 'alt_sequence_old', 'label', 'tissue_id']

    fasta = Fasta(fasta_path, key_function=lambda k: k.split(' ')[0])

    rebuilt_rows = []
    skipped_parse = 0
    skipped_ctx = 0
    skipped_ref_mismatch = 0
    processed = 0

    for _, row in tqdm(df.iterrows(), total=len(df), desc="Rebuilding"):
        processed += 1
        variant_id = str(row['variant_id'])
        parsed = parse_variant_id(variant_id)
        if not parsed:
            skipped_parse += 1
            continue
        chrom, pos, ref, alt = parsed
        if not (len(ref) == 1 and len(alt) == 1):
            skipped_parse += 1
            continue

        ctx = extract_context_sequence(fasta, chrom, pos, seq_len, strip_chr)
        if not ctx:
            skipped_ctx += 1
            continue

        center_idx = (seq_len // 2) - 1  # left-of-center index for even lengths
        if len(ctx) != seq_len:
            skipped_ctx += 1
            continue
        if ctx[center_idx] != ref:
            skipped_ref_mismatch += 1
            continue

        ref_seq = ctx
        alt_seq = ctx[:center_idx] + alt + ctx[center_idx + 1:]

        rebuilt_rows.append([
            variant_id,
            ref_seq,
            alt_seq,
            int(row['label']),
            int(row['tissue_id']) if pd.notnull(row['tissue_id']) else 15
        ])

    if not rebuilt_rows:
        logging.error("No records were rebuilt. Please check input formats and FASTA/assembly version.")
        logging.error(f"Processed: {processed:,}; skipped_parse={skipped_parse:,}; skipped_context={skipped_ctx:,}; skipped_ref_mismatch={skipped_ref_mismatch:,}")
        try:
            samp = df['variant_id'].head(5).tolist()
            logging.error(f"Sample variant_id head: {samp}")
        except Exception:
            pass
        sys.exit(1)

    out_df = pd.DataFrame(rebuilt_rows, columns=['variant_id', 'ref_sequence', 'alt_sequence', 'label', 'tissue_id'])
    os.makedirs(os.path.dirname(output_tsv), exist_ok=True)
    out_df.to_csv(output_tsv, sep='\t', index=False, header=False)

    logging.info("--- Rebuild complete ---")
    logging.info(f"Processed: {processed:,}")
    logging.info(f"Written: {len(out_df):,}")
    logging.info(f"Skipped (parse): {skipped_parse:,}")
    logging.info(f"Skipped (context): {skipped_ctx:,}")
    logging.info(f"Skipped (ref mismatch): {skipped_ref_mismatch:,}")
    logging.info(f"Output saved to: {output_tsv}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Rebuild ref/alt sequences from variant coordinates for AlphaGenome (RatGTEx).")
    parser.add_argument('--input_tsv', type=str, required=True,
                        help='Balanced TSV with columns: variant_id, ref_sequence, alt_sequence, label, tissue_id')
    parser.add_argument('--fasta_file', type=str, required=True,
                        help='Reference genome FASTA (e.g., Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa or rn7)')
    parser.add_argument('--seq_len', type=int, default=16384, help='Sequence length (must be supported by AlphaGenome)')
    parser.add_argument('--output_tsv', type=str, required=True,
                        help='Path to write rebuilt TSV (same 5 columns, new sequences)')
    parser.add_argument('--no_strip_chr', action='store_true', help='Do not normalise/strip chr prefix in chromosome names')
    args = parser.parse_args()

    rebuild_sequences(
        args.input_tsv,
        args.fasta_file,
        args.output_tsv,
        args.seq_len,
        strip_chr=not args.no_strip_chr
    )
