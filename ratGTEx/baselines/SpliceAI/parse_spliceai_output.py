import argparse
import sys

def parse_spliceai_vcf(vcf_file_path, output_tsv_path):
    """
    Parses a SpliceAI VCF output file to extract the maximum delta score for each variant.

    Args:
        vcf_file_path (str): Path to the input VCF file with SpliceAI annotations.
        output_tsv_path (str): Path to the output TSV file.
    """
    print(f"Reading SpliceAI VCF: {vcf_file_path}")
    print(f"Writing parsed scores TSV: {output_tsv_path}")

    with open(vcf_file_path, 'r') as vcf_file, open(output_tsv_path, 'w') as tsv_file:
        # Write header for the output file
        tsv_file.write("CHROM\tPOS\tID\tREF\tALT\tGENE\tDS_AG\tDS_AL\tDS_DG\tDS_DL\tMAX_DS\n")

        for line in vcf_file:
            # Skip VCF header lines
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            
            # Basic VCF columns
            chrom, pos, rs_id, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]
            info = fields[7]

            # Find the SpliceAI annotation in the INFO field
            if 'SpliceAI=' not in info:
                continue

            # Extract the SpliceAI string from semicolon-delimited INFO
            spliceai_info_str = ''
            for info_part in info.split(';'):
                if info_part.startswith('SpliceAI='):
                    spliceai_info_str = info_part.replace('SpliceAI=', '')
                    break
            
            if not spliceai_info_str:
                continue

            # Parse the SpliceAI annotation string: ALLELE|GENE|DS_AG|DS_AL|DS_DG|DS_DL|...
            parts = spliceai_info_str.split('|')
            
            # Ensure the format is as expected before trying to access indices
            if len(parts) < 6:
                continue
                
            gene = parts[1]
            try:
                ds_ag = float(parts[2])
                ds_al = float(parts[3])
                ds_dg = float(parts[4])
                ds_dl = float(parts[5])
            except (ValueError, IndexError):
                continue

            max_ds = max(ds_ag, ds_al, ds_dg, ds_dl)

            # Write the parsed information to the output TSV
            output_line = f"{chrom}\t{pos}\t{rs_id}\t{ref}\t{alt}\t{gene}\t{ds_ag}\t{ds_al}\t{ds_dg}\t{ds_dl}\t{max_ds}\n"
            tsv_file.write(output_line)

    print("Parsing complete.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse a SpliceAI VCF file to extract max delta scores.")
    parser.add_argument("-I", "--input_vcf", required=True, help="Path to the input VCF file annotated by SpliceAI.")
    parser.add_argument("-O", "--output_tsv", required=True, help="Path for the output TSV file.")
    
    args = parser.parse_args()

    parse_spliceai_vcf(args.input_vcf, args.output_tsv)
