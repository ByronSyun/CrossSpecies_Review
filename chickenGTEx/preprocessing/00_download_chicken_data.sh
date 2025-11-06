#!/bin/bash
set -euo pipefail

RAW_DATA_DIR="./chicken_raw_data"
REFERENCE_DIR="./reference_genome"
mkdir -p "${RAW_DATA_DIR}" "${REFERENCE_DIR}"

GENOME_URL="https://ftp.ensembl.org/pub/release-102/fasta/gallus_gallus/dna/Gallus_gallus.GRCg6a.dna.toplevel.fa.gz"
GENOME_FILE="${REFERENCE_DIR}/Gallus_gallus.GRCg6a.dna.toplevel.fa.gz"

if [ ! -f "$GENOME_FILE" ]; then
    wget -c "${GENOME_URL}" -O "${GENOME_FILE}"
fi

SQTL_URL="https://ngdc.cncb.ac.cn/chickengtex/static/download/chicken_gtex_v0.significant_cis_sqtl.txt.tar.gz"
SQTL_FILE_NAME="chicken_gtex_v0.significant_cis_sqtl.txt.tar.gz"
SQTL_FILE="${RAW_DATA_DIR}/${SQTL_FILE_NAME}"

if [ ! -f "$SQTL_FILE" ]; then
    if ! wget -c "${SQTL_URL}" -O "${SQTL_FILE}"; then
        echo "Download failed. Visit: https://ngdc.cncb.ac.cn/chickengtex/download"
        exit 1
    fi
fi

GITHUB_DIR="${RAW_DATA_DIR}/ChickenGTEx_pilot_phase"
if [ -d "$GITHUB_DIR" ]; then
    cd "$GITHUB_DIR" && git pull && cd - > /dev/null
else
    git clone https://github.com/guandailu/ChickenGTEx_pilot_phase.git "$GITHUB_DIR"
fi

VCF_URL="https://ftp.ensembl.org/pub/release-102/variation/vcf/gallus_gallus/gallus_gallus.vcf.gz"
VCF_FILE="${RAW_DATA_DIR}/gallus_gallus_variation.vcf.gz"
VCF_ALT_URL="https://ftp.ncbi.nlm.nih.gov/snp/organisms/chicken_9031/VCF/00-All.vcf.gz"
VCF_ALT_FILE="${RAW_DATA_DIR}/chicken_dbSNP_all.vcf.gz"

if [ ! -f "$VCF_FILE" ] && [ ! -f "$VCF_ALT_FILE" ]; then
    if ! wget -c "${VCF_URL}" -O "${VCF_FILE}" 2>/dev/null; then
        if ! wget -c "${VCF_ALT_URL}" -O "${VCF_ALT_FILE}" 2>/dev/null; then
            echo "VCF download failed. Manual download required."
            exit 1
        fi
    fi
fi

SQTL_UNPACK_DIR="${RAW_DATA_DIR}/sQTL_summary_stats"
if [ -f "$SQTL_FILE" ] && [ ! -d "${SQTL_UNPACK_DIR}" ]; then
    mkdir -p "${SQTL_UNPACK_DIR}"
    tar -zxf "${SQTL_FILE}" -C "${SQTL_UNPACK_DIR}" --strip-components=1
fi

echo "ChickenGTEx data download completed."