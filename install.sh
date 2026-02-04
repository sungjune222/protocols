#!/bin/bash

set -euo pipefail

sudo apt install pip
pip install -U --no-deps "git+ssh://git@github.com/broadinstitute/CellBender.git@v0.3.2"

echo "Start: Cell Ranger Installation & Reference Building"

# ==========================================
# 1. Download and install Cell Ranger
# ==========================================

TOOL_DIR="$ROOT_DIR/.tools"

mkdir -p "$TOOL_DIR"

if ls "$TOOL_DIR"/cellranger-*/cellranger 1> /dev/null 2>&1; then
    echo "Cell Ranger already installed. Skipping download."
    CELLRANGER_PATH=$(ls -d "$TOOL_DIR"/cellranger-* | head -n 1)
else
    read -p "Paste Cell Ranger download URL and press Enter: " DOWNLOAD_URL
    DOWNLOAD_URL=$(echo "$DOWNLOAD_URL" | tr -d '"' | tr -d "'")
    wget -O cellranger.tar.gz "$DOWNLOAD_URL"
    tar -xzvf cellranger.tar.gz -C "$TOOL_DIR"
    rm cellranger.tar.gz
    CELLRANGER_PATH=$(ls -d "$TOOL_DIR"/cellranger-* | head -n 1)
fi
export PATH="$CELLRANGER_PATH:$PATH"

if ! command -v cellranger &> /dev/null; then
    echo "Error: Cell Ranger path setup failed."
    exit 1
fi
echo "Cell Ranger is ready: $(which cellranger)"

# ==========================================
# 2. T2T-CHM13 v2.0 Reference Builder
# ==========================================
BASE_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0"
FASTA_FILE="GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz"
GTF_FILE="GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf.gz"

REF_DIR="$ROOT_DIR/references/raw"
OUTPUT_GENOME="T2T_CHM13_v2_0"

mkdir -p "$REF_DIR"

echo "Checking Reference Files..."

if [ ! -f "$REF_DIR/$FASTA_FILE" ]; then
    echo "Downloading FASTA..."
    wget -O "$REF_DIR/$FASTA_FILE" "$BASE_URL/$FASTA_FILE"
fi

if [ ! -f "$REF_DIR/$GTF_FILE" ]; then
    echo "Downloading GTF..."
    wget -O "$REF_DIR/$GTF_FILE" "$BASE_URL/$GTF_FILE"
fi

echo "Unzipping (Keep originals)..."
if [ ! -f "$REF_DIR/${FASTA_FILE%.gz}" ]; then
    gunzip -k "$REF_DIR/$FASTA_FILE"
fi
if [ ! -f "$REF_DIR/${GTF_FILE%.gz}" ]; then
    gunzip -k "$REF_DIR/$GTF_FILE"
fi

FASTA_UNZIPPED="$REF_DIR/${FASTA_FILE%.gz}"
GTF_UNZIPPED="$REF_DIR/${GTF_FILE%.gz}"

sed -i 's/Curated Genomic/Curated_Genomic/g' "$GTF_UNZIPPED"

echo "Running cellranger mkref..."
cd references

cellranger telemetry disable
cellranger mkref \
    --genome="$OUTPUT_GENOME" \
    --fasta="$FASTA_UNZIPPED" \
    --genes="$GTF_UNZIPPED" \
    --nthreads="$N_THREADS"