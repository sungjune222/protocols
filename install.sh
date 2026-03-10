#!/bin/bash
set -euo pipefail

sudo apt update
sudo apt install -y docker.io
curl -fsSL https://nvidia.github.io/libnvidia-container/gpgkey | sudo gpg --dearmor -o /usr/share/keyrings/nvidia-container-toolkit-keyring.gpg \
  && curl -s -L https://nvidia.github.io/libnvidia-container/stable/deb/nvidia-container-toolkit.list | \
    sed 's#deb https://#deb [signed-by=/usr/share/keyrings/nvidia-container-toolkit-keyring.gpg] https://#' | \
    sudo tee /etc/apt/sources.list.d/nvidia-container-toolkit.list

sudo apt install -y aria2
sudo apt install -y nvidia-container-toolkit
sudo nvidia-ctk runtime configure --runtime=docker
sudo systemctl restart docker

sudo apt install -y acl
sudo usermod -aG docker $USER
newgrp docker
sudo chown -R $USER:$USER .
sudo setfacl -R -d -m u:$USER:rwx .
sudo setfacl -R -m u:$USER:rwx .

ENV_FILE="./.env"
if [ -f "$ENV_FILE" ]; then
    source "$ENV_FILE"
else
    echo "Error: Cannot find ($ENV_FILE)"
    exit 1
fi

# ==========================================
# 1. Download and install Cell Ranger
# ==========================================

TOOL_DIR="$ROOT_DIR/.tools"

mkdir -p "$TOOL_DIR"

if ls "$TOOL_DIR"/cellranger-*/cellranger 1> /dev/null 2>&1; then
    echo "Cell Ranger already installed. Skipping download."
    CELLRANGER_PATH=$(ls -d "$TOOL_DIR"/cellranger-* | head -n 1)
else
    echo -ne "\033[1;32mPaste Cell Ranger download URL and press Enter: \033[0m"
    read DOWNLOAD_URL
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
# 2. Reference Builder
# ==========================================
function build_cellranger_ref() {
    local BASE_URL="$1"
    local FASTA_FILE="$2"
    local GTF_FILE="$3"
    local OUTPUT_GENOME="$4"

    local REF_DIR="$ROOT_DIR/references/raw"
    local FASTA_UNZIPPED="$REF_DIR/${FASTA_FILE%.gz}"
    local GTF_UNZIPPED="$REF_DIR/${GTF_FILE%.gz}"

    echo " Building Reference: $OUTPUT_GENOME"
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
    if [ ! -f "$FASTA_UNZIPPED" ]; then
        gunzip -k "$REF_DIR/$FASTA_FILE"
    fi
    if [ ! -f "$GTF_UNZIPPED" ]; then
        gunzip -k "$REF_DIR/$GTF_FILE"
    fi

    sed -i 's/Curated Genomic/Curated_Genomic/g' "$GTF_UNZIPPED"

    mkdir -p "$ROOT_DIR/references"
    cd "$ROOT_DIR/references" || exit 1

    if [ -d "$OUTPUT_GENOME" ]; then
        echo "Reference '$OUTPUT_GENOME' already exists. Skipping mkref."
    else
        echo "Running cellranger mkref..."
        cellranger telemetry disable
        cellranger mkref \
            --genome="$OUTPUT_GENOME" \
            --fasta="$FASTA_UNZIPPED" \
            --genes="$GTF_UNZIPPED" \
            --nthreads="$N_THREADS"
    fi
}

# ==========================================
# 3. Building References
# ==========================================

# Homo sapiens (human) T2T-CHM13 v2.0
build_cellranger_ref \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0" \
    "GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz" \
    "GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf.gz" \
    "T2T_CHM13_v2_0"

# Mus musculus (house mouse) GRCm39
build_cellranger_ref \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39" \
    "GCF_000001635.27_GRCm39_genomic.fna.gz" \
    "GCF_000001635.27_GRCm39_genomic.gtf.gz" \
    "GRCm39"

# Callithrix jacchus (white-tufted-ear marmoset) calJac240_pri
build_cellranger_ref \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/049/354/715/GCF_049354715.1_calJac240_pri" \
    "GCF_049354715.1_calJac240_pri_genomic.fna.gz" \
    "GCF_049354715.1_calJac240_pri_genomic.gtf.gz" \
    "calJac240_pri"

# ==========================================
# 4. Download Homo Sapiens Ensembl GTF for later use
# ==========================================
BASE_URL="https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens"
GTF_FILE="Homo_sapiens.GRCh38.115.gtf.gz"

REF_DIR="$ROOT_DIR/references/raw"

if [ ! -f "$REF_DIR/$GTF_FILE" ]; then
    echo "Downloading GTF..."
    wget -O "$REF_DIR/$GTF_FILE" "$BASE_URL/$GTF_FILE"
fi

if [ ! -f "$REF_DIR/${GTF_FILE%.gz}" ]; then
    gunzip -k "$REF_DIR/$GTF_FILE"
fi

GTF_UNZIPPED="$REF_DIR/${GTF_FILE%.gz}"

# ==========================================
# 5. Download Mus Musculus Ensembl GTF for later use
# ==========================================
BASE_URL="https://ftp.ensembl.org/pub/release-115/gtf/mus_musculus"
GTF_FILE="Mus_musculus.GRCm39.115.gtf.gz"

REF_DIR="$ROOT_DIR/references/raw"

if [ ! -f "$REF_DIR/$GTF_FILE" ]; then
    echo "Downloading GTF..."
    wget -O "$REF_DIR/$GTF_FILE" "$BASE_URL/$GTF_FILE"
fi

if [ ! -f "$REF_DIR/${GTF_FILE%.gz}" ]; then
    gunzip -k "$REF_DIR/$GTF_FILE"
fi

GTF_UNZIPPED="$REF_DIR/${GTF_FILE%.gz}"
