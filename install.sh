#!/bin/bash
set -euo pipefail

wget https://developer.download.nvidia.com/compute/cuda/repos/wsl-ubuntu/x86_64/cuda-keyring_1.1-1_all.deb
sudo dpkg -i cuda-keyring_1.1-1_all.deb
rm cuda-keyring_1.1-1_all.deb
rm cuda-keyring_1.1-1_all.deb.1

sudo apt update
sudo apt install -y cuda-toolkit-13-0

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

REF_DIR="$ROOT_DIR/references/raw"
mkdir -p "$REF_DIR"

function download_reference() {
    local BASE_URL="$1"
    local FILE="$2"

    local FILE_UNZIPPED="$REF_DIR/${FILE%.gz}"

    if [ ! -f "$REF_DIR/$FILE" ]; then
        echo "Downloading '$FILE'"
        wget -O "$REF_DIR/$FILE" "$BASE_URL/$FILE"
    fi

    # Unzipping (Keep originals)
    if [ ! -f "$FILE_UNZIPPED" ]; then
        gunzip -k "$REF_DIR/$FILE"
    fi
}

# ==========================================
# 1. Download and install 10x Genomics tools
# ==========================================

TOOL_DIR="$ROOT_DIR/.tools"
mkdir -p "$TOOL_DIR"

install_10x_tool() {
    local PROGRAM="$1" 

    if ls "$TOOL_DIR"/${PROGRAM}-*/${PROGRAM} 1> /dev/null 2>&1; then
        echo "$PROGRAM already installed. Skipping download."
        local PROGRAM_DIR=$(ls -d "$TOOL_DIR"/${PROGRAM}-* | head -n 1)
    else
        echo -ne "\033[1;32mPaste $PROGRAM download URL and press Enter: \033[0m"
        read DOWNLOAD_URL
        DOWNLOAD_URL=$(echo "$DOWNLOAD_URL" | tr -d '"' | tr -d "'")

        local ARCHIVE_NAME="${PROGRAM}.tar.gz"
        wget -O "$ARCHIVE_NAME" "$DOWNLOAD_URL"
        tar -xzvf "$ARCHIVE_NAME" -C "$TOOL_DIR"
        rm "$ARCHIVE_NAME"
        local PROGRAM_DIR=$(ls -d "$TOOL_DIR"/${PROGRAM}-* | head -n 1)
    fi
    export PATH="$PROGRAM_DIR:$PATH"

    if ! command -v "$PROGRAM" &> /dev/null; then
        echo "Error: $PROGRAM path setup failed."
        exit 1
    fi
    echo "$PROGRAM is ready: $(which $PROGRAM)"
}

install_10x_tool "cellranger"
install_10x_tool "cellranger-arc"

# ==========================================
# 2. Building STAR index
# ==========================================
function build_star_index() {
    local BASE_URL="$1"
    local FASTA_FILE="$2"
    local GTF_FILE="$3"
    local OUTPUT_REFERENCE="$4"

    echo "Building STAR index: $OUTPUT_REFERENCE"
    echo "Checking Reference Files..."

    # Downloading fasta and gtf file if needed
    download_reference "$BASE_URL" "$FASTA_FILE"
    download_reference "$BASE_URL" "$GTF_FILE" 

    FASTA_UNZIPPED="$REF_DIR/${FASTA_FILE%.gz}"
    GTF_UNZIPPED="$REF_DIR/${GTF_FILE%.gz}"

    sed -i 's/Curated Genomic/Curated_Genomic/g' "$GTF_UNZIPPED"

    BULK_REFERENCE_DIR="$ROOT_DIR/references/bulk"
    mkdir -p "$BULK_REFERENCE_DIR"
    cd "$BULK_REFERENCE_DIR" || exit 1

    if [ -d "$OUTPUT_REFERENCE" ]; then
        echo "Index '$OUTPUT_REFERENCE' already exists"
    else
        echo "Running STAR"
        STAR \
            --runThreadN "$N_THREADS" \
            --runMode genomeGenerate \
            --genomeDir "$BULK_REFERENCE_DIR/$OUTPUT_REFERENCE" \
            --genomeFastaFiles "$FASTA_UNZIPPED" \
            --sjdbGTFfile "$GTF_UNZIPPED" 
    fi
}

# Homo sapiens (human) GRCh38.p14
build_star_index \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14" \
    "GCF_000001405.40_GRCh38.p14_genomic.fna.gz" \
    "GCF_000001405.40_GRCh38.p14_genomic.gtf.gz" \
    "GRCh38"

# Mus musculus (house mouse) GRCm39
build_star_index \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39" \
    "GCF_000001635.27_GRCm39_genomic.fna.gz" \
    "GCF_000001635.27_GRCm39_genomic.gtf.gz" \
    "GRCm39"

# Callithrix jacchus (white-tufted-ear marmoset) calJac240_pri
build_star_index \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/049/354/715/GCF_049354715.1_calJac240_pri" \
    "GCF_049354715.1_calJac240_pri_genomic.fna.gz" \
    "GCF_049354715.1_calJac240_pri_genomic.gtf.gz" \
    "calJac240_pri"

# ==========================================
# 3. Building cellranger References
# ==========================================
function build_cellranger_ref() {
    local BASE_URL="$1"
    local FASTA_FILE="$2"
    local GTF_FILE="$3"
    local OUTPUT_REFERENCE="$4"

    echo "Building Cellranger Reference: $OUTPUT_REFERENCE"
    echo "Checking Reference Files..."

    # Downloading fasta and gtf file if needed
    download_reference "$BASE_URL" "$FASTA_FILE"
    download_reference "$BASE_URL" "$GTF_FILE" 

    FASTA_UNZIPPED="$REF_DIR/${FASTA_FILE%.gz}"
    GTF_UNZIPPED="$REF_DIR/${GTF_FILE%.gz}"

    sed -i 's/Curated Genomic/Curated_Genomic/g' "$GTF_UNZIPPED"

    local SC_REF_ROOT="$ROOT_DIR/references/sc"
    mkdir -p "$SC_REF_ROOT"
    cd "$SC_REF_ROOT" || exit 1

    if [ -d "$OUTPUT_REFERENCE" ]; then
        echo "Reference '$OUTPUT_REFERENCE' already exists, skipping cellranger mkref"
    else
        echo "Running cellranger mkref..."
        cellranger telemetry disable
        cellranger mkref \
            --genome="$OUTPUT_REFERENCE" \
            --fasta="$FASTA_UNZIPPED" \
            --genes="$GTF_UNZIPPED" \
            --nthreads="$N_THREADS"
    fi
}

# Homo sapiens (human) GRCh38.p14
build_cellranger_ref \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14" \
    "GCF_000001405.40_GRCh38.p14_genomic.fna.gz" \
    "GCF_000001405.40_GRCh38.p14_genomic.gtf.gz" \
    "GRCh38"
    
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
# 5. Downloading Ensembl GTF 
# ==========================================
# Homo Sapiens
download_reference \
    "https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens" \
    "Homo_sapiens.GRCh38.115.gtf.gz" 

# Mus Musculus
download_reference \
    "https://ftp.ensembl.org/pub/release-115/gtf/mus_musculus" \
    "Mus_musculus.GRCm39.115.gtf.gz" 
