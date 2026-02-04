#!/bin/bash
ENV_FILE="./.env"
if [ -f "$ENV_FILE" ]; then
    source "$ENV_FILE"
else
    echo "Error: Cannot find ($ENV_FILE)"
    exit 1
fi

TOOL_PATH="$ROOT_DIR/.tools"
CELLRANGER_PATH=$(ls -d "$TOOL_PATH"/cellranger-* | head -n 1)
export PATH="$CELLRANGER_PATH:$PATH"

# ==========================================
# Please set your parameters here
# ==========================================
PROJECT_ID="PRJNA544731"
# ==========================================

SRA_DATA_DIR="$ROOT_DIR/$SRA_DATA"
FASTQ_DATA_DIR="$ROOT_DIR/$FASTQ_DATA"
CELLRANGER_DIR="$ROOT_DIR/$CELLRANGER"

mkdir -p $SRA_DATA_DIR
mkdir -p $FASTQ_DATA_DIR
mkdir -p $CELLRANGER_DIR

esearch -db sra -query $PROJECT_ID | efetch -format runinfo > runinfo.csv
cut -d',' -f1 runinfo.csv | tail -n +2 > srr_list.txt

TOTAL_COUNT=$(wc -l < srr_list.txt)
echo "Found $TOTAL_COUNT samples"

count=1
while read SRR_ID; do
    echo "Processing $count / $TOTAL_COUNT : $SRR_ID (Mode: $MODE)"

    if [ -f "$SRA_DATA_DIR/$SRR_ID/$SRR_ID.sra" ]; then
        echo "$SRR_ID.sra already exists. Skipping download."
    else
        echo "Downloading .sra..."
        prefetch $SRR_ID --output-directory $SRA_DATA_DIR
    fi

    echo "Converting to FastQ..."
    fasterq-dump "$SRA_DATA_DIR/$SRR_ID/$SRR_ID.sra" \
        --split-files \
        --threads $N_THREADS \
        --outdir "$FASTQ_DATA_DIR/$SRR_ID" \
        --temp "$FASTQ_DATA_DIR/temp"

    mv "$FASTQ_DATA_DIR/$SRR_ID/${SRR_ID}_1.fastq" "$FASTQ_DATA_DIR/$SRR_ID/${SRR_ID}_S1_L001_R1_001.fastq"
    mv "$FASTQ_DATA_DIR/$SRR_ID/${SRR_ID}_2.fastq" "$FASTQ_DATA_DIR/$SRR_ID/${SRR_ID}_S1_L001_R2_001.fastq"
    cd "$CELLRANGER_DATA_DIR"
    cellranger count --id="$SRR_ID" \
        --transcriptome="$ROOT_DIR/references/T2T_CHM13_v2_0" \
        --fastqs="$FASTQ_DATA_DIR/$SRR_ID" \
        --sample="$SRR_ID" \
        --localcores=$N_THREADS \
        --create-bam false
    cd - > /dev/null

    ((count++))

done < srr_list.txt

