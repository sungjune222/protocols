#!/bin/bash
set -euo pipefail

ENV_FILE="./.env"
if [ -f "$ENV_FILE" ]; then
    source "$ENV_FILE"
else
    echo "Error: Cannot find ($ENV_FILE)"
    exit 1
fi
PIXI_EXEC=$(command -v pixi)

# ==========================================
# Please set your parameters here
# ==========================================
PROJECT_ID="PRJNA1108209"
REFERENCE_GENOME="calJac240_pri"
# ==========================================

TOOL_PATH="$ROOT_DIR/.tools"
CELLRANGER_PATH=$(ls -d "$TOOL_PATH"/cellranger-* | head -n 1)
export PATH="$CELLRANGER_PATH:$PATH"
BAMTOFASTQ="$CELLRANGER_PATH/lib/bin/bamtofastq"

FASTQ_DATA_DIR="$ROOT_DIR/$FASTQ_DATA"

CELLRANGER_DIR="$ROOT_DIR/$CELLRANGER"
CELLRANGER_PROJECT_DIR="$CELLRANGER_DIR/$PROJECT_ID"

CELLBENDER_DIR="$ROOT_DIR/$CELLBENDER"
CELLBENDER_PROJECT_DIR="$CELLBENDER_DIR/$PROJECT_ID"

mkdir -p "$FASTQ_DATA_DIR"
mkdir -p "$CELLRANGER_PROJECT_DIR"
mkdir -p "$CELLBENDER_PROJECT_DIR"

echo "=== [PART 1] Detecting Project Type and Download ==="

SRA_PROJECT_DIR="$ROOT_DIR/$SRA_DATA/$PROJECT_ID"
SRA_META_DIR="$SRA_PROJECT_DIR/_meta"
mkdir -p "$SRA_META_DIR"

RUNINFO_CSV="$SRA_META_DIR/runinfo.csv"
SRR_LIST="$SRA_META_DIR/srr_list.txt"
SRR_URL_MAP="$SRA_META_DIR/srr_url_map.txt"
DOWNLOAD_INPUT_LIST="$SRA_META_DIR/download_input_list.txt"

esearch -db sra -query "$PROJECT_ID" | efetch -format runinfo > "$RUNINFO_CSV"
cut -d',' -f1 "$RUNINFO_CSV" | tail -n +2 > "$SRR_LIST"
awk -F',' 'NR>1 {print $1 "\t" $10}' "$RUNINFO_CSV" > "$SRR_URL_MAP"

TOTAL_COUNT=$(wc -l < "$SRR_LIST")
echo "Found $TOTAL_COUNT samples in $PROJECT_ID"

echo "=== [PART 1] Starting High-Speed Download (aria2c) ==="

rm -f "$DOWNLOAD_INPUT_LIST"
: > "$DOWNLOAD_INPUT_LIST"

while IFS=$'\t' read -r SRR_ID DOWNLOAD_URL; do
    if [ -n "$DOWNLOAD_URL" ]; then
        if [ -f "$SRA_PROJECT_DIR/$SRR_ID.sra" ]; then
            echo "Skipping $SRR_ID (Already exists)"
        else
            echo "$DOWNLOAD_URL" >> "$DOWNLOAD_INPUT_LIST"
            echo "  out=$SRR_ID.sra" >> "$DOWNLOAD_INPUT_LIST"
            echo "  dir=$SRA_PROJECT_DIR" >> "$DOWNLOAD_INPUT_LIST"
        fi
    fi
done < "$SRR_URL_MAP"

if [ -s "$DOWNLOAD_INPUT_LIST" ]; then
    aria2c -i "$DOWNLOAD_INPUT_LIST" \
        -x 16 -s 16 -j 4 \
        --file-allocation=none \
        --summary-interval=0
else
    echo "All files seem to be downloaded."
fi


echo "=== [PART 2] Processing Pipeline ==="
count=1
while IFS= read -r SRR_ID; do
    echo "[Process] $count / $TOTAL_COUNT : $SRR_ID"

    SRA_FILE="$SRA_PROJECT_DIR/$SRR_ID.sra"

    if [ ! -f "$SRA_FILE" ]; then
        echo "  Error: SRA file not found at $SRA_FILE. Skipping."
        ((count++))
        continue
    fi

    OUTPUT_PATH="$CELLRANGER_PROJECT_DIR/$SRR_ID"

    if [ -d "$OUTPUT_PATH" ]; then
        echo "  -> Output directory exists. Skipping analysis."
    else
        echo "  -> Converting SRA to FastQ..."
        mkdir -p "$FASTQ_DATA_DIR/$SRR_ID"
        
        fasterq-dump "$SRA_FILE" \
            --split-files \
            --threads "$N_THREADS" \
            --outdir "$FASTQ_DATA_DIR/$SRR_ID" \
            --temp "$FASTQ_DATA_DIR/temp" > /dev/null 2>&1

        mv "$FASTQ_DATA_DIR/$SRR_ID/${SRR_ID}_1.fastq" "$FASTQ_DATA_DIR/$SRR_ID/${SRR_ID}_S1_L001_R1_001.fastq"
        mv "$FASTQ_DATA_DIR/$SRR_ID/${SRR_ID}_2.fastq" "$FASTQ_DATA_DIR/$SRR_ID/${SRR_ID}_S1_L001_R2_001.fastq"

        echo "  -> Running Cell Ranger..."
        cd "$CELLRANGER_PROJECT_DIR" || exit
        
        cellranger count --id="$SRR_ID" \
            --transcriptome="$ROOT_DIR/references/$REFERENCE_GENOME" \
            --fastqs="$FASTQ_DATA_DIR/$SRR_ID" \
            --sample="$SRR_ID" \
            --localcores="$N_THREADS" \
            --create-bam false > "$SRR_ID.log" 2>&1
        
        EXIT_CODE=$?
        if [ -d "$SRR_ID" ]; then
            mv "$SRR_ID.log" "$SRR_ID/"
        fi
        cd - > /dev/null

        if [ $EXIT_CODE -eq 0 ]; then
            echo "  -> Cell Ranger finished successfully. Cleaning up FASTQ..."
            rm -rf "$FASTQ_DATA_DIR/$SRR_ID"
        else
            echo "  -> Warning: Cell Ranger failed. Keeping FASTQ"
        fi
    fi

    ((count++))
done < "$SRR_LIST"

echo "=== [PART 3] Starting CellBender Denoising (Docker) ==="

DOCKER_DIR="$ROOT_DIR/docker"
CB_DOCKERFILE="cellbender.Dockerfile"
CB_IMAGE_NAME="cellbender:pipeline"

function prepare_docker_image() {
    local IMG_NAME=$1
    local DOCKERFILE=$2

    echo "=== [Docker Check] Image: $IMG_NAME ==="
    if [[ "$(docker images -q $IMG_NAME 2> /dev/null)" == "" ]]; then
        if [ ! -f "$DOCKER_DIR/$DOCKERFILE" ]; then
            echo "Error: Dockerfile ($DOCKERFILE) missing!"
            exit 1
        fi

        docker build -t "$IMG_NAME" -f "$DOCKER_DIR/$DOCKERFILE" "$DOCKER_DIR"
        if [ $? -ne 0 ]; then
            echo "Error: Docker build failed!"
            exit 1
        fi
        echo "  -> Build Complete!"
    else
        echo "  -> Image FOUND. Skipping build"
    fi
}
prepare_docker_image "$CB_IMAGE_NAME" "$CB_DOCKERFILE"

count=1
while IFS= read -r SRR_ID; do
    echo "[Process: CellBender] $count / $TOTAL_COUNT : $SRR_ID"

    CELLRANGER_OUTS="$CELLRANGER_PROJECT_DIR/$SRR_ID/outs"
    INPUT_H5="$CELLRANGER_OUTS/raw_feature_bc_matrix.h5"
    METRICS_CSV="$CELLRANGER_OUTS/metrics_summary.csv"

    SAMPLE_CB_DIR="$CELLBENDER_PROJECT_DIR/$SRR_ID"
    CB_OUTPUT_FILE="$SAMPLE_CB_DIR/$SRR_ID.h5"
    mkdir -p "$SAMPLE_CB_DIR"

    if [ ! -f "$INPUT_H5" ]; then
        echo "  Warning: Input not found. Skipping."
        ((count++)); continue
    fi

    if [ -f "$CB_OUTPUT_FILE" ]; then
        echo "  -> Output exists. Skipping."
        ((count++)); continue
    fi

    if [ ! -f "$METRICS_CSV" ]; then
        echo "  Error: Critical file missing -> $METRICS_CSV"
        echo "  Cannot determine expected cells. Exiting."
        exit 1
    fi

    REAL_CELLS=$("$PIXI_EXEC" run python3 -c "
import csv
import sys

try:
    with open('$METRICS_CSV', 'r') as f:
        reader = csv.reader(f)
        next(reader)
        raw_val = next(reader)[0].replace(',', '')
        print(raw_val)
except Exception:
    print('ERROR')
")

    if [[ "$REAL_CELLS" =~ ^[0-9]+$ ]]; then
        AUTO_EXPECTED_CELLS=$REAL_CELLS
        echo "  -> Auto-detected Expected Cells: $AUTO_EXPECTED_CELLS"
    else
        echo "  Error: Failed to parse cell count from CSV (Got: '$REAL_CELLS')"
        echo "  Exiting to prevent bad analysis."
        exit 1
    fi

    AUTO_TOTAL_DROPLETS=$(( AUTO_EXPECTED_CELLS * 5 ))
    if [ "$AUTO_TOTAL_DROPLETS" -lt 15000 ]; then AUTO_TOTAL_DROPLETS=15000; fi
    
    echo "  -> Setting Total Droplets: $AUTO_TOTAL_DROPLETS"
    docker run --rm --gpus all \
        -v "$CELLRANGER_OUTS":/input:ro \
        -v "$SAMPLE_CB_DIR":/output \
        -u $(id -u):$(id -g) \
        -e MPLCONFIGDIR=/tmp \
        -e HOME=/tmp \
        "$CB_IMAGE_NAME" \
        remove-background \
        --input /input/raw_feature_bc_matrix.h5 \
        --output /output/$SRR_ID.h5 \
        --cuda \
        --expected-cells $AUTO_EXPECTED_CELLS \
        --total-droplets-included $AUTO_TOTAL_DROPLETS \
        --fpr 0.01 \
        --epochs 150

    if [ $? -eq 0 ]; then
        echo "  -> Success."
    else
        echo "  -> Failed."
        exit 1
    fi

    ((count++))
done < "$SRR_LIST"

echo "=== [PART 4] Starting COMPOSITE Doublet Detection ==="

PYTHON_SCRIPT="$ROOT_DIR/pipeline/utils/doublet_composite.py"
SCCOMPOSITE_DIR="$ROOT_DIR/$SCCOMPOSITE"
SCCOMPOSITE_PROJECT_DIR="$SCCOMPOSITE_DIR/$PROJECT_ID"
SUMMARY_CSV="$SCCOMPOSITE_PROJECT_DIR/summary_gof.csv"

mkdir -p "$SCCOMPOSITE_PROJECT_DIR"

echo "SampleID,GOF_Score,Status,Note" > "$SUMMARY_CSV"

count=1
while IFS= read -r SRR_ID; do
    echo "[Process: COMPOSITE] $count / $TOTAL_COUNT : $SRR_ID"
    
    CB_OUTPUT_DIR="$CELLBENDER_PROJECT_DIR/$SRR_ID"
    INPUT_TARGET="$CB_OUTPUT_DIR/${SRR_ID}_filtered.h5"
    
    SCCOMPOSITE_OUT="$SCCOMPOSITE_PROJECT_DIR/$SRR_ID"
    mkdir -p "$SCCOMPOSITE_OUT"
    
    if [ ! -f "$INPUT_TARGET" ]; then
        echo "  Warning: Input ($INPUT_TARGET) not found. Skipping."
        echo "$SRR_ID,N/A,MISSING_INPUT,Input H5 not found" >> "$SUMMARY_CSV"
        ((count++)); continue
    fi

    "$PIXI_EXEC" run python3 "$PYTHON_SCRIPT" \
        --input "$INPUT_TARGET" \
        --outdir "$SCCOMPOSITE_OUT" \
        --sample_id "$SRR_ID" > "$SCCOMPOSITE_OUT/run.log" 2>&1
    
    PY_EXIT_CODE=$?
    echo "RUN FINISHED"

    if [ $PY_EXIT_CODE -eq 0 ]; then
        
        GOF_FILE="$SCCOMPOSITE_OUT/${SRR_ID}_gof_score.txt"
        
        if [ -f "$GOF_FILE" ]; then
            GOF_VAL=$(cat "$GOF_FILE")
        else
            GOF_VAL="0"
        fi
        IS_GOOD_FIT=$(echo "$GOF_VAL >= 3.0" | bc -l)
        
        if [ "$IS_GOOD_FIT" -eq 1 ]; then
            STATUS="PASS"
            NOTE="Reliable"
        else
            STATUS="WARNING"
            NOTE="Low GOF (<3.0)"
        fi
        
        echo "  -> Success. GOF: $GOF_VAL ($STATUS)"
        echo "$SRR_ID,$GOF_VAL,$STATUS,$NOTE" >> "$SUMMARY_CSV"
        
    else
        echo "  -> Error. Check logs at $SCCOMPOSITE_OUT/run.log"
        echo "$SRR_ID,N/A,FAILED,Script Error" >> "$SUMMARY_CSV"
    fi
    
    ((count++))
done < "$SRR_LIST"

echo "=== [PART 5] Merging all clean .h5ad files ==="

MERGE_SCRIPT="$ROOT_DIR/pipeline/utils/merge_h5ad.py"
H5AD_MATRIX_DIR="$ROOT_DIR/$H5AD_MATRIX"

"$PIXI_EXEC" run python3 "$MERGE_SCRIPT" \
    --input_dir "$SCCOMPOSITE_PROJECT_DIR" \
    --project_id "$PROJECT_ID" \
    --output_dir "$H5AD_MATRIX_DIR"

MERGE_EXIT_CODE=$?
if [ $MERGE_EXIT_CODE -eq 0 ]; then
    echo "  -> Successfully created: $H5AD_MATRIX_DIR/$PROJECT_ID.h5ad"
else
    echo "  -> Error: Merging failed."
    exit 1
fi

echo "Pipeline completed!!!"
