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
PROJECT_ID="PRJNA1103232"  
REFERENCE_GENOME="GRCm39"  # "GRCh38"

# Input mode:
#   sra   : use SRA accession, download .sra, convert with fasterq-dump
#   local : use local 10x FASTQ zip files
#   immport : use local ImmPort downloaded FASTQ files without any folder
INPUT_MODE="sra"
IMMPORT_KEY="primary_gut"
# ==========================================

TOOL_PATH="$ROOT_DIR/.tools"
CELLRANGER_PATH=$(ls -d "$TOOL_PATH"/cellranger-* | head -n 1)
if [[ -z "${CELLRANGER_PATH:-}" || ! -d "$CELLRANGER_PATH" ]]; then
    echo "Error: Cell Ranger directory not found under $TOOL_PATH"
    exit 1
fi
export PATH="$CELLRANGER_PATH:$PATH"

PROJECT_ROOT_DIR="$DATA_DIR/$PROJECT_ID"

CELLRANGER_DIR="$PROJECT_ROOT_DIR/cellranger"
CELLBENDER_DIR="$PROJECT_ROOT_DIR/cellbender"
SCDBLFINDER_DIR="$PROJECT_ROOT_DIR/scdblfinder"
META_DIR="$PROJECT_ROOT_DIR/meta_data"

mkdir -p "$CELLRANGER_DIR"
mkdir -p "$CELLBENDER_DIR"
mkdir -p "$SCDBLFINDER_DIR"
mkdir -p "$META_DIR"

LIBRARY_LIST="$META_DIR/library_list.txt"
: > "$LIBRARY_LIST"

# ==========================================
# [PART 1] Download
# ==========================================

if [[ "$INPUT_MODE" == "sra" ]]; then
    echo "=== [PART 1] SRA mode: Detecting Project Type and Download ==="
    
    SRA_DIR="$PROJECT_ROOT_DIR/sra"
    mkdir -p "$SRA_DIR"

    RUNINFO_CSV="$META_DIR/runinfo.csv"
    SRA_XML="$META_DIR/sra_metadata.xml"

    LIBRARY_MANIFEST="$META_DIR/library_manifest.txt"
    RUN_MANIFEST="$META_DIR/run_manifest.txt"
    DOWNLOAD_INPUT_LIST="$META_DIR/download_input_list.txt"
    
    echo "=== Fetching SRA IDs ==="
    SRA_IDS=$(curl -sS "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=$PROJECT_ID&retmax=5000" | grep -o '<Id>[^<]*' | sed 's/<Id>//' | paste -sd "," -)
    echo "=== Downloading RunInfo ==="
    curl -sS -X POST -d "db=sra&id=$SRA_IDS" "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?rettype=runinfo&retmode=text" > "$RUNINFO_CSV"
    echo "=== Downloading SRA XML ==="
    curl -sS -X POST -d "db=sra&id=$SRA_IDS" "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?rettype=docsum&retmode=xml" > "$SRA_XML"

    SRA_RUNINFO_PARSE_PYTHON_SCRIPT="$ROOT_DIR/pipeline/utils/single_cell/sra_runinfo_parse.py"
    "$PIXI_EXEC" run python3 "$SRA_RUNINFO_PARSE_PYTHON_SCRIPT" \
        --runinfo_csv "$RUNINFO_CSV" \
        --sra_dir "$SRA_DIR" \
        --library_manifest "$LIBRARY_MANIFEST" \
        --run_manifest "$RUN_MANIFEST" \
        --download_input_list "$DOWNLOAD_INPUT_LIST"

    if [[ -f "$LIBRARY_MANIFEST" ]]; then
        tail -n +2 "$LIBRARY_MANIFEST" | awk -F'\t' '{print $1}' >> "$LIBRARY_LIST"
    fi

    TOTAL_COUNT=$(wc -l < "$LIBRARY_LIST")
    echo "Found $TOTAL_COUNT experiment-level samples in $PROJECT_ID"
    if [[ "$TOTAL_COUNT" -eq 0 ]]; then
        exit 1
    fi

    echo "=== [PART 1] Starting High-Speed Download (aria2c) ==="

    if [ -s "$DOWNLOAD_INPUT_LIST" ]; then
        aria2c -i "$DOWNLOAD_INPUT_LIST" \
            -x 16 -s 16 -j 8 \
            --file-allocation=none \
            --disk-cache=2048M \
            --summary-interval=0
    else
        echo "All files seem to be downloaded."
    fi

elif [[ "$INPUT_MODE" == "local" ]]; then
    echo "=== [PART 1] Local mode: Detecting local zip files ==="

    find "$PROJECT_ROOT_DIR" -maxdepth 1 -type f -name "*.zip" | sort | while IFS= read -r ZIP_FILE; do
        LIBRARY_ID=$(basename "$ZIP_FILE")
        LIBRARY_ID="${LIBRARY_ID%.zip}"

        echo -e "${LIBRARY_ID}" >> "$LIBRARY_LIST"
    done

    TOTAL_COUNT=$(wc -l < "$LIBRARY_LIST")
    echo "Found $TOTAL_COUNT samples in $PROJECT_ID"
    if [[ "$TOTAL_COUNT" -eq 0 ]]; then
        exit 1
    fi

elif [[ "$INPUT_MODE" == "immport" ]]; then
    echo "=== [PART 1] ImmPort mode: Selecting local ImmPort downloaded FASTQ files ==="
    PROJECT_FASTQ_DIR="$PROJECT_ROOT_DIR/fastq"
    if [[ ! -d "$PROJECT_FASTQ_DIR" ]]; then
        PROJECT_FASTQ_DIR="$PROJECT_ROOT_DIR"
    fi
    FASTQ_METADATA_FILE="$META_DIR/meta_data.txt"
    echo "ImmPort key: ${IMMPORT_KEY:-<all>}"

    IMMPORT_FASTQ_MANIFEST_PYTHON_SCRIPT="$ROOT_DIR/pipeline/utils/single_cell/immport_fastq_manifest.py"
    "$PIXI_EXEC" run python3 "$IMMPORT_FASTQ_MANIFEST_PYTHON_SCRIPT" \
        --metadata_file "$FASTQ_METADATA_FILE" \
        --key "$IMMPORT_KEY" \
        --library_list "$LIBRARY_LIST"

    TOTAL_COUNT=$(wc -l < "$LIBRARY_LIST")
    echo "Found $TOTAL_COUNT ImmPort-mode libraries in $PROJECT_ID"
    if [[ "$TOTAL_COUNT" -eq 0 ]]; then
        exit 1
    fi

else
    echo "Error: unknown INPUT_MODE=$INPUT_MODE"
    echo "Allowed values: sra, local, immport"
    exit 1
fi

# ==========================================
# [PART 2] Cellranger
# ==========================================

run_cellranger_count() {
    local LIBRARY_ID="$1"
    local FASTQ_LIBRARY_DIR="$2"
    local CHEMISTRY="${3:-auto}"
    local CHEMISTRY_ARGS=()

    if [[ "$CHEMISTRY" != "auto" ]]; then
        CHEMISTRY_ARGS=(--chemistry="$CHEMISTRY")
    fi

    if [[ -d "$CELLRANGER_DIR/$LIBRARY_ID" && ! -d "$CELLRANGER_DIR/$LIBRARY_ID/outs" ]]; then
        rm -rf "$CELLRANGER_DIR/$LIBRARY_ID"
    fi

    (
        cd "$CELLRANGER_DIR" || exit 1
        echo "  -> Running Cell Ranger (chemistry=$CHEMISTRY)..."

        cellranger count \
            --id="$LIBRARY_ID" \
            --transcriptome="$ROOT_DIR/references/sc/$REFERENCE_GENOME" \
            --fastqs="$FASTQ_LIBRARY_DIR" \
            --sample="$LIBRARY_ID" \
            --localcores="$N_THREADS" \
            --create-bam=false \
            "${CHEMISTRY_ARGS[@]}" \
            > "${LIBRARY_ID}.log"
    )
}

if [[ "$INPUT_MODE" == "sra" ]]; then
    echo "=== [PART 2] SRA mode: per-library FASTQ conversion + Cell Ranger ==="

    convert_sra_one_run() {
        local RUN_ID="$1"
        local LIBRARY_ID="$2"
        local LANE_ID="$3"
        local SRA_FILE=$(echo "$4" | tr -d '\r' | xargs)
        local FASTQ_LIBRARY_DIR="$5"

        local TMP_RUN_DIR="$FASTQ_DATA/_tmp_${RUN_ID}"

        rm -rf "$TMP_RUN_DIR"
        mkdir -p "$TMP_RUN_DIR" "$FASTQ_LIBRARY_DIR"

        if [[ ! -f "$SRA_FILE" ]]; then
            echo "  Error: SRA file not found: $SRA_FILE"
            exit 1
        fi

        echo "  -> Converting $RUN_ID to FastQ..."

        fasterq-dump "$SRA_FILE" \
            --split-files \
            --include-technical \
            --threads "$N_THREADS" \
            --mem 15G \
            --outdir "$TMP_RUN_DIR" \
            --temp "$DUMP_DIR"

        declare -a indices=()
        declare -a reads=()

        while IFS= read -r f; do
            LEN=$(awk 'NR==2 {print length($0); exit}' "$f")
            LEN=${LEN:-0}

            if (( LEN < 20 )); then
                indices+=("$f")
            else
                reads+=("$f")
            fi
        done < <(find "$TMP_RUN_DIR" -maxdepth 1 -type f -name "*.fastq" | sort -V)

        if [ ${#reads[@]} -eq 2 ]; then
            i_cnt=1
            for f in "${indices[@]}"; do
                mv "$f" "$FASTQ_LIBRARY_DIR/${LIBRARY_ID}_S1_${LANE_ID}_I${i_cnt}_001.fastq"
                i_cnt=$((i_cnt + 1))
            done

            rA="${reads[0]}"
            rB="${reads[1]}"
            
            lenA=$(awk 'NR==2 {print length($0); exit}' "$rA")
            lenB=$(awk 'NR==2 {print length($0); exit}' "$rB")

            if [ "$lenA" -lt "$lenB" ]; then
                mv "$rA" "$FASTQ_LIBRARY_DIR/${LIBRARY_ID}_S1_${LANE_ID}_R1_001.fastq"
                mv "$rB" "$FASTQ_LIBRARY_DIR/${LIBRARY_ID}_S1_${LANE_ID}_R2_001.fastq"
            elif [ "$lenA" -gt "$lenB" ]; then
                mv "$rB" "$FASTQ_LIBRARY_DIR/${LIBRARY_ID}_S1_${LANE_ID}_R1_001.fastq"
                mv "$rA" "$FASTQ_LIBRARY_DIR/${LIBRARY_ID}_S1_${LANE_ID}_R2_001.fastq"
            else
                mv "$rA" "$FASTQ_LIBRARY_DIR/${LIBRARY_ID}_S1_${LANE_ID}_R1_001.fastq"
                mv "$rB" "$FASTQ_LIBRARY_DIR/${LIBRARY_ID}_S1_${LANE_ID}_R2_001.fastq"
            fi
        elif [ ${#reads[@]} -eq 1 ] && [ ${#indices[@]} -eq 3 ]; then
            local BARCODE_FASTQ=""
            local UMI_FASTQ=""
            local INDEX_FASTQ=""
            local INDEX_COUNT=0

            for f in "${indices[@]}"; do
                LEN=$(awk 'NR==2 {print length($0); exit}' "$f")

                if [ "$LEN" -eq 14 ]; then
                    BARCODE_FASTQ="$f"
                elif [ "$LEN" -eq 10 ]; then
                    UMI_FASTQ="$f"
                elif [ "$LEN" -eq 8 ]; then
                    INDEX_FASTQ="$f"
                    INDEX_COUNT=$((INDEX_COUNT + 1))
                fi
            done

            if [[ -n "$BARCODE_FASTQ" && -n "$UMI_FASTQ" && -n "$INDEX_FASTQ" && "$INDEX_COUNT" -eq 1 ]]; then
                echo "  -> Detected split 10x v1 reads (14bp barcode + 10bp UMI + 8bp sample index + cDNA). Renaming as SC3Pv1 FASTQs..."

                mv "${reads[0]}" "$FASTQ_LIBRARY_DIR/${LIBRARY_ID}_S1_${LANE_ID}_R1_001.fastq"
                mv "$UMI_FASTQ" "$FASTQ_LIBRARY_DIR/${LIBRARY_ID}_S1_${LANE_ID}_R2_001.fastq"
                mv "$BARCODE_FASTQ" "$FASTQ_LIBRARY_DIR/${LIBRARY_ID}_S1_${LANE_ID}_I1_001.fastq"
                mv "$INDEX_FASTQ" "$FASTQ_LIBRARY_DIR/${LIBRARY_ID}_S1_${LANE_ID}_I2_001.fastq"
                CELLRANGER_CHEMISTRY="SC3Pv1"
            else
                echo "  Error: Expected split 10x v1 lengths 14/10/8 plus one cDNA read for $RUN_ID"
                for f in "${indices[@]}" "${reads[@]}"; do
                    LEN=$(awk 'NR==2 {print length($0); exit}' "$f")
                    echo "    $(basename "$f"): ${LEN:-0}bp"
                done
                exit 1
            fi
        else
            echo "  Error: Expected 2 read files (>=20bp), found ${#reads[@]} for $RUN_ID"
            for f in "${indices[@]}" "${reads[@]}"; do
                LEN=$(awk 'NR==2 {print length($0); exit}' "$f")
                echo "    $(basename "$f"): ${LEN:-0}bp"
            done
            exit 1
        fi

        rm -rf "$TMP_RUN_DIR"
    }

    count=1
    tail -n +2 "$LIBRARY_MANIFEST" | while IFS=$'\t' read -r LIBRARY_ID SAMPLE_ID EXPERIMENT RUNS; do
        echo "[SRA] $count / $TOTAL_COUNT : $LIBRARY_ID"

        OUTPUT_PATH="$CELLRANGER_DIR/$LIBRARY_ID"
        FASTQ_LIBRARY_DIR="$FASTQ_DATA/$LIBRARY_ID"

        if [[ -d "$OUTPUT_PATH/outs" ]]; then
            echo "  -> Cell Ranger output exists. Skipping."
            count=$((count + 1))
            continue
        fi

        rm -rf "$FASTQ_LIBRARY_DIR"
        mkdir -p "$FASTQ_LIBRARY_DIR"
        CELLRANGER_CHEMISTRY="auto"

        while IFS=$'\t' read -r RUN_ID LANE_ID SRA_FILE; do
            convert_sra_one_run \
                "$RUN_ID" \
                "$LIBRARY_ID" \
                "$LANE_ID" \
                "$SRA_FILE" \
                "$FASTQ_LIBRARY_DIR"
        done < <(
            awk -F'\t' -v lid="$LIBRARY_ID" '
                NR > 1 && $2 == lid {
                    print $1 "\t" $5 "\t" $6
                }
            ' "$RUN_MANIFEST"
        )

        if run_cellranger_count "$LIBRARY_ID" "$FASTQ_LIBRARY_DIR" "$CELLRANGER_CHEMISTRY"; then
            echo "  -> Cell Ranger finished"
            rm -rf "$FASTQ_LIBRARY_DIR"
        else
            echo "  -> Warning: Cell Ranger failed"
            rm -rf "$FASTQ_LIBRARY_DIR"
        fi

        count=$((count + 1))
    done

elif [[ "$INPUT_MODE" == "local" ]]; then
    echo "=== [PART 2] Local mode: per-library unzip + Cell Ranger ==="

    count=1
    while IFS= read -r LIBRARY_ID; do
        echo "[Local] $count / $TOTAL_COUNT : $LIBRARY_ID"

        ZIP_FILE="$PROJECT_ROOT_DIR/${LIBRARY_ID}.zip"
        FASTQ_LIBRARY_DIR="$FASTQ_DATA/$LIBRARY_ID"
        OUTPUT_PATH="$CELLRANGER_DIR/$LIBRARY_ID"

        if [[ -d "$OUTPUT_PATH/outs" ]]; then
            echo "  -> Cell Ranger output exists. Skipping."
            count=$((count + 1))
            continue
        fi

        if [[ ! -f "$ZIP_FILE" ]]; then
            echo "  Error: zip file not found: $ZIP_FILE"
            exit 1
        fi

        rm -rf "$FASTQ_LIBRARY_DIR"
        mkdir -p "$FASTQ_LIBRARY_DIR"

        echo "  -> Extracting FASTQ files from zip..."
        EXTRACT_FASTQ_PYTHON_SCRIPT="$ROOT_DIR/pipeline/utils/single_cell/extract_fastq_from_zip.py"

        "$PIXI_EXEC" run python3 "$EXTRACT_FASTQ_PYTHON_SCRIPT" \
            --zip_file "$ZIP_FILE" \
            --outdir "$FASTQ_LIBRARY_DIR"

        if ! ls "$FASTQ_LIBRARY_DIR"/*_R1_001.fastq.gz >/dev/null 2>&1; then
            echo "  Error: no R1 FASTQ found after unzip: $ZIP_FILE"
            exit 1
        fi

        if run_cellranger_count "$LIBRARY_ID" "$FASTQ_LIBRARY_DIR"; then
            echo "  -> Cell Ranger finished. Cleaning FASTQ..."
            rm -rf "$FASTQ_LIBRARY_DIR"
        else
            echo "  -> Warning: Cell Ranger failed"
            exit 1
        fi

        count=$((count + 1))
    done < "$LIBRARY_LIST"

elif [[ "$INPUT_MODE" == "immport" ]]; then
    echo "=== [PART 2] ImmPort mode: per-library symlink staging + Cell Ranger ==="

    stage_fastq_library() {
        local LIBRARY_ID="$1"
        local FASTQ_LIBRARY_DIR="$2"
        local STAGED_COUNT=0

        rm -rf "$FASTQ_LIBRARY_DIR"
        mkdir -p "$FASTQ_LIBRARY_DIR"

        while IFS=$'\t' read -r SOURCE_PATH CELLRANGER_FASTQ_NAME; do
            if [[ ! -f "$SOURCE_PATH" ]]; then
                echo "  Error: FASTQ file not found: $SOURCE_PATH"
                exit 1
            fi

            ln -sf "$SOURCE_PATH" "$FASTQ_LIBRARY_DIR/$CELLRANGER_FASTQ_NAME"
            STAGED_COUNT=$((STAGED_COUNT + 1))
        done < <(
            awk -F'\t' -v lid="$LIBRARY_ID" -v key="$IMMPORT_KEY" -v fastq_dir="$PROJECT_FASTQ_DIR" '
                NR == 1 {
                    for (i = 1; i <= NF; i++) h[$i] = i
                    next
                }
                (key == "" || $(h["ARM Name"]) == key) {
                    source = $(h["File Name"])
                    dest = $(h["Original File Name"])
                    if (dest == "") dest = source
                    sub(/\.fastq\.[0-9]+\.gz$/, ".fastq.gz", dest)
                    sub(/\.fq\.[0-9]+\.gz$/, ".fq.gz", dest)

                    lib = dest
                    sub(/_S[0-9]+_L[0-9]+_[RI][0-9]_001\.(fastq|fq)\.gz$/, "", lib)
                    if (lib == lid) print fastq_dir "/" source "\t" dest
                }
            ' "$FASTQ_METADATA_FILE"
        )

        if [[ "$STAGED_COUNT" -eq 0 ]]; then
            echo "  Error: no FASTQ records found in $FASTQ_METADATA_FILE for $LIBRARY_ID"
            exit 1
        fi
    }

    count=1
    while IFS= read -r LIBRARY_ID; do
        echo "[ImmPort] $count / $TOTAL_COUNT : $LIBRARY_ID"

        FASTQ_LIBRARY_DIR="$FASTQ_DATA/$LIBRARY_ID"
        OUTPUT_PATH="$CELLRANGER_DIR/$LIBRARY_ID"

        if [[ -d "$OUTPUT_PATH/outs" ]]; then
            echo "  -> Cell Ranger output exists. Skipping."
            count=$((count + 1))
            continue
        fi

        echo "  -> Staging Cell Ranger-compatible FASTQ symlinks..."
        stage_fastq_library "$LIBRARY_ID" "$FASTQ_LIBRARY_DIR"

        if ! ls "$FASTQ_LIBRARY_DIR"/*_R1_001.fastq.gz >/dev/null 2>&1; then
            echo "  Error: no R1 FASTQ found after staging: $FASTQ_LIBRARY_DIR"
            exit 1
        fi

        if run_cellranger_count "$LIBRARY_ID" "$FASTQ_LIBRARY_DIR"; then
            echo "  -> Cell Ranger finished. Cleaning staged FASTQ symlinks..."
            rm -rf "$FASTQ_LIBRARY_DIR"
        else
            echo "  -> Warning: Cell Ranger failed"
            exit 1
        fi

        count=$((count + 1))
    done < "$LIBRARY_LIST"
fi

# ==========================================
# [PART 3] Cellbender
# ==========================================

echo "=== [PART 3] Starting CellBender Denoising ==="

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

SUCCESS_SAMPLES_CSV="$CELLBENDER_DIR/success_samples.csv"
CB_STATUS_CSV="$CELLBENDER_DIR/cellbender_status.csv"

echo -e "LibraryID,FilteredH5" > "$SUCCESS_SAMPLES_CSV"
echo "LibraryID,Status,Reason,ExpectedCells,TotalBarcodes,MeanReadsPerCell,MedianGenesPerCell,MedianUMIPerCell,OutputH5" > "$CB_STATUS_CSV"

count=1
while IFS= read -r LIBRARY_ID; do
    echo "[Process: CellBender] $count / $TOTAL_COUNT : $LIBRARY_ID"

    CELLRANGER_OUTS="$CELLRANGER_DIR/$LIBRARY_ID/outs"
    INPUT_H5="$CELLRANGER_OUTS/raw_feature_bc_matrix.h5"
    METRICS_CSV="$CELLRANGER_OUTS/metrics_summary.csv"

    SAMPLE_CB_DIR="$CELLBENDER_DIR/$LIBRARY_ID"
    CB_FILTERED_FILE="$SAMPLE_CB_DIR/${LIBRARY_ID}_filtered.h5"
    mkdir -p "$SAMPLE_CB_DIR"
    
    if [ -f "$CB_FILTERED_FILE" ]; then
        echo "  -> Existing filtered output found. Reusing."
        echo "$LIBRARY_ID,$CB_FILTERED_FILE" >> "$SUCCESS_SAMPLES_CSV"
        echo "$LIBRARY_ID,REUSE_EXISTING,OK,N/A,N/A,N/A,N/A,N/A,$CB_FILTERED_FILE" >> "$CB_STATUS_CSV"
        count=$((count + 1)); continue
    fi

    if [ ! -f "$INPUT_H5" ]; then
        echo "  Warning: Input not found. Skipping."
        echo "$LIBRARY_ID,SKIP,RAW_INPUT_MISSING,N/A,N/A,N/A,N/A,N/A," >> "$CB_STATUS_CSV"
        count=$((count + 1)); continue
    fi

    if [ ! -f "$METRICS_CSV" ]; then
        echo "  Warning: metrics_summary.csv missing. Skipping."
        echo "$LIBRARY_ID,SKIP,METRICS_MISSING,N/A,N/A,N/A,N/A,N/A," >> "$CB_STATUS_CSV"
        count=$((count + 1)); continue
    fi

    METRICS=$("$PIXI_EXEC" run python3 -c '
import csv
import sys

path = sys.argv[1]
with open(path, "r", newline="") as f:
    reader = csv.DictReader(f)
    row = next(reader)

def parse_intlike(x):
    x = str(x).replace(",", "").strip()
    if x == "":
        return "NA"
    try:
        return str(int(float(x)))
    except Exception:
        return "NA"

print("|".join([
    parse_intlike(row.get("Estimated Number of Cells", "")),
    parse_intlike(row.get("Mean Reads per Cell", "")),
    parse_intlike(row.get("Median Genes per Cell", "")),
    parse_intlike(row.get("Median UMI Counts per Cell", "")),
]))
' "$METRICS_CSV")
    IFS='|' read -r REAL_CELLS MEAN_READS MEDIAN_GENES MEDIAN_UMI <<< "$METRICS"

    if ! [[ "$REAL_CELLS" =~ ^[0-9]+$ ]]; then
        echo "  Warning: failed to parse metrics. Skipping."
        echo "$LIBRARY_ID,SKIP,METRICS_PARSE_FAILED,$REAL_CELLS,N/A,$MEAN_READS,$MEDIAN_GENES,$MEDIAN_UMI," >> "$CB_STATUS_CSV"
        count=$((count + 1))
        continue
    fi

    TOTAL_BARCODES=$("$PIXI_EXEC" run python3 -c '
import h5py
import sys

with h5py.File(sys.argv[1], "r") as f:
    print(len(f["matrix"]["barcodes"]))
' "$INPUT_H5")

    AUTO_EXPECTED_CELLS="$REAL_CELLS"
    AUTO_TOTAL_DROPLETS=$(( AUTO_EXPECTED_CELLS * 3 ))
    if [ "$AUTO_TOTAL_DROPLETS" -lt 15000 ]; then AUTO_TOTAL_DROPLETS=15000; fi
    if [ "$AUTO_TOTAL_DROPLETS" -gt "$TOTAL_BARCODES" ]; then AUTO_TOTAL_DROPLETS="$TOTAL_BARCODES"; fi
    
    echo "  -> Metrics: expected_cells=$REAL_CELLS, mean_reads=$MEAN_READS, median_genes=$MEDIAN_GENES, median_umi=$MEDIAN_UMI, total_barcodes=$TOTAL_BARCODES"
    echo "  -> Running CellBender with expected_cells=$AUTO_EXPECTED_CELLS total_droplets=$AUTO_TOTAL_DROPLETS"

    if docker run --rm --gpus all \
        -v "$CELLRANGER_OUTS":/input:ro \
        -v "$SAMPLE_CB_DIR":/output \
        -u "$(id -u):$(id -g)" \
        -e MPLCONFIGDIR=/tmp \
        -e HOME=/tmp \
        "$CB_IMAGE_NAME" \
        remove-background \
        --input /input/raw_feature_bc_matrix.h5 \
        --output /output/$LIBRARY_ID.h5 \
        --cuda \
        --expected-cells "$AUTO_EXPECTED_CELLS" \
        --total-droplets-included "$AUTO_TOTAL_DROPLETS" \
        --fpr 0.01 \
        --epochs 150
    then
        if [ -f "$CB_FILTERED_FILE" ]; then
            echo "  -> CellBender success."
            echo "$LIBRARY_ID,$CB_FILTERED_FILE" >> "$SUCCESS_SAMPLES_CSV"
            echo "$LIBRARY_ID,SUCCESS,OK,$AUTO_EXPECTED_CELLS,$TOTAL_BARCODES,$MEAN_READS,$MEDIAN_GENES,$MEDIAN_UMI,$CB_FILTERED_FILE" >> "$CB_STATUS_CSV"
        else
            echo "  -> CellBender finished but filtered output missing."
            echo "$LIBRARY_ID,EXCLUDE,FILTERED_OUTPUT_MISSING,$AUTO_EXPECTED_CELLS,$TOTAL_BARCODES,$MEAN_READS,$MEDIAN_GENES,$MEDIAN_UMI," >> "$CB_STATUS_CSV"
        fi
    else
        echo "  -> CellBender failed. Sample excluded."
        echo "$LIBRARY_ID,EXCLUDE,CELLBENDER_FAILED,$AUTO_EXPECTED_CELLS,$TOTAL_BARCODES,$MEAN_READS,$MEDIAN_GENES,$MEDIAN_UMI," >> "$CB_STATUS_CSV"
    fi

    count=$((count + 1))
done < "$LIBRARY_LIST"

# ==========================================
# [PART 4] scDblFinder
# ==========================================

echo "=== [PART 4] Starting scDblFinder Doublet Detection ==="

SUMMARY_CSV="$SCDBLFINDER_DIR/summary_scdblfinder.csv"

SCDBLFINDER_R_SCRIPT="$ROOT_DIR/R/scdblfinder.R"
SCDBLFINDER_PYTHON_SCRIPT="$ROOT_DIR/pipeline/utils/single_cell/scdblfinder.py"

"$PIXI_EXEC" run Rscript "$SCDBLFINDER_R_SCRIPT" \
    --manifest "$SUCCESS_SAMPLES_CSV" \
    --outdir "$SCDBLFINDER_DIR" \
    --threads 2 \
    --dbr_per1k 0.008 \
    --seed 1

while IFS=, read -r LIBRARY_ID INPUT_TARGET || [[ -n "${LIBRARY_ID:-}" ]]; do
    [[ -z "${LIBRARY_ID:-}" || "$LIBRARY_ID" == "LibraryID" ]] && continue

    OUTDIR="$SCDBLFINDER_DIR/$LIBRARY_ID"
    SINGLET_CSV="$OUTDIR/${LIBRARY_ID}_singlet_barcodes.csv"

    if [[ ! -f "$INPUT_TARGET" ]]; then
        echo "  -> Skip $LIBRARY_ID: input missing"
        continue
    fi

    if [[ ! -f "$SINGLET_CSV" ]]; then
        echo "  -> Skip $LIBRARY_ID: singlet barcode file missing"
        continue
    fi

    "$PIXI_EXEC" run python3 "$SCDBLFINDER_PYTHON_SCRIPT" \
        --input "$INPUT_TARGET" \
        --singlets "$SINGLET_CSV" \
        --outdir "$OUTDIR" \
        --sample_id "$LIBRARY_ID"
done < "$SUCCESS_SAMPLES_CSV"

echo "LibraryID,TotalCells,RemovedZeroCountCells,PredictedDoublets,PredictedSinglets,ObservedDoubletFraction,ExpectedDoubletFraction,dbr_per1k,scDblFinderVersion,Status,Note" > "$SUMMARY_CSV"

for f in "$SCDBLFINDER_DIR"/*/*_summary.csv; do
    [[ -f "$f" ]] && tail -n +2 "$f" >> "$SUMMARY_CSV"
done

# ==========================================
# [PART 5] Merging each h5ad files
# ==========================================

echo "=== [PART 5] Merging all clean .h5ad files ==="

MERGE_SCRIPT="$ROOT_DIR/pipeline/utils/single_cell/merge_h5ad.py"
H5AD_MATRIX_DIR="$ROOT_DIR/$MERGED_H5AD"

MERGE_CMD=(
    "$PIXI_EXEC" run python3 "$MERGE_SCRIPT"
    --input_dir "$SCDBLFINDER_DIR"
    --project_id "$PROJECT_ID"
    --output_dir "$H5AD_MATRIX_DIR"
)

if [[ "$INPUT_MODE" == "sra" ]]; then
    MERGE_CMD+=(--library_manifest "$LIBRARY_MANIFEST")
    MERGE_CMD+=(--sra_xml "$SRA_XML")
fi

"${MERGE_CMD[@]}"

MERGE_EXIT_CODE=$?
if [ $MERGE_EXIT_CODE -eq 0 ]; then
    echo "  -> Successfully created: $H5AD_MATRIX_DIR/$PROJECT_ID.h5ad"
else
    echo "  -> Error: Merging failed."
    exit 1
fi

echo "Pipeline completed!!!"
