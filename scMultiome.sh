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
PROJECT_ID="PRJNA1054583"
REFERENCE_GENOME="GRCm39"

# Input mode:
#   sra   : use paired SRA scRNAseq/scATACseq runs for 10x Multiome
INPUT_MODE="sra"
# ==========================================

TOOL_PATH="$ROOT_DIR/.tools"
CELLRANGER_ARC_PATH=$(ls -d "$TOOL_PATH"/cellranger-arc-* | head -n 1)
if [[ -z "${CELLRANGER_ARC_PATH:-}" || ! -d "$CELLRANGER_ARC_PATH" ]]; then
    echo "Error: Cell Ranger ARC directory not found under $TOOL_PATH"
    exit 1
fi
export PATH="$CELLRANGER_ARC_PATH:$PATH"

PROJECT_ROOT_DIR="$DATA_DIR/$PROJECT_ID"

CELLRANGER_ARC_DIR="$PROJECT_ROOT_DIR/cellranger-arc"
CELLBENDER_DIR="$PROJECT_ROOT_DIR/cellbender_gex"
SCDBLFINDER_DIR="$PROJECT_ROOT_DIR/scdblfinder_gex"
SNAPATAC_DIR="$PROJECT_ROOT_DIR/snapatac2"
META_DIR="$PROJECT_ROOT_DIR/meta_data"

ARC_REFERENCE_DIR="$ROOT_DIR/references/sc_multiomics/$REFERENCE_GENOME"
CHROM_SIZES="$ARC_REFERENCE_DIR/fasta/genome.fa.fai"
GENE_ANNOTATION="$ARC_REFERENCE_DIR/genes/genes.gtf.gz"

mkdir -p "$CELLRANGER_ARC_DIR"
mkdir -p "$CELLBENDER_DIR"
mkdir -p "$SCDBLFINDER_DIR"
mkdir -p "$META_DIR"

[[ -f "$CHROM_SIZES" && -f "$GENE_ANNOTATION" ]] || {
    echo "Error: SnapATAC2 reference files are missing under $ARC_REFERENCE_DIR"
    exit 1
}

LIBRARY_LIST="$META_DIR/library_list.txt"
: > "$LIBRARY_LIST"

# ==========================================
# [PART 1] Download
# ==========================================

if [[ "$INPUT_MODE" == "sra" ]]; then
    echo "=== [PART 1] SRA mode: Detecting Paired Multiome Libraries and Download ==="

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

    SRA_MULTIOME_RUNINFO_PARSE_PYTHON_SCRIPT="$ROOT_DIR/pipeline/utils/single_cell/sra_multiome_runinfo_parse.py"
    "$PIXI_EXEC" run python3 "$SRA_MULTIOME_RUNINFO_PARSE_PYTHON_SCRIPT" \
        --runinfo_csv "$RUNINFO_CSV" \
        --sra_xml "$SRA_XML" \
        --sra_dir "$SRA_DIR" \
        --library_manifest "$LIBRARY_MANIFEST" \
        --run_manifest "$RUN_MANIFEST" \
        --download_input_list "$DOWNLOAD_INPUT_LIST"

    if [[ -f "$LIBRARY_MANIFEST" ]]; then
        tail -n +2 "$LIBRARY_MANIFEST" | awk -F'\t' '{print $1}' >> "$LIBRARY_LIST"
    fi

    TOTAL_COUNT=$(wc -l < "$LIBRARY_LIST")
    echo "Found $TOTAL_COUNT paired multiome samples in $PROJECT_ID"
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
else
    echo "Error: unknown INPUT_MODE=$INPUT_MODE"
    echo "Allowed values: sra"
    exit 1
fi

# ==========================================
# [PART 2] Cellranger ARC
# ==========================================

arc_outs_complete() {
    local OUTS_DIR="$1"
    [[ -f "$OUTS_DIR/raw_feature_bc_matrix.h5" \
        && -f "$OUTS_DIR/filtered_feature_bc_matrix.h5" \
        && -f "$OUTS_DIR/atac_fragments.tsv.gz" \
        && -f "$OUTS_DIR/summary.csv" ]]
}

run_cellranger_arc_count() {
    local LIBRARY_ID="$1"
    local LIBRARIES_CSV="$2"

    if [[ -d "$CELLRANGER_ARC_DIR/$LIBRARY_ID" ]] && ! arc_outs_complete "$CELLRANGER_ARC_DIR/$LIBRARY_ID/outs"; then
        echo "Error: incomplete Cell Ranger ARC output exists; refusing to delete it:"
        echo "  $CELLRANGER_ARC_DIR/$LIBRARY_ID"
        return 1
    fi

    (
        cd "$CELLRANGER_ARC_DIR" || exit 1
        echo "  -> Running Cell Ranger ARC..."

        cellranger-arc count \
            --id="$LIBRARY_ID" \
            --reference="$ROOT_DIR/references/sc_multiomics/$REFERENCE_GENOME" \
            --libraries="$LIBRARIES_CSV" \
            --localcores="$N_THREADS" \
            --create-bam=false \
            > "${LIBRARY_ID}.log" 2>&1
    )
}

if [[ "$INPUT_MODE" == "sra" ]]; then
    echo "=== [PART 2] SRA mode: per-library FASTQ conversion + Cell Ranger ARC ==="

    convert_sra_one_run() {
        local RUN_ID="$1"
        local LIBRARY_ID="$2"
        local ASSAY="$3"
        local LIBRARY_SAMPLE="$4"
        local LANE_ID="$5"
        local SRA_FILE=$(echo "$6" | tr -d '\r' | xargs)
        local FASTQ_ASSAY_DIR="$7"

        local TMP_RUN_DIR
        TMP_RUN_DIR=$(mktemp -d "$FASTQ_DATA/.${RUN_ID}.XXXXXX")
        mkdir -p "$FASTQ_ASSAY_DIR"

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

        if [[ "$ASSAY" == "gex" ]]; then
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
                    mv "$f" "$FASTQ_ASSAY_DIR/${LIBRARY_SAMPLE}_S1_${LANE_ID}_I${i_cnt}_001.fastq"
                    i_cnt=$((i_cnt + 1))
                done

                rA="${reads[0]}"
                rB="${reads[1]}"

                lenA=$(awk 'NR==2 {print length($0); exit}' "$rA")
                lenB=$(awk 'NR==2 {print length($0); exit}' "$rB")

                if [ "$lenA" -lt "$lenB" ]; then
                    mv "$rA" "$FASTQ_ASSAY_DIR/${LIBRARY_SAMPLE}_S1_${LANE_ID}_R1_001.fastq"
                    mv "$rB" "$FASTQ_ASSAY_DIR/${LIBRARY_SAMPLE}_S1_${LANE_ID}_R2_001.fastq"
                elif [ "$lenA" -gt "$lenB" ]; then
                    mv "$rB" "$FASTQ_ASSAY_DIR/${LIBRARY_SAMPLE}_S1_${LANE_ID}_R1_001.fastq"
                    mv "$rA" "$FASTQ_ASSAY_DIR/${LIBRARY_SAMPLE}_S1_${LANE_ID}_R2_001.fastq"
                else
                    mv "$rA" "$FASTQ_ASSAY_DIR/${LIBRARY_SAMPLE}_S1_${LANE_ID}_R1_001.fastq"
                    mv "$rB" "$FASTQ_ASSAY_DIR/${LIBRARY_SAMPLE}_S1_${LANE_ID}_R2_001.fastq"
                fi
            else
                echo "  Error: Expected 2 GEX read files (>=20bp), found ${#reads[@]} for $RUN_ID"
                for f in "${indices[@]}" "${reads[@]}"; do
                    LEN=$(awk 'NR==2 {print length($0); exit}' "$f")
                    echo "    $(basename "$f"): ${LEN:-0}bp"
                done
                exit 1
            fi
        elif [[ "$ASSAY" == "atac" ]]; then
            local I1_FASTQ=""
            local I2_FASTQ=""
            declare -a reads=()

            while IFS= read -r f; do
                LEN=$(awk 'NR==2 {print length($0); exit}' "$f")
                LEN=${LEN:-0}

                if (( LEN < 15 )); then
                    I1_FASTQ="$f"
                elif (( LEN < 35 )); then
                    I2_FASTQ="$f"
                else
                    reads+=("$f")
                fi
            done < <(find "$TMP_RUN_DIR" -maxdepth 1 -type f -name "*.fastq" | sort -V)

            if [ ${#reads[@]} -eq 2 ] && [[ -n "$I2_FASTQ" ]]; then
                rA="${reads[0]}"
                rB="${reads[1]}"

                lenA=$(awk 'NR==2 {print length($0); exit}' "$rA")
                lenB=$(awk 'NR==2 {print length($0); exit}' "$rB")

                if [ "$lenA" -lt "$lenB" ]; then
                    mv "$rB" "$FASTQ_ASSAY_DIR/${LIBRARY_SAMPLE}_S1_${LANE_ID}_R1_001.fastq"
                    mv "$rA" "$FASTQ_ASSAY_DIR/${LIBRARY_SAMPLE}_S1_${LANE_ID}_R2_001.fastq"
                else
                    mv "$rA" "$FASTQ_ASSAY_DIR/${LIBRARY_SAMPLE}_S1_${LANE_ID}_R1_001.fastq"
                    mv "$rB" "$FASTQ_ASSAY_DIR/${LIBRARY_SAMPLE}_S1_${LANE_ID}_R2_001.fastq"
                fi

                if [[ -n "$I1_FASTQ" ]]; then
                    mv "$I1_FASTQ" "$FASTQ_ASSAY_DIR/${LIBRARY_SAMPLE}_S1_${LANE_ID}_I1_001.fastq"
                fi
                mv "$I2_FASTQ" "$FASTQ_ASSAY_DIR/${LIBRARY_SAMPLE}_S1_${LANE_ID}_I2_001.fastq"
            else
                echo "  Error: Expected 2 ATAC read files (>=35bp) and I2 barcode for $RUN_ID"
                for f in "$TMP_RUN_DIR"/*.fastq; do
                    [[ -f "$f" ]] || continue
                    LEN=$(awk 'NR==2 {print length($0); exit}' "$f")
                    echo "    $(basename "$f"): ${LEN:-0}bp"
                done
                exit 1
            fi
        else
            echo "  Error: unknown assay=$ASSAY"
            exit 1
        fi

        rmdir "$TMP_RUN_DIR"
    }

    write_libraries_csv() {
        local LIBRARY_ID="$1"
        local FASTQ_LIBRARY_DIR="$2"
        local LIBRARIES_CSV="$3"

        {
            echo "fastqs,sample,library_type"
            echo "$FASTQ_LIBRARY_DIR/gex,${LIBRARY_ID}_GEX,Gene Expression"
            echo "$FASTQ_LIBRARY_DIR/atac,${LIBRARY_ID}_ATAC,Chromatin Accessibility"
        } > "$LIBRARIES_CSV"
    }

    SUCCESS_SAMPLES_CSV="$CELLRANGER_ARC_DIR/success_samples.csv"
    ARC_STATUS_CSV="$CELLRANGER_ARC_DIR/cellranger_arc_status.csv"

    echo "LibraryID,ArcOuts" > "$SUCCESS_SAMPLES_CSV"
    echo "LibraryID,Status,Reason,ArcOuts" > "$ARC_STATUS_CSV"

    count=1
    tail -n +2 "$LIBRARY_MANIFEST" | while IFS=$'\t' read -r LIBRARY_ID GEX_EXPERIMENT GEX_RUNS ATAC_EXPERIMENT ATAC_RUNS; do
        echo "[SRA] $count / $TOTAL_COUNT : $LIBRARY_ID"

        OUTPUT_PATH="$CELLRANGER_ARC_DIR/$LIBRARY_ID"
        FASTQ_LIBRARY_DIR="$FASTQ_DATA/$LIBRARY_ID"
        GEX_FASTQ_DIR="$FASTQ_LIBRARY_DIR/gex"
        ATAC_FASTQ_DIR="$FASTQ_LIBRARY_DIR/atac"
        LIBRARIES_CSV="$META_DIR/${LIBRARY_ID}_libraries.csv"

        if arc_outs_complete "$OUTPUT_PATH/outs"; then
            echo "  -> Complete Cell Ranger ARC output exists. Skipping."
            echo "$LIBRARY_ID,$OUTPUT_PATH/outs" >> "$SUCCESS_SAMPLES_CSV"
            echo "$LIBRARY_ID,REUSE_EXISTING,OK,$OUTPUT_PATH/outs" >> "$ARC_STATUS_CSV"
            count=$((count + 1))
            continue
        fi

        rm -rf "$FASTQ_LIBRARY_DIR"
        mkdir -p "$GEX_FASTQ_DIR" "$ATAC_FASTQ_DIR"

        while IFS=$'\t' read -r RUN_ID SAMPLE_ID ASSAY EXPERIMENT LANE_ID SRA_FILE; do
            LIBRARY_SAMPLE="${SAMPLE_ID}_${ASSAY^^}"

            convert_sra_one_run \
                "$RUN_ID" \
                "$SAMPLE_ID" \
                "$ASSAY" \
                "$LIBRARY_SAMPLE" \
                "$LANE_ID" \
                "$SRA_FILE" \
                "$GEX_FASTQ_DIR"
        done < <(
            awk -F'\t' -v lid="$LIBRARY_ID" '
                NR > 1 && $2 == lid && $3 == "gex" {
                    print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6
                }
            ' "$RUN_MANIFEST"
        )

        while IFS=$'\t' read -r RUN_ID SAMPLE_ID ASSAY EXPERIMENT LANE_ID SRA_FILE; do
            LIBRARY_SAMPLE="${SAMPLE_ID}_${ASSAY^^}"

            convert_sra_one_run \
                "$RUN_ID" \
                "$SAMPLE_ID" \
                "$ASSAY" \
                "$LIBRARY_SAMPLE" \
                "$LANE_ID" \
                "$SRA_FILE" \
                "$ATAC_FASTQ_DIR"
        done < <(
            awk -F'\t' -v lid="$LIBRARY_ID" '
                NR > 1 && $2 == lid && $3 == "atac" {
                    print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6
                }
            ' "$RUN_MANIFEST"
        )

        write_libraries_csv "$LIBRARY_ID" "$FASTQ_LIBRARY_DIR" "$LIBRARIES_CSV"

        if run_cellranger_arc_count "$LIBRARY_ID" "$LIBRARIES_CSV"; then
            echo "  -> Cell Ranger ARC finished. Cleaning FASTQ..."
            echo "$LIBRARY_ID,$OUTPUT_PATH/outs" >> "$SUCCESS_SAMPLES_CSV"
            echo "$LIBRARY_ID,SUCCESS,OK,$OUTPUT_PATH/outs" >> "$ARC_STATUS_CSV"
            rm -rf "$FASTQ_LIBRARY_DIR"
        else
            echo "  -> Warning: Cell Ranger ARC failed"
            echo "$LIBRARY_ID,EXCLUDE,CELLRANGER_ARC_FAILED,$OUTPUT_PATH/outs" >> "$ARC_STATUS_CSV"
            exit 1
        fi

        count=$((count + 1))
    done
fi

# ==========================================
# [PART 3] CellBender on ARC raw H5
# ==========================================

echo "=== [PART 3] Starting CellBender on ARC raw H5 ==="

DOCKER_DIR="$ROOT_DIR/docker"
CB_DOCKERFILE="cellbender.Dockerfile"
CB_IMAGE_NAME="cellbender:pipeline"
EXTRACT_GEX_SCRIPT="$ROOT_DIR/pipeline/utils/single_cell/extract_gex_from_h5.py"

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

CB_SUCCESS_CSV="$CELLBENDER_DIR/success_samples.csv"
CB_STATUS_CSV="$CELLBENDER_DIR/cellbender_status.csv"

echo "LibraryID,FilteredH5" > "$CB_SUCCESS_CSV"
echo "LibraryID,Status,Reason,ExpectedCells,TotalBarcodes,OutputH5" > "$CB_STATUS_CSV"

count=1
while IFS= read -r LIBRARY_ID; do
    echo "[Process: CellBender] $count / $TOTAL_COUNT : $LIBRARY_ID"

    ARC_OUTS="$CELLRANGER_ARC_DIR/$LIBRARY_ID/outs"
    ARC_RAW_H5="$ARC_OUTS/raw_feature_bc_matrix.h5"
    ARC_SUMMARY_CSV="$ARC_OUTS/summary.csv"

    SAMPLE_CB_DIR="$CELLBENDER_DIR/$LIBRARY_ID"
    GEX_RAW_H5="$SAMPLE_CB_DIR/${LIBRARY_ID}_gex_raw.h5"
    GEX_CB_FILTERED="$SAMPLE_CB_DIR/${LIBRARY_ID}_gex_filtered.h5"
    mkdir -p "$SAMPLE_CB_DIR"

    if [[ -s "$GEX_CB_FILTERED" ]]; then
        echo "  -> Existing GEX-only CellBender output found. Reusing."
        echo "$LIBRARY_ID,$GEX_CB_FILTERED" >> "$CB_SUCCESS_CSV"
        echo "$LIBRARY_ID,REUSE_EXISTING,OK,N/A,N/A,$GEX_CB_FILTERED" >> "$CB_STATUS_CSV"
        count=$((count + 1)); continue
    fi

    if [[ ! -s "$ARC_RAW_H5" ]]; then
        echo "$LIBRARY_ID,SKIP,ARC_RAW_INPUT_MISSING,N/A,N/A," >> "$CB_STATUS_CSV"
        count=$((count + 1)); continue
    fi

    METRICS=$("$PIXI_EXEC" run python3 - "$ARC_SUMMARY_CSV" <<'PY'
import csv, os, sys

def as_int(x):
    try:
        return str(int(float(str(x).replace(",", "").strip())))
    except Exception:
        return "NA"

row = {}
if os.path.isfile(sys.argv[1]):
    with open(sys.argv[1], newline="") as f:
        row = next(csv.DictReader(f), {})

lookup = {str(k).strip().lower(): v for k, v in row.items()}
for key in (
    "estimated number of cells - gene expression",
    "estimated number of cells",
):
    value = lookup.get(key)
    if value not in (None, ""):
        print(as_int(value))
        break
else:
    print("NA")
PY
)

    if ! [[ "$METRICS" =~ ^[0-9]+$ ]]; then
        echo "$LIBRARY_ID,SKIP,EXPECTED_CELLS_UNKNOWN,$METRICS,N/A," >> "$CB_STATUS_CSV"
        count=$((count + 1)); continue
    fi

    TOTAL_BARCODES=$("$PIXI_EXEC" run python3 - "$ARC_RAW_H5" <<'PY'
import h5py, sys
with h5py.File(sys.argv[1], "r") as f:
    print(f["matrix/barcodes"].shape[0])
PY
)

    AUTO_TOTAL_DROPLETS=$((METRICS * 3))
    (( AUTO_TOTAL_DROPLETS < 15000 )) && AUTO_TOTAL_DROPLETS=15000
    (( AUTO_TOTAL_DROPLETS > TOTAL_BARCODES )) && AUTO_TOTAL_DROPLETS="$TOTAL_BARCODES"

    "$PIXI_EXEC" run python3 "$EXTRACT_GEX_SCRIPT" \
        --input "$ARC_RAW_H5" \
        --output "$GEX_RAW_H5"

    echo "  -> expected_cells=$METRICS total_droplets=$AUTO_TOTAL_DROPLETS"

    if [[ ! -s "$GEX_CB_FILTERED" ]]; then
        docker run --rm --gpus all \
            -v "$SAMPLE_CB_DIR":/work \
            -w /work \
            -u "$(id -u):$(id -g)" \
            -e MPLCONFIGDIR=/tmp \
            -e HOME=/tmp \
            "$CB_IMAGE_NAME" \
            remove-background \
            --input "/work/${LIBRARY_ID}_gex_raw.h5" \
            --output "/work/${LIBRARY_ID}_gex.h5" \
            --cuda \
            --expected-cells "$METRICS" \
            --total-droplets-included "$AUTO_TOTAL_DROPLETS" \
            --fpr 0.01 \
            --epochs 150 || {
                echo "$LIBRARY_ID,EXCLUDE,CELLBENDER_FAILED,$METRICS,$TOTAL_BARCODES," >> "$CB_STATUS_CSV"
                exit 1
            }
    fi

    [[ -s "$GEX_CB_FILTERED" ]] || {
        echo "$LIBRARY_ID,EXCLUDE,GEX_FILTERED_OUTPUT_MISSING,$METRICS,$TOTAL_BARCODES," >> "$CB_STATUS_CSV"
        exit 1
    }

    echo "$LIBRARY_ID,$GEX_CB_FILTERED" >> "$CB_SUCCESS_CSV"
    echo "$LIBRARY_ID,SUCCESS,OK,$METRICS,$TOTAL_BARCODES,$GEX_CB_FILTERED" >> "$CB_STATUS_CSV"

    count=$((count + 1))
done < "$LIBRARY_LIST"

# ==========================================
# [PART 4] scDblFinder
# ==========================================

echo "=== [PART 4] Starting scDblFinder Doublet Detection ==="

SUMMARY_CSV="$SCDBLFINDER_DIR/summary_scdblfinder.csv"
SCDBLFINDER_R_SCRIPT="$ROOT_DIR/R/scdblfinder.R"
SCDBLFINDER_PYTHON_SCRIPT="$ROOT_DIR/pipeline/utils/single_cell/scdblfinder.py"

SCDBL_COMPLETE=true
while IFS= read -r LIBRARY_ID; do
    [[ -f "$SCDBLFINDER_DIR/$LIBRARY_ID/${LIBRARY_ID}_singlet_barcodes.csv" \
        && -f "$SCDBLFINDER_DIR/$LIBRARY_ID/${LIBRARY_ID}_summary.csv" \
        && -f "$SCDBLFINDER_DIR/$LIBRARY_ID/${LIBRARY_ID}_clean.h5ad" ]] \
        || SCDBL_COMPLETE=false
done < "$LIBRARY_LIST"

if [[ "$SCDBL_COMPLETE" == true ]]; then
    echo "  -> Complete scDblFinder outputs found. Reusing."
else
    "$PIXI_EXEC" run Rscript "$SCDBLFINDER_R_SCRIPT" \
        --manifest "$CB_SUCCESS_CSV" \
        --outdir "$SCDBLFINDER_DIR" \
        --threads 2 \
        --dbr_per1k 0.008 \
        --seed 1

    while IFS=, read -r LIBRARY_ID INPUT_TARGET || [[ -n "${LIBRARY_ID:-}" ]]; do
        [[ -z "${LIBRARY_ID:-}" || "$LIBRARY_ID" == "LibraryID" ]] && continue

        OUTDIR="$SCDBLFINDER_DIR/$LIBRARY_ID"
        SINGLET_CSV="$OUTDIR/${LIBRARY_ID}_singlet_barcodes.csv"

        [[ -f "$INPUT_TARGET" && -f "$SINGLET_CSV" ]] || {
            echo "Error: input or singlet barcode file missing for $LIBRARY_ID"
            exit 1
        }

        "$PIXI_EXEC" run python3 "$SCDBLFINDER_PYTHON_SCRIPT" \
            --input "$INPUT_TARGET" \
            --singlets "$SINGLET_CSV" \
            --outdir "$OUTDIR" \
            --sample_id "$LIBRARY_ID"
    done < "$CB_SUCCESS_CSV"

    echo "LibraryID,TotalCells,RemovedZeroCountCells,PredictedDoublets,PredictedSinglets,ObservedDoubletFraction,ExpectedDoubletFraction,dbr_per1k,scDblFinderVersion,Status,Note" > "$SUMMARY_CSV"
    for f in "$SCDBLFINDER_DIR"/*/*_summary.csv; do
        [[ -f "$f" ]] && tail -n +2 "$f" >> "$SUMMARY_CSV"
    done
fi

# ==========================================
# [PART 5] Merge GEX-only h5ad files
# ==========================================

echo "=== [PART 5] Merging clean GEX-only h5ad files ==="

MERGE_SCRIPT="$ROOT_DIR/pipeline/utils/single_cell/merge_h5ad.py"
H5AD_MATRIX_DIR="$ROOT_DIR/$MERGED_H5AD"
RNA_H5AD="$H5AD_MATRIX_DIR/${PROJECT_ID}_rna.h5ad"

"$PIXI_EXEC" run python3 "$MERGE_SCRIPT" \
    --input_dir "$SCDBLFINDER_DIR" \
    --project_id "$PROJECT_ID" \
    --output_dir "$H5AD_MATRIX_DIR" \
    --output_name "${PROJECT_ID}_rna.h5ad" \
    --library_manifest "$LIBRARY_MANIFEST" \
    --sra_xml "$SRA_XML" \
    --multiome

# ==========================================
# [PART 6] SnapATAC2 de novo peaks + paired SCENIC+ files
# ==========================================

echo "=== [PART 6] Calling de novo peaks with SnapATAC2 ==="
echo $LIBRARY_LIST
SNAPATAC_SCRIPT="$ROOT_DIR/pipeline/utils/single_cell/snapatac.py"
if [[ -f "$SNAPATAC_DIR/_SUCCESS" ]]; then
    echo "  -> Completed SnapATAC2 output found. Reusing."
else
    "$PIXI_EXEC" run python3 "$SNAPATAC_SCRIPT" \
        --cellranger-arc-dir "$CELLRANGER_ARC_DIR" \
        --project-id "$PROJECT_ID" \
        --library-list "$LIBRARY_LIST" \
        --scdblfinder-dir "$SCDBLFINDER_DIR" \
        --rna-h5ad "$RNA_H5AD" \
        --outdir "$SNAPATAC_DIR" \
        --chrom-sizes "$CHROM_SIZES" \
        --gene-annotation "$GENE_ANNOTATION"
fi

# ==========================================
# [PART 7] SCENIC+
# ==========================================

echo "=== [PART 7] Starting pycisTopic & SCENIC+ ==="

SCENICPLUS_DOCKER_IMAGE="scenicplus-local:1.0a2"
SCENICPLUS_DOCKERFILE="$ROOT_DIR/docker/scenicplus.Dockerfile"

SCENICPLUS_INPUT_DIR="$SNAPATAC_DIR/scenic_input"
SCENICPLUS_DIR="$PROJECT_ROOT_DIR/scenicplus"
SCENICPLUS_SCRIPT="$ROOT_DIR/pipeline/utils/single_cell/scenicplus_prepare.py"

echo "=== [Docker Check] Image: $SCENICPLUS_DOCKER_IMAGE ==="
if ! docker image inspect "$SCENICPLUS_DOCKER_IMAGE" >/dev/null 2>&1 
then
    echo "  -> Building local SCENIC+ image..."
    docker buildx build \
        --load \
        --build-arg SCENICPLUS_REF=v1.0a2 \
        -t "$SCENICPLUS_DOCKER_IMAGE" \
        -f "$SCENICPLUS_DOCKERFILE" \
        "$ROOT_DIR" || {
            echo "Error: SCENIC+ Docker build failed!"
            exit 1
        }
    echo "  -> Build complete."
else
    echo "  -> Image found. Skipping build."
fi

docker run --rm \
    -u "$(id -u):$(id -g)" \
    -v "$ROOT_DIR:$ROOT_DIR" \
    -v "$SNAPATAC_DIR:$SNAPATAC_DIR" \
    -v "$PROJECT_ROOT_DIR:$PROJECT_ROOT_DIR" \
    -w "$ROOT_DIR" \
    -e MPLCONFIGDIR=/tmp \
    -e HOME=/tmp \
    -e NUMBA_CACHE_DIR=/tmp \
    -e MALLET_MEMORY=48g \
    "$SCENICPLUS_DOCKER_IMAGE" \
    python3 "$SCENICPLUS_SCRIPT" \
        --input-dir "$SCENICPLUS_INPUT_DIR" \
        --outdir "$SCENICPLUS_DIR" \
        --threads "$N_THREADS" \
        --project-id "$PROJECT_ID" || {
            echo "Error: SCENIC+ execution failed in Docker!"
            exit 1
        }

echo "  -> GEX-only RNA: $RNA_H5AD"
echo "  -> Paired SCENIC+ input: $SCENICPLUS_INPUT_DIR"
echo "Pipeline completed."
