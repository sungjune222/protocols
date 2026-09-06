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
PROJECT="il_11" 
REFERENCE_GENOME="GRCm39" 
REFERENCE_GTF="GCF_000001635.27_GRCm39_genomic.gtf"
# ==========================================

ROOT_PROJECT_DIR="$DATA_DIR/$PROJECT"
SAMPLE_LIST="$ROOT_PROJECT_DIR/sample_list.txt"
FASTQ_DIR="$ROOT_PROJECT_DIR/fastq"

# ==========================================
# [PART 1] Sample Detection
# ==========================================

echo "=== [PART 1] FASTQ Sample Detection ==="
find "$FASTQ_DIR" -maxdepth 1 -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" \) \
    | sed 's|.*/||' \
    | sed -E 's/\.fastq\.gz$//' \
    | sed -E 's/_[12]$//' \
    | sort -u > "$SAMPLE_LIST"

TOTAL_COUNT=$(wc -l < "$SAMPLE_LIST")
echo "Found $TOTAL_COUNT local samples"

while IFS= read -r SAMPLE_ID; do
    R1="$FASTQ_DIR/${SAMPLE_ID}_1.fastq.gz"
    R2="$FASTQ_DIR/${SAMPLE_ID}_2.fastq.gz"

    if [[ ! -f "$R1" ]]; then
        echo "Error: missing R1 for $SAMPLE_ID -> $R1"
        exit 1
    fi
    if [[ ! -f "$R2" ]]; then
        echo "Error: missing R2 for $SAMPLE_ID -> $R2"
        exit 1
    fi
done < "$SAMPLE_LIST"
echo "All paired FASTQ files are present"

# ==========================================
# [PART 2] STAR Alignment
# ==========================================

echo "=== [PART 2] STAR Alignment from FASTQ ==="
STAR_DIR="$ROOT_PROJECT_DIR/$STAR"
mkdir -p "$STAR_DIR"

BAMS=()
count=1
while IFS= read -r SAMPLE_ID; do
    echo "[Process] $count / $TOTAL_COUNT : $SAMPLE_ID"
    R1="$FASTQ_DIR/${SAMPLE_ID}_1.fastq.gz"
    R2="$FASTQ_DIR/${SAMPLE_ID}_2.fastq.gz"
    
    SAMPLE_DIR="$STAR_DIR/$SAMPLE_ID"
    mkdir -p "$SAMPLE_DIR"
    FINAL_BAM="$SAMPLE_DIR/${SAMPLE_ID}.bam"

    if [[ -s "$FINAL_BAM" ]]; then
        echo "Skipping $SAMPLE_ID (Already exists)"
        BAMS+=("$FINAL_BAM")
        count=$((count + 1)); continue
    fi

    STAR \
        --runThreadN "$N_THREADS" \
        --genomeDir "$ROOT_DIR/references/bulk/$REFERENCE_GENOME" \
        --readFilesIn "$R1" "$R2" \
        --readFilesCommand zcat \
        --twopassMode Basic \
        --outFileNamePrefix "$SAMPLE_DIR/" \
        --outSAMtype BAM SortedByCoordinate 

    mv "$SAMPLE_DIR/Aligned.sortedByCoord.out.bam" "$FINAL_BAM"
    samtools index "$FINAL_BAM"
    samtools flagstat "$FINAL_BAM" > "$SAMPLE_DIR/${SAMPLE_ID}.flagstat.txt"

    BAMS+=("$FINAL_BAM")
    count=$((count + 1))
done < "$SAMPLE_LIST"

# ==========================================
# [PART 3] featureCounts
# ==========================================

echo "=== [PART 3] Running featureCounts ==="
TARGET_FEATURES=("exon" "gene")

for FEATURE in "${TARGET_FEATURES[@]}"; do 
    echo " -> Processing feature : [$FEATURE]"
    FEATURE_TSV="$ROOT_PROJECT_DIR/$FEATURE.tsv"

    if [[ -s "$FEATURE_TSV" && -s "$FEATURE_TSV.summary" ]]; then
        echo "Skipping $FEATURE (Already exists)"
        continue
    fi

    # Paired-end mode
    featureCounts \
        -T "$N_THREADS" \
        -a "$ROOT_DIR/references/raw/$REFERENCE_GTF" \
        -o "$FEATURE_TSV" \
        -t "$FEATURE" \
        -g gene_id \
        -s 0 \
        -p \
        --countReadPairs \
        -B \
        -C \
        "${BAMS[@]}"
done

# ==========================================
# [PART 4] Sample QC
# ==========================================

echo "=== [PART 4] Sample Quality Control ==="
QC_SCRIPT="$ROOT_DIR/pipeline/utils/bulk_qc.py"
QC_DIR="$ROOT_PROJECT_DIR/qc"

"$PIXI_EXEC" run python3 "$QC_SCRIPT" \
    --counts "$ROOT_PROJECT_DIR/exon.tsv" \
    --gtf "$ROOT_DIR/references/raw/$REFERENCE_GTF" \
    --output_dir "$QC_DIR" \
    --project "$PROJECT"

echo "=== DONE ==="