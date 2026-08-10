#!/usr/bin/env bash
# =============================================================================
# Refined RNA-seq Preprocessing Pipeline
#
# Combines Destiny and Micaiah's preprocessing scripts
#
# AI Disclaimer: AI tools were used to assist with code commenting and error handling.
# All logic was reviewed by the authors.
#
# Usage:
#   Single sample test:  bash rnaseq_pipeline_refined.sh --test
#   Full run (all):      bash rnaseq_pipeline_refined.sh
# =============================================================================

set -euo pipefail   # exit on error, undefined var, or pipe failure

# =============================================================================
# CONFIGURATION — edit these paths before running
# =============================================================================

BASE_DIR="$(pwd)"                          # set to your project root

# directories will be created if they do not exist
DATA_DIR="${BASE_DIR}/data"                 # data Folder
FASTQ_DIR="${DATA_DIR}/fastq_files"         #  raw Fastq files
TRIM_DIR="${DATA_DIR}/trimmed_fastq"        # trimmed files
BAM_DIR="${BASE_DIR}/bam"                   # BAM files
QUANT_DIR="${BASE_DIR}/quants"              # count files
LOG_DIR="${BASE_DIR}/logs"                  # log files
ANNOTATION_DIR="${BASE_DIR}/annotation"     # annotation files
HISAT2_IDX="${BASE_DIR}/HISAT2/hg38/genome"   # path to HISAT2 genome index
GTF="${ANNOTATION_DIR}/Homo_sapiens.GRCh38.115.gtf" # path to gtf

# path to adapters; will be unique to your environment
ADAPTERS="/opt/anaconda3/envs/rnaseq_env/share/trimmomatic/adapters/TruSeq3-PE.fa" 


THREADS=6                                       # Specific to your device; 4 commonly used

# All 10 SRR accessions
ALL_SAMPLES=(
    SRR31443248
    SRR31443249
    SRR31443250
    SRR31443251
    SRR31443252
    SRR31443253
    SRR31443254
    SRR31443255
    SRR31443256
    SRR31443257
)

# =============================================================================
# PARSE FLAGS
# --test runs only the first sample so you can validate the pipeline quickly
# =============================================================================

TEST_MODE=false
if [[ "${1:-}" == "--test" ]]; then
    TEST_MODE=true
    SAMPLES=("${ALL_SAMPLES[0]}")
    echo ">>> TEST MODE: running on ${SAMPLES[0]} only <<<"
else
    SAMPLES=("${ALL_SAMPLES[@]}")
    echo ">>> FULL MODE: running on ${#SAMPLES[@]} samples <<<"
fi

# =============================================================================
# SETUP — create all output directories upfront
# =============================================================================

mkdir -p \
    "$FASTQ_DIR" \
    "$TRIM_DIR" \
    "$BAM_DIR" \
    "$QUANT_DIR" \
    "$LOG_DIR" \
    "$ANNOTATION_DIR" \
    data/qc_files/fastqc \
    data/qc_files/multiqc \
    data/qc_files/post_trim_fastqc \
    data/qc_files/post_trim_multiqc \
    HISAT2

# =============================================================================
# STEP 0: Download raw FASTQ files from SRA 
# SRA Toolkit must already be installed
# (SKIP THIS STEP!! if files already exist)
# =============================================================================

echo ""
echo "===== STEP 0: Downloading FASTQ files from SRA ====="

cd "$FASTQ_DIR"
for SAMPLE in "${SAMPLES[@]}"; do
    R1="${FASTQ_DIR}/${SAMPLE}_1.fastq.gz"
    if [[ -f "$R1" ]]; then
        echo "  [skip] ${SAMPLE} already downloaded"
    else
        echo "  Downloading ${SAMPLE}..."
        fastq-dump --split-files --gzip "$SAMPLE"
    fi
done
cd "$BASE_DIR"

# =============================================================================
# STEP 1: Raw QC — FastQC + MultiQC
# =============================================================================

echo ""
echo "===== STEP 1: Raw QC (FastQC + MultiQC) ====="

# Build a list of all R1+R2 files for the samples we're running
RAW_FASTQS=()
for SAMPLE in "${SAMPLES[@]}"; do
    RAW_FASTQS+=( "${FASTQ_DIR}/${SAMPLE}_1.fastq.gz" "${FASTQ_DIR}/${SAMPLE}_2.fastq.gz" )
done

fastqc \
    "${RAW_FASTQS[@]}" \
    -o data/qc_files/fastqc/ \
    -t "$THREADS"

multiqc data/qc_files/fastqc/ \
    -o data/qc_files/multiqc/ \
    --force

# =============================================================================
# STEP 2: Trimming — Trimmomatic (paired-end)
# + post-trim FastQC per sample
# =============================================================================

echo ""
echo "===== STEP 2: Trimmomatic + post-trim FastQC ====="

for SAMPLE in "${SAMPLES[@]}"; do
    R1="${FASTQ_DIR}/${SAMPLE}_1.fastq.gz"
    R2="${FASTQ_DIR}/${SAMPLE}_2.fastq.gz"

    OUT_R1P="${TRIM_DIR}/${SAMPLE}_1_paired.fastq.gz"
    OUT_R1U="${TRIM_DIR}/${SAMPLE}_1_unpaired.fastq.gz"
    OUT_R2P="${TRIM_DIR}/${SAMPLE}_2_paired.fastq.gz"
    OUT_R2U="${TRIM_DIR}/${SAMPLE}_2_unpaired.fastq.gz"
    TRIM_LOG="${LOG_DIR}/${SAMPLE}_trimmomatic.log"

    echo "  Trimming ${SAMPLE}..."

    trimmomatic PE \
        -threads "$THREADS" \
        -phred33 \
        "$R1" "$R2" \
        "$OUT_R1P" "$OUT_R1U" \
        "$OUT_R2P" "$OUT_R2U" \
        ILLUMINACLIP:"${ADAPTERS}":2:30:10 \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:20 \
        MINLEN:36 \
        2> "$TRIM_LOG"           

    echo "  Post-trim FastQC for ${SAMPLE}..."

    fastqc \
        -t "$THREADS" \
        -o data/qc_files/post_trim_fastqc/ \
        "$OUT_R1P" \
        "$OUT_R2P"
done

# =============================================================================
# STEP 3: Post-trim MultiQC 
# =============================================================================

echo ""
echo "===== STEP 3: Post-trim MultiQC ====="

multiqc data/qc_files/post_trim_fastqc/ \
    -o data/qc_files/post_trim_multiqc/ \
    --force

# =============================================================================
# STEP 4: Alignment — HISAT2 → samtools sort + index
# =============================================================================

echo ""
echo "===== STEP 4: Alignment (HISAT2 + samtools) ====="

# Download genome index if missing
if [[ ! -f "${HISAT2_IDX}.1.ht2" ]]; then
    echo "  HISAT2 index not found — downloading hg38..."
    wget -q https://genome-idx.s3.amazonaws.com/hisat/hg38_genome.tar.gz -P HISAT2/
    tar -xzf HISAT2/hg38_genome.tar.gz -C HISAT2/
    rm HISAT2/hg38_genome.tar.gz
fi

for SAMPLE in "${SAMPLES[@]}"; do
    R1="${TRIM_DIR}/${SAMPLE}_1_paired.fastq.gz"
    R2="${TRIM_DIR}/${SAMPLE}_2_paired.fastq.gz"
    BAM="${BAM_DIR}/${SAMPLE}.bam"
    ALIGN_LOG="${LOG_DIR}/${SAMPLE}_hisat2.log"

    echo "  Aligning ${SAMPLE}..."

    hisat2 \
        -p "$THREADS" \
        -q \
        -x "$HISAT2_IDX" \
        -1 "$R1" \
        -2 "$R2" \
        2> "$ALIGN_LOG" \
    | samtools sort \
        -@ "$THREADS" \
        -T "${BAM_DIR}/${SAMPLE}.sorttmp" \
        -o "$BAM"

    samtools index "$BAM"      

    echo "  Done: ${BAM}"
done

# =============================================================================
# STEP 5: Quantification — featureCounts
# (skipped in test mode since a single-sample count matrix is usually not needed
#  until all samples are aligned; set SKIP_QUANT=false to override)
# =============================================================================

SKIP_QUANT=false
if [[ "$TEST_MODE" == true ]]; then
    SKIP_QUANT=true
    echo ""
    echo "===== STEP 5: Quantification — SKIPPED in test mode ====="
fi

if [[ "$SKIP_QUANT" == false ]]; then
    echo ""
    echo "===== STEP 5: Quantification (featureCounts) ====="

    # Download GTF if missing
    if [[ ! -f "$GTF" ]]; then
        echo "  GTF not found — downloading..."
        wget -q \
            https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz \
            -P "$ANNOTATION_DIR/"
        gunzip "${ANNOTATION_DIR}/Homo_sapiens.GRCh38.115.gtf.gz"
    fi

    featureCounts \
        -T "$THREADS" \
        -p \
        -B \
        -C \
        -a "$GTF" \
        -o "${QUANT_DIR}/featurecounts_counts.txt" \
        "${BAM_DIR}"/*.bam \
        > "${LOG_DIR}/featurecounts.log" 2>&1

    #Check Counts

    echo ""
    echo "  Total genes detected:"
    grep -v '^#' "${QUANT_DIR}/featurecounts_counts.txt" | tail -n +2 | wc -l

    echo ""
    echo "  Preview (first 5 rows):"
    grep -v '^#' "${QUANT_DIR}/featurecounts_counts.txt" | head -n 6
fi

# =============================================================================
# DONE
# =============================================================================

echo ""
echo "===== Pipeline complete ====="
echo "  QC reports  : data/qc_files/"
echo "  Trimmed     : ${TRIM_DIR}/"
echo "  BAM files   : ${BAM_DIR}/"
echo "  Counts      : ${QUANT_DIR}/"
echo "  Logs        : ${LOG_DIR}/"
