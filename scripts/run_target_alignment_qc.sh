#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   bash run_target_alignment_qc.sh /data/ch-w/final_data /data/ch-w/genefasta2337.txt /data/ch-w/shuffled_20.fasta /data/ch-w/target_capture_qc_out 8
#
# Arguments:
#   1. FASTQ root directory
#   2. genefasta2337.txt
#   3. shuffled_20.fasta
#   4. output directory
#   5. number of threads

FASTQ_ROOT="${1:-/data/ch-w/final_data}"
GENE_FASTA="${2:-/data/ch-w/genefasta2337.txt}"
SPIKE_FASTA="${3:-/data/ch-w/shuffled_20.fasta}"
OUTDIR="${4:-/data/ch-w/target_capture_qc_out}"
THREADS="${5:-8}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "${OUTDIR}"/{ref,bam,stats,summary,logs}

echo "[1/5] Checking required tools..."
command -v python3 >/dev/null 2>&1 || { echo "ERROR: python3 not found"; exit 1; }
command -v bwa >/dev/null 2>&1 || { echo "ERROR: bwa not found. Install with: conda install -c bioconda bwa"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not found. Install with: conda install -c bioconda samtools"; exit 1; }

echo "[2/5] Creating sanitized combined reference..."
python3 "${SCRIPT_DIR}/sanitize_reference.py" \
  --genes "${GENE_FASTA}" \
  --spikeins "${SPIKE_FASTA}" \
  --out-fasta "${OUTDIR}/ref/glycosidase2337_plus_spike20.sanitized.fa" \
  --out-map "${OUTDIR}/ref/reference_id_mapping.tsv"

echo "[3/5] Indexing reference with BWA and samtools..."
bwa index "${OUTDIR}/ref/glycosidase2337_plus_spike20.sanitized.fa" \
  > "${OUTDIR}/logs/bwa_index.log" 2>&1

samtools faidx "${OUTDIR}/ref/glycosidase2337_plus_spike20.sanitized.fa"

echo "[4/5] Aligning paired FASTQ files..."
R1_LIST="${OUTDIR}/logs/r1_files.txt"
find "${FASTQ_ROOT}" -type f \( -name "*_R1.fastq.gz" -o -name "*_R1.fq.gz" -o -name "*_1.fastq.gz" -o -name "*_1.fq.gz" \) | sort > "${R1_LIST}"

N=$(wc -l < "${R1_LIST}" | tr -d ' ')
echo "Found ${N} R1 files."

if [[ "${N}" -eq 0 ]]; then
  echo "ERROR: no R1 FASTQ files found under ${FASTQ_ROOT}"
  exit 1
fi

while read -r R1; do
  R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"
  R2="${R2/_R1.fq.gz/_R2.fq.gz}"
  R2="${R2/_1.fastq.gz/_2.fastq.gz}"
  R2="${R2/_1.fq.gz/_2.fq.gz}"

  if [[ ! -f "${R2}" ]]; then
    echo "WARNING: R2 file not found for ${R1}; skipping." | tee -a "${OUTDIR}/logs/missing_pairs.log"
    continue
  fi

  BASENAME="$(basename "${R1}")"
  SAMPLE="${BASENAME%_R1.fastq.gz}"
  SAMPLE="${SAMPLE%_R1.fq.gz}"
  SAMPLE="${SAMPLE%_1.fastq.gz}"
  SAMPLE="${SAMPLE%_1.fq.gz}"

  BAM="${OUTDIR}/bam/${SAMPLE}.sorted.bam"

  if [[ -f "${BAM}" && -f "${BAM}.bai" ]]; then
    echo "Already done: ${SAMPLE}"
    continue
  fi

  echo "Aligning ${SAMPLE}"

  bwa mem -t "${THREADS}" "${OUTDIR}/ref/glycosidase2337_plus_spike20.sanitized.fa" "${R1}" "${R2}" \
    2> "${OUTDIR}/logs/${SAMPLE}.bwa.log" \
    | samtools sort -@ "${THREADS}" -o "${BAM}" -

  samtools index "${BAM}"
  samtools flagstat "${BAM}" > "${OUTDIR}/stats/${SAMPLE}.flagstat.txt"
  samtools idxstats "${BAM}" > "${OUTDIR}/stats/${SAMPLE}.idxstats.tsv"

done < "${R1_LIST}"

echo "[5/5] Summarizing mapping and spike-in detection..."
python3 "${SCRIPT_DIR}/summarize_target_alignment.py" \
  --stats-dir "${OUTDIR}/stats" \
  --out-dir "${OUTDIR}/summary" \
  --detect-threshold 5

echo "Done."
echo "Main outputs:"
echo "  ${OUTDIR}/summary/sample_mapping_summary.tsv"
echo "  ${OUTDIR}/summary/spikein_reads_matrix.tsv"
echo "  ${OUTDIR}/summary/spikein_expected_check.tsv"
echo "  ${OUTDIR}/summary/reference_mapped_reads_long.tsv"
