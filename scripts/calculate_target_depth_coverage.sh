#!/usr/bin/env bash
# calculate_target_depth_coverage.sh
#
# Calculate mean target depth and target coverage breadth for a targeted-capture
# reference in which BAM reference contigs are named gene_0001, gene_0002, ...
#
# Assumptions used in this project:
#   - gene_0001 to gene_2337 are glycoside hydrolase target sequences.
#   - gene_2338 to gene_2357 are synthetic spike-in sequences.
#   - BAM files are sorted BAM files aligned to the combined target reference.
#
# Outputs:
#   1. gh2337_targets_full.bed
#      BED intervals for the 2,337 glycoside hydrolase target contigs.
#   2. combined_reference_full.bed
#      BED intervals for all reference contigs in the BAM.
#   3. target_depth_coverage_summary.tsv
#      Per-sample mean target depth and coverage breadth at >=1x, >=5x and >=10x.
#   4. target_depth_coverage_cohort_summary.tsv
#      Cohort-level mean, median, minimum and maximum values.
#
# Example:
#   bash scripts/calculate_target_depth_coverage.sh \
#     --bam-dir /data/ch-w/糖苷酶/target_capture_qc_out_submission_order_v2/bam \
#     --out-dir /data/ch-w/glycosidase_omix_upload/processed_data/target_capture_qc
#
# Requirements:
#   samtools
#   python3 with pandas is optional; if unavailable, the per-sample table is still generated.

set -euo pipefail

BAM_DIR=""
OUT_DIR=""
TARGET_COUNT=2337
BAM_GLOB="*.sorted.bam"

usage() {
    cat <<EOF
Usage:
  bash calculate_target_depth_coverage.sh --bam-dir BAM_DIR --out-dir OUT_DIR [options]

Required:
  --bam-dir DIR       Directory containing sorted BAM files.
  --out-dir DIR       Output directory.

Optional:
  --target-count N    Number of glycoside hydrolase target contigs. Default: 2337.
  --bam-glob PATTERN  BAM filename pattern inside BAM_DIR. Default: *.sorted.bam.
  -h, --help          Show this help message.

Outputs:
  OUT_DIR/gh2337_targets_full.bed
  OUT_DIR/combined_reference_full.bed
  OUT_DIR/target_depth_coverage_summary.tsv
  OUT_DIR/target_depth_coverage_cohort_summary.tsv
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --bam-dir)
            BAM_DIR="$2"
            shift 2
            ;;
        --out-dir)
            OUT_DIR="$2"
            shift 2
            ;;
        --target-count)
            TARGET_COUNT="$2"
            shift 2
            ;;
        --bam-glob)
            BAM_GLOB="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "ERROR: Unknown argument: $1" >&2
            usage
            exit 1
            ;;
    esac
done

if [[ -z "$BAM_DIR" || -z "$OUT_DIR" ]]; then
    echo "ERROR: --bam-dir and --out-dir are required." >&2
    usage
    exit 1
fi

if [[ ! -d "$BAM_DIR" ]]; then
    echo "ERROR: BAM directory does not exist: $BAM_DIR" >&2
    exit 1
fi

if ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: samtools is not available in PATH." >&2
    exit 1
fi

mkdir -p "$OUT_DIR"

shopt -s nullglob
BAM_FILES=( "$BAM_DIR"/$BAM_GLOB )
shopt -u nullglob

if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No BAM files found using pattern: $BAM_DIR/$BAM_GLOB" >&2
    exit 1
fi

FIRST_BAM="${BAM_FILES[0]}"

GH_BED="$OUT_DIR/gh${TARGET_COUNT}_targets_full.bed"
COMBINED_BED="$OUT_DIR/combined_reference_full.bed"
OUT_TSV="$OUT_DIR/target_depth_coverage_summary.tsv"
COHORT_TSV="$OUT_DIR/target_depth_coverage_cohort_summary.tsv"

echo "BAM directory: $BAM_DIR"
echo "BAM files found: ${#BAM_FILES[@]}"
echo "First BAM: $FIRST_BAM"
echo "Output directory: $OUT_DIR"
echo "Target contig count: $TARGET_COUNT"
echo

echo "Generating BED files from BAM idxstats..."

samtools idxstats "$FIRST_BAM" | \
awk 'BEGIN{OFS="\t"} $1!="*" {print $1,0,$2}' \
> "$COMBINED_BED"

samtools idxstats "$FIRST_BAM" | \
awk -v target_count="$TARGET_COUNT" 'BEGIN{OFS="\t"}
$1 ~ /^gene_/ {
    n=substr($1,6)+0
    if(n>=1 && n<=target_count) print $1,0,$2
}' \
> "$GH_BED"

combined_n=$(wc -l < "$COMBINED_BED" | tr -d ' ')
gh_n=$(wc -l < "$GH_BED" | tr -d ' ')

echo "Combined BED: $COMBINED_BED ($combined_n intervals)"
echo "GH target BED: $GH_BED ($gh_n intervals)"

if [[ "$gh_n" -ne "$TARGET_COUNT" ]]; then
    echo "WARNING: GH target BED contains $gh_n intervals, expected $TARGET_COUNT." >&2
    echo "         Please check whether BAM reference names follow gene_0001... naming." >&2
fi

echo
echo "Calculating per-sample target depth and coverage breadth..."

printf "sample_id\ttarget_bases_total\tmean_target_depth\ttarget_bases_covered_1x\ttarget_bases_covered_5x\ttarget_bases_covered_10x\tcoverage_breadth_1x\tcoverage_breadth_5x\tcoverage_breadth_10x\n" > "$OUT_TSV"

for BAM in "${BAM_FILES[@]}"; do
    file=$(basename "$BAM")

    sample=$(echo "$file" | sed -E 's/.*-(F[0-9]+)-.*/\1/')
    if [[ "$sample" == "$file" ]]; then
        sample="$file"
        sample="${sample%.bam}"
        sample="${sample%.sorted}"
        sample="${sample%.dedup}"
        sample="${sample%.rmdup}"
    fi

    echo "Processing $sample"

    if [[ ! -f "${BAM}.bai" ]]; then
        samtools index "$BAM"
    fi

    samtools depth -a -b "$GH_BED" "$BAM" | \
    awk -v sample="$sample" 'BEGIN{
        n=0; sum=0; c1=0; c5=0; c10=0
    }
    {
        d=$3
        n++
        sum+=d
        if(d>=1)c1++
        if(d>=5)c5++
        if(d>=10)c10++
    }
    END{
        if(n>0){
            print sample"\t"n"\t"sum/n"\t"c1"\t"c5"\t"c10"\t"c1/n"\t"c5/n"\t"c10/n
        }else{
            print sample"\t0\tNA\t0\t0\t0\tNA\tNA\tNA"
        }
    }' >> "$OUT_TSV"
done

echo
echo "Per-sample output written to: $OUT_TSV"
echo "Rows:"
wc -l "$OUT_TSV"

echo
echo "Generating cohort-level summary..."

if command -v python3 >/dev/null 2>&1; then
python3 - <<PY
from pathlib import Path
import sys

try:
    import pandas as pd
except Exception:
    print("WARNING: pandas is not available; cohort summary table was not generated.", file=sys.stderr)
    sys.exit(0)

out_tsv = Path("$OUT_TSV")
cohort_tsv = Path("$COHORT_TSV")

df = pd.read_csv(out_tsv, sep="\t")
cols = [
    "mean_target_depth",
    "coverage_breadth_1x",
    "coverage_breadth_5x",
    "coverage_breadth_10x",
]

rows = []
for col in cols:
    rows.append({
        "metric": col,
        "mean": df[col].mean(),
        "median": df[col].median(),
        "min": df[col].min(),
        "max": df[col].max(),
    })

summary = pd.DataFrame(rows)
summary.to_csv(cohort_tsv, sep="\t", index=False)
print(summary.to_string(index=False))
print(f"Cohort summary written to: {cohort_tsv}")
PY
else
    echo "WARNING: python3 is not available; cohort summary table was not generated." >&2
fi

echo
echo "Done."
