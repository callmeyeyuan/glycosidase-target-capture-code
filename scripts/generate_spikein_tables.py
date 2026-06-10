#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Generate spike-in summary tables for the target-capture dataset.

This script converts a long-format reference-level mapping table into:
  1) spikein_reads_matrix.tsv
  2) spikein_expected_check.tsv
  3) spikein_detection_summary.tsv

Expected input:
  reference_mapped_reads_long.tsv

Required columns in reference_mapped_reads_long.tsv:
  sample
  reference_id
  reference_type
  spikein_index
  mapped_reads

The script assumes that spike-in reference IDs are named spikein_01 ... spikein_20
and that the final FASTA is ordered according to the submission design:
  spikein_01-10 = common spike-ins
  spikein_11-20 = sample-specific spike-ins

Optional:
  --design spikein_addition_design.tsv

Accepted design table formats:
  Format A:
    known_sample_id    expected_spikeins
    F77                1,2,3,4,5,6,7,8,9,10,11

  Format B:
    known_sample_id    common_spikeins    sample_specific_spikein
    F77                1,2,3,4,5,6,7,8,9,10    11

If --design is not supplied, the default 10 positive-sample design used in this
study is applied.

Example:
python generate_spikein_tables.py \
  --reference-mapped-long target_capture_qc/reference_mapped_reads_long.tsv \
  --out-dir target_capture_qc \
  --threshold 5
"""

import argparse
import csv
import re
from pathlib import Path
from typing import Dict, List, Set, Optional

import pandas as pd


DEFAULT_EXPECTED_BY_SAMPLE = {
    "F77":  list(range(1, 11)) + [11],
    "F153": list(range(1, 11)) + [12],
    "F88":  list(range(1, 11)) + [13],
    "F154": list(range(1, 11)) + [14],
    "F135": list(range(1, 11)) + [15],
    "F168": list(range(1, 11)) + [16],
    "F175": list(range(1, 11)) + [17],
    "F48":  list(range(1, 11)) + [18],
    "F49":  list(range(1, 11)) + [19],
    "F179": list(range(1, 11)) + [20],
}


def parse_int_list(value) -> List[int]:
    """
    Parse spike-in index strings.

    Supported examples:
      1,2,3
      spikein01,spikein02
      spikein_01,spikein_02
      spikein01-spikein10
      spikein01-spikein10,spikein11
    """
    if pd.isna(value):
        return []
    s = str(value).strip()
    if not s:
        return []

    out = []
    items = re.split(r"[,;\s]+", s)

    def one_token_to_int(tok: str) -> int:
        tok = tok.strip()
        tok = tok.replace("spikein_", "")
        tok = tok.replace("spikein", "")
        tok = tok.replace("Spikein", "")
        tok = tok.replace("SPIKEIN", "")
        return int(tok)

    for item in items:
        item = item.strip()
        if not item:
            continue

        # handle ranges such as spikein01-spikein10
        if "-" in item:
            left, right = item.split("-", 1)
            start = one_token_to_int(left)
            end = one_token_to_int(right)
            if start <= end:
                out.extend(range(start, end + 1))
            else:
                out.extend(range(start, end - 1, -1))
        else:
            out.append(one_token_to_int(item))

    return out


def load_design(path: Optional[str]) -> Dict[str, List[int]]:
    if not path:
        return DEFAULT_EXPECTED_BY_SAMPLE.copy()

    df = pd.read_csv(path, sep="\t")
    cols = set(df.columns)

    sample_col = None
    for c in ["known_sample_id", "sample_id", "sample"]:
        if c in cols:
            sample_col = c
            break
    if sample_col is None:
        raise SystemExit(
            "Design file must contain one of: known_sample_id, sample_id, sample"
        )

    expected = {}
    if "expected_spikeins" in cols:
        for _, row in df.iterrows():
            sid = str(row[sample_col]).strip()
            vals = parse_int_list(row["expected_spikeins"])
            if sid and vals:
                expected[sid] = vals
    elif "common_spikeins" in cols and "sample_specific_spikein" in cols:
        for _, row in df.iterrows():
            sid = str(row[sample_col]).strip()
            common = parse_int_list(row["common_spikeins"])
            specific = parse_int_list(row["sample_specific_spikein"])
            vals = common + specific
            if sid and vals:
                expected[sid] = vals
    else:
        raise SystemExit(
            "Design file must contain either expected_spikeins or "
            "common_spikeins + sample_specific_spikein."
        )

    return expected


def identify_sample_id(sample_name: str, expected_by_sample: Dict[str, List[int]]) -> str:
    sample_name = str(sample_name)
    for sid in expected_by_sample:
        if re.search(rf"(^|[^A-Za-z0-9]){re.escape(sid)}([^A-Za-z0-9]|$)", sample_name):
            return sid
    m = re.search(r"(F\d+)", sample_name)
    if m and m.group(1) in expected_by_sample:
        return m.group(1)
    return ""


def spike_index_from_reference(row) -> Optional[int]:
    # Prefer explicit spikein_index if available
    if "spikein_index" in row.index and not pd.isna(row["spikein_index"]):
        try:
            return int(float(row["spikein_index"]))
        except Exception:
            pass

    rid = str(row.get("reference_id", ""))
    m = re.search(r"spikein[_-]?(\d+)", rid, flags=re.IGNORECASE)
    if m:
        return int(m.group(1))
    return None


def comma_join(vals) -> str:
    vals = sorted(int(v) for v in vals)
    return ",".join(str(v) for v in vals)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-mapped-long", required=True,
                        help="reference_mapped_reads_long.tsv")
    parser.add_argument("--design", default=None,
                        help="Optional spikein_addition_design.tsv")
    parser.add_argument("--out-dir", required=True,
                        help="Output directory")
    parser.add_argument("--threshold", type=int, default=5,
                        help="Read threshold for detection")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    expected_by_sample = load_design(args.design)

    df = pd.read_csv(args.reference_mapped_long, sep="\t")

    required = {"sample", "reference_id", "mapped_reads"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"Missing required columns in reference_mapped_long: {missing}")

    if "reference_type" in df.columns:
        spike = df[df["reference_type"].astype(str).str.lower().str.contains("spike")].copy()
    else:
        spike = df[df["reference_id"].astype(str).str.contains("spikein", case=False)].copy()

    if spike.empty:
        raise SystemExit("No spike-in rows found in reference_mapped_long.")

    spike["spikein_index"] = spike.apply(spike_index_from_reference, axis=1)
    spike = spike.dropna(subset=["spikein_index"]).copy()
    spike["spikein_index"] = spike["spikein_index"].astype(int)
    spike = spike[(spike["spikein_index"] >= 1) & (spike["spikein_index"] <= 20)]

    spike["spikein_col"] = spike["spikein_index"].map(lambda x: f"spikein_{x:02d}")
    spike["known_sample_id"] = spike["sample"].map(
        lambda x: identify_sample_id(x, expected_by_sample)
    )

    matrix = (
        spike.pivot_table(
            index=["sample", "known_sample_id"],
            columns="spikein_col",
            values="mapped_reads",
            aggfunc="sum",
            fill_value=0,
        )
        .reset_index()
    )

    for i in range(1, 21):
        col = f"spikein_{i:02d}"
        if col not in matrix.columns:
            matrix[col] = 0

    matrix = matrix[["sample", "known_sample_id"] + [f"spikein_{i:02d}" for i in range(1, 21)]]
    matrix_path = out_dir / "spikein_reads_matrix.tsv"
    matrix.to_csv(matrix_path, sep="\t", index=False)

    check_rows = []
    pos_samples = 0
    pos_all_expected = 0
    pos_with_missing = 0
    unexpected_events_pos = 0
    neg_samples = 0
    neg_with_any = 0

    for _, row in matrix.iterrows():
        sample = row["sample"]
        sid = str(row["known_sample_id"]).strip()
        expected = set(expected_by_sample.get(sid, [])) if sid else set()
        detected = set(
            i for i in range(1, 21)
            if int(row[f"spikein_{i:02d}"]) >= args.threshold
        )

        missing = expected - detected
        unexpected = detected - expected

        if sid:
            pos_samples += 1
            if not missing:
                pos_all_expected += 1
            else:
                pos_with_missing += 1
            unexpected_events_pos += len(unexpected)
        else:
            neg_samples += 1
            if detected:
                neg_with_any += 1

        check_rows.append({
            "sample": sample,
            "known_sample_id": sid,
            "expected_spikeins": comma_join(expected),
            "detected_spikeins": comma_join(detected),
            "missing_expected_spikeins": comma_join(missing),
            "unexpected_detected_spikeins": comma_join(unexpected),
            "detect_threshold_reads": args.threshold,
        })

    check_path = out_dir / "spikein_expected_check.tsv"
    pd.DataFrame(check_rows).to_csv(check_path, sep="\t", index=False)

    summary_path = out_dir / "spikein_detection_summary.tsv"
    with open(summary_path, "w", encoding="utf-8", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["metric", "value"])
        w.writerow(["threshold_reads", args.threshold])
        w.writerow(["positive_samples", pos_samples])
        w.writerow(["positive_samples_with_all_expected_detected", pos_all_expected])
        w.writerow(["positive_samples_with_missing_expected", pos_with_missing])
        w.writerow(["total_unexpected_detected_events_in_positive_samples", unexpected_events_pos])
        w.writerow(["negative_samples", neg_samples])
        w.writerow(["negative_samples_with_any_detected_spikein", neg_with_any])

    print(f"Wrote: {matrix_path}")
    print(f"Wrote: {check_path}")
    print(f"Wrote: {summary_path}")


if __name__ == "__main__":
    main()
