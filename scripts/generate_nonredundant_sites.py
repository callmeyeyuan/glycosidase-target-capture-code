#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Generate a cohort-level non-redundant variant table from sample-level VCF files.

This script extracts CHROM, POS, REF and ALT from all VCF/VCF.GZ files in a
directory, removes duplicate records across samples and writes:
  1) nonredundant_sites.tsv
  2) nonredundant_sites_summary.tsv

By default, multi-allelic ALT values are split into separate records.

Example:
python generate_nonredundant_sites.py \
  --vcf-dir merged_vcf \
  --output nonredundant_variants/nonredundant_sites.tsv \
  --summary nonredundant_variants/nonredundant_sites_summary.tsv
"""

import argparse
import gzip
from pathlib import Path
from typing import Iterable, Tuple

import pandas as pd


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def iter_vcf_records(path: Path, split_alt: bool = True) -> Iterable[Tuple[str, int, str, str]]:
    with open_text(path) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            chrom, pos, _vid, ref, alt = parts[:5]
            if not chrom or not pos or not ref or not alt or alt == ".":
                continue
            try:
                pos_int = int(pos)
            except ValueError:
                continue

            if split_alt:
                alts = [a for a in alt.split(",") if a and a != "."]
            else:
                alts = [alt]

            for a in alts:
                yield chrom, pos_int, ref, a


def sample_name_from_vcf(path: Path) -> str:
    name = path.name
    for suffix in [".vcf.gz", ".vcf", ".gz"]:
        if name.endswith(suffix):
            name = name[:-len(suffix)]
    return name


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf-dir", required=True, help="Directory containing sample-level VCF files")
    parser.add_argument("--pattern", default="*.vcf*",
                        help="Glob pattern for VCF files; default: *.vcf*")
    parser.add_argument("--output", required=True,
                        help="Output nonredundant_sites.tsv")
    parser.add_argument("--summary", default=None,
                        help="Optional output summary TSV")
    parser.add_argument("--no-header", action="store_true",
                        help="Write output without header to match older table style")
    parser.add_argument("--no-split-alt", action="store_true",
                        help="Do not split comma-separated ALT alleles")
    args = parser.parse_args()

    vcf_dir = Path(args.vcf_dir)
    files = sorted(vcf_dir.glob(args.pattern))
    files = [p for p in files if p.is_file() and (str(p).endswith(".vcf") or str(p).endswith(".vcf.gz"))]

    if not files:
        raise SystemExit(f"No VCF files found in {vcf_dir} with pattern {args.pattern}")

    unique_records = set()
    total_records = 0
    per_sample_counts = []

    for fp in files:
        sample_count = 0
        for rec in iter_vcf_records(fp, split_alt=not args.no_split_alt):
            unique_records.add(rec)
            total_records += 1
            sample_count += 1
        per_sample_counts.append((sample_name_from_vcf(fp), sample_count))

    records = sorted(unique_records, key=lambda x: (x[0], x[1], x[2], x[3]))
    out_df = pd.DataFrame(records, columns=["CHROM", "POS", "REF", "ALT"])

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(output, sep="\t", index=False, header=not args.no_header)

    if args.summary:
        summary = {
            "vcf_files": len(files),
            "total_sample_level_variant_records": total_records,
            "nonredundant_variant_records": len(out_df),
            "unique_variant_positions_CHROM_POS": out_df[["CHROM", "POS"]].drop_duplicates().shape[0],
            "target_references_with_variants": out_df["CHROM"].nunique(),
            "min_records_per_sample": min(c for _, c in per_sample_counts),
            "max_records_per_sample": max(c for _, c in per_sample_counts),
        }
        summary_path = Path(args.summary)
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(
            [{"metric": k, "value": v} for k, v in summary.items()]
        ).to_csv(summary_path, sep="\t", index=False)

    print(f"VCF files: {len(files)}")
    print(f"Total sample-level variant records: {total_records}")
    print(f"Non-redundant variant records: {len(out_df)}")
    print(f"Wrote: {output}")
    if args.summary:
        print(f"Wrote: {args.summary}")


if __name__ == "__main__":
    main()
