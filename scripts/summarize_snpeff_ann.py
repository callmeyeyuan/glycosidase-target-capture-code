#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import csv
import gzip
from pathlib import Path
from collections import Counter

IMPACT_SCORE = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1, "": 0}

def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")

def get_ann(info):
    for item in info.split(";"):
        if item.startswith("ANN="):
            return item[4:]
    return ""

def parse_ann(info):
    ann = get_ann(info)
    if not ann:
        return []
    out = []
    for raw in ann.split(","):
        f = raw.split("|")
        out.append({
            "annotation": f[1] if len(f) > 1 else "",
            "impact": f[2] if len(f) > 2 else "",
            "gene": f[3] if len(f) > 3 else "",
            "raw": raw
        })
    return out

def primary(entries):
    if not entries:
        return None
    return sorted(entries, key=lambda x: IMPACT_SCORE.get(x["impact"], 0), reverse=True)[0]

def simple_category(annotation):
    parts = set(annotation.split("&"))
    if "stop_gained" in parts:
        return "stop_gained"
    if "missense_variant" in parts:
        return "missense_variant"
    if "synonymous_variant" in parts:
        return "synonymous_variant"
    if "frameshift_variant" in parts:
        return "frameshift_variant"
    if "start_lost" in parts:
        return "start_lost"
    if "stop_lost" in parts:
        return "stop_lost"
    if "upstream_gene_variant" in parts:
        return "upstream_gene_variant"
    if "downstream_gene_variant" in parts:
        return "downstream_gene_variant"
    return annotation if annotation else "no_annotation"

def summarize_vcf(path):
    sample = Path(path).name
    for suf in [".ann.vcf.gz", ".ann.vcf"]:
        if sample.endswith(suf):
            sample = sample[:-len(suf)]
    total = 0
    with_ann = 0
    eff = Counter()
    impact = Counter()
    category = Counter()
    all_eff = Counter()
    gene = Counter()
    warnings = Counter()
    with open_text(path) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 8:
                continue
            total += 1
            entries = parse_ann(p[7])
            if entries:
                with_ann += 1
            pr = primary(entries)
            if pr:
                eff[pr["annotation"]] += 1
                impact[pr["impact"]] += 1
                category[simple_category(pr["annotation"])] += 1
                if pr["gene"]:
                    gene[pr["gene"]] += 1
            for e in entries:
                all_eff[e["annotation"]] += 1
                raw = e["raw"]
                if "WARNING_" in raw or "ERROR_" in raw:
                    for token in raw.split("|")[-1].split("&"):
                        if token.startswith("WARNING_") or token.startswith("ERROR_"):
                            warnings[token] += 1
    return sample, total, with_ann, eff, impact, category, all_eff, gene, warnings

def write_counter(path, header, counter):
    with open(path, "w", encoding="utf-8", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(header)
        for k, v in counter.most_common():
            w.writerow([k, v])

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ann-dir", default="/data/ch-w/callSNP/snpEff_result")
    ap.add_argument("--out-dir", default="/data/ch-w/callSNP/snpEff_summary")
    args = ap.parse_args()
    ann_dir = Path(args.ann_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    files = sorted(list(ann_dir.glob("*.ann.vcf")) + list(ann_dir.glob("*.ann.vcf.gz")))
    if not files:
        raise SystemExit(f"No annotated VCF files found in {ann_dir}")
    total_eff = Counter()
    total_impact = Counter()
    total_category = Counter()
    total_all_eff = Counter()
    total_gene = Counter()
    total_warn = Counter()
    overview_rows = []
    for fp in files:
        sample, total, with_ann, eff, impact, category, all_eff, gene, warnings = summarize_vcf(fp)
        overview_rows.append([sample, total, with_ann, round(with_ann/total, 6) if total else 0])
        total_eff.update(eff)
        total_impact.update(impact)
        total_category.update(category)
        total_all_eff.update(all_eff)
        total_gene.update(gene)
        total_warn.update(warnings)
    with open(out_dir/"snpeff_sample_overview.tsv", "w", encoding="utf-8", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["sample", "total_variant_records", "records_with_ANN", "annotation_rate"])
        w.writerows(overview_rows)
    write_counter(out_dir/"snpeff_primary_effect_total.tsv", ["effect", "count"], total_eff)
    write_counter(out_dir/"snpeff_primary_impact_total.tsv", ["impact", "count"], total_impact)
    write_counter(out_dir/"snpeff_primary_category_total.tsv", ["category", "count"], total_category)
    write_counter(out_dir/"snpeff_all_effect_total.tsv", ["effect", "count"], total_all_eff)
    write_counter(out_dir/"snpeff_gene_total.tsv", ["gene", "count"], total_gene)
    write_counter(out_dir/"snpeff_warning_total.tsv", ["warning_or_error", "count"], total_warn)
    print(f"Annotated VCF files: {len(files)}")
    print(f"Output directory: {out_dir}")

if __name__ == "__main__":
    main()
