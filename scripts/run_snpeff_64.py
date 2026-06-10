#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import multiprocessing
from pathlib import Path

JAVA = "/data/ch-w/java/jdk1.8.0_381/bin/java"
SNPEFF = "/data/ch-w/SnpEff/snpEff/snpEff.jar"
DATA_DIR = "/data/ch-w/SnpEff/snpEff/data"
IN_DIR = "/data/ch-w/callSNP/merge"
OUT_DIR = "/data/ch-w/callSNP/snpEff_result"
LIST_FILE = "/data/ch-w/callSNP/list64.txt"
THREADS = 4

Path(OUT_DIR).mkdir(parents=True, exist_ok=True)

def run_one(sample):
    sample = sample.strip()
    if not sample:
        return "empty sample"
    in_vcf = f"{IN_DIR}/{sample}.unann.vcf"
    out_vcf = f"{OUT_DIR}/{sample}.ann.vcf"
    log_file = f"{OUT_DIR}/{sample}.snpeff.log"
    if not os.path.exists(in_vcf):
        return f"{sample}: MISSING input"
    if os.path.exists(out_vcf) and os.path.getsize(out_vcf) > 0:
        return f"{sample}: SKIP existing"
    cmd = (
        f'{JAVA} -jar {SNPEFF} ann '
        f'-dataDir {DATA_DIR} '
        f'Ref2337 {in_vcf} '
        f'> {out_vcf} 2> {log_file}'
    )
    ret = os.system(cmd)
    size = os.path.getsize(out_vcf) if os.path.exists(out_vcf) else 0
    return f"{sample}: exit={ret}, output_size={size}"

def main():
    with open(LIST_FILE, "r", encoding="utf-8", errors="replace") as f:
        samples = [x.strip() for x in f if x.strip()]
    print(f"Samples: {len(samples)}")
    print(f"Processes: {THREADS}")
    with multiprocessing.Pool(processes=THREADS) as pool:
        for result in pool.imap_unordered(run_one, samples):
            print(result, flush=True)

if __name__ == "__main__":
    main()
