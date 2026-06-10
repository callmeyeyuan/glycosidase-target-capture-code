#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import gzip
import os
import re
from pathlib import Path
from multiprocessing import Pool, cpu_count


ILLUMINA_RE = re.compile(
    r"^@(?P<instrument>[^:]+):(?P<run_id>[^:]+):(?P<flowcell_id>[^:]+):"
    r"(?P<lane>[^:]+):(?P<tile>[^:]+):(?P<x_pos>[^:]+):(?P<y_pos>\S+)"
)


def open_fastq(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def infer_read_direction(filename):
    name = Path(filename).name
    if re.search(r"(_R1_|\.R1\.|_1\.fq|_1\.fastq|_1\.fq\.gz|_1\.fastq\.gz)", name):
        return "R1"
    if re.search(r"(_R2_|\.R2\.|_2\.fq|_2\.fastq|_2\.fq\.gz|_2\.fastq\.gz)", name):
        return "R2"
    return "unknown"


def infer_platform_from_header(instrument):
    """
    This is only a rough inference based on the instrument ID in FASTQ headers.
    It cannot replace the official sequencing platform information from the provider.
    """
    if not instrument:
        return "unknown"

    inst = instrument.upper()

    if inst.startswith("A"):
        return "Illumina-like, possibly NovaSeq"
    if inst.startswith("NS"):
        return "Illumina NextSeq-like"
    if inst.startswith("MN"):
        return "Illumina MiniSeq-like"
    if inst.startswith("M"):
        return "Illumina MiSeq-like"
    if inst.startswith("D"):
        return "Illumina HiSeq 2500-like"
    if inst.startswith("HWI") or inst.startswith("HISEQ"):
        return "Illumina HiSeq-like"
    if inst.startswith("K") or inst.startswith("ST"):
        return "Illumina-like"

    return "Illumina-like or unknown"


def summarize_fastq(args_tuple):
    path, max_reads = args_tuple
    path = Path(path)

    reads = 0
    bases = 0
    gc = 0
    q20_bases = 0
    q30_bases = 0
    n_bases = 0

    min_len = None
    max_len = 0
    first_header = ""
    instrument = ""
    run_id = ""
    flowcell_id = ""
    lane = ""

    try:
        with open_fastq(path) as fh:
            while True:
                header = fh.readline().rstrip()
                if not header:
                    break

                seq = fh.readline().rstrip().upper()
                plus = fh.readline()
                qual = fh.readline().rstrip()

                if not first_header:
                    first_header = header
                    m = ILLUMINA_RE.match(header)
                    if m:
                        instrument = m.group("instrument")
                        run_id = m.group("run_id")
                        flowcell_id = m.group("flowcell_id")
                        lane = m.group("lane")

                reads += 1
                length = len(seq)
                bases += length

                if min_len is None or length < min_len:
                    min_len = length
                if length > max_len:
                    max_len = length

                gc += seq.count("G") + seq.count("C")
                n_bases += seq.count("N")

                for ch in qual:
                    q = ord(ch) - 33
                    if q >= 20:
                        q20_bases += 1
                    if q >= 30:
                        q30_bases += 1

                if max_reads and reads >= max_reads:
                    break

        avg_len = bases / reads if reads else 0
        gc_pct = gc / bases * 100 if bases else 0
        n_pct = n_bases / bases * 100 if bases else 0
        q20_pct = q20_bases / bases * 100 if bases else 0
        q30_pct = q30_bases / bases * 100 if bases else 0

        return {
            "file": str(path),
            "file_name": path.name,
            "read_direction": infer_read_direction(path.name),
            "file_size_GiB": round(path.stat().st_size / 1024**3, 4),
            "reads_counted": reads,
            "bases_counted": bases,
            "min_read_len": min_len if min_len is not None else 0,
            "max_read_len": max_len,
            "avg_read_len": round(avg_len, 2),
            "GC_percent": round(gc_pct, 2),
            "N_percent": round(n_pct, 4),
            "Q20_percent": round(q20_pct, 2),
            "Q30_percent": round(q30_pct, 2),
            "first_header": first_header,
            "instrument": instrument,
            "run_id": run_id,
            "flowcell_id": flowcell_id,
            "lane": lane,
            "platform_inferred": infer_platform_from_header(instrument),
            "status": "OK",
            "error": "",
        }

    except Exception as e:
        return {
            "file": str(path),
            "file_name": path.name,
            "read_direction": infer_read_direction(path.name),
            "file_size_GiB": round(path.stat().st_size / 1024**3, 4) if path.exists() else 0,
            "reads_counted": 0,
            "bases_counted": 0,
            "min_read_len": 0,
            "max_read_len": 0,
            "avg_read_len": 0,
            "GC_percent": 0,
            "N_percent": 0,
            "Q20_percent": 0,
            "Q30_percent": 0,
            "first_header": "",
            "instrument": "",
            "run_id": "",
            "flowcell_id": "",
            "lane": "",
            "platform_inferred": "unknown",
            "status": "ERROR",
            "error": str(e),
        }


def find_fastqs(input_path):
    input_path = Path(input_path)
    if input_path.is_file():
        return [input_path]

    patterns = ["*.fastq.gz", "*.fq.gz", "*.fastq", "*.fq"]
    files = []
    for pattern in patterns:
        files.extend(input_path.rglob(pattern))
    return sorted(files)


def main():
    parser = argparse.ArgumentParser(
        description="Summarize FASTQ metrics: reads, bases, GC, Q20, Q30, read length and header information."
    )
    parser.add_argument("input", help="FASTQ file or directory containing FASTQ files")
    parser.add_argument("-o", "--output", default="fastq_qc_summary.tsv", help="Output TSV file")
    parser.add_argument("-t", "--threads", type=int, default=min(4, cpu_count()), help="Number of parallel processes")
    parser.add_argument(
        "--max-reads",
        type=int,
        default=0,
        help="Maximum reads to scan per FASTQ file. Use 0 for full file. For quick estimation, use 100000.",
    )

    args = parser.parse_args()

    fastqs = find_fastqs(args.input)
    if not fastqs:
        raise SystemExit(f"No FASTQ files found in {args.input}")

    print(f"Found {len(fastqs)} FASTQ files.")
    print(f"Output: {args.output}")
    if args.max_reads:
        print(f"Sampling mode: first {args.max_reads} reads per file.")
    else:
        print("Full scan mode: all reads will be counted.")

    worker_args = [(str(p), args.max_reads) for p in fastqs]

    with Pool(processes=args.threads) as pool:
        results = pool.map(summarize_fastq, worker_args)

    fields = [
        "file",
        "file_name",
        "read_direction",
        "file_size_GiB",
        "reads_counted",
        "bases_counted",
        "min_read_len",
        "max_read_len",
        "avg_read_len",
        "GC_percent",
        "N_percent",
        "Q20_percent",
        "Q30_percent",
        "first_header",
        "instrument",
        "run_id",
        "flowcell_id",
        "lane",
        "platform_inferred",
        "status",
        "error",
    ]

    with open(args.output, "w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in results:
            writer.writerow(row)

    total_files = len(results)
    ok_files = sum(1 for r in results if r["status"] == "OK")
    total_reads = sum(int(r["reads_counted"]) for r in results if r["status"] == "OK")
    total_bases = sum(int(r["bases_counted"]) for r in results if r["status"] == "OK")

    print("Done.")
    print(f"FASTQ files processed: {ok_files}/{total_files}")
    print(f"Total reads counted: {total_reads}")
    print(f"Total bases counted: {total_bases}")


if __name__ == "__main__":
    main()
