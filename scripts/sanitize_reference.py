#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
from pathlib import Path

VALID = set("ACGTN")

def read_fasta(path):
    header = None
    seq_parts = []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n\r")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts).upper()
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if header is not None:
            yield header, "".join(seq_parts).upper()

def clean_seq(seq):
    invalid = 0
    out = []
    for ch in seq.upper():
        if ch in VALID:
            out.append(ch)
        else:
            out.append("N")
            invalid += 1
    return "".join(out), invalid

def wrap(seq, width=80):
    for i in range(0, len(seq), width):
        yield seq[i:i+width]

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--genes", required=True)
    p.add_argument("--spikeins", required=True)
    p.add_argument("--out-fasta", required=True)
    p.add_argument("--out-map", required=True)
    args = p.parse_args()

    records = []
    total_invalid = 0

    for i, (header, raw_seq) in enumerate(read_fasta(args.genes), start=1):
        seq, invalid = clean_seq(raw_seq)
        total_invalid += invalid
        records.append((f"gene_{i:04d}", "gene", i, len(seq), invalid, header, seq))

    for i, (header, raw_seq) in enumerate(read_fasta(args.spikeins), start=1):
        seq, invalid = clean_seq(raw_seq)
        total_invalid += invalid
        records.append((f"spikein_{i:02d}", "spikein", i, len(seq), invalid, header, seq))

    Path(args.out_fasta).parent.mkdir(parents=True, exist_ok=True)

    with open(args.out_fasta, "w", encoding="utf-8") as out:
        for new_id, typ, idx, length, invalid, header, seq in records:
            out.write(f">{new_id}\n")
            for part in wrap(seq):
                out.write(part + "\n")

    with open(args.out_map, "w", encoding="utf-8", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["new_id", "record_type", "index", "length_bp", "invalid_chars_replaced_by_N", "original_header"])
        for new_id, typ, idx, length, invalid, header, seq in records:
            w.writerow([new_id, typ, idx, length, invalid, header])

    print(f"Genes: {sum(1 for r in records if r[1] == 'gene')}")
    print(f"Spike-ins: {sum(1 for r in records if r[1] == 'spikein')}")
    print(f"Total records: {len(records)}")
    print(f"Total non-ACGTN characters replaced by N: {total_invalid}")
    print(f"Wrote FASTA: {args.out_fasta}")
    print(f"Wrote map: {args.out_map}")

if __name__ == "__main__":
    main()
