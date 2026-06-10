#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import gzip
from pathlib import Path
from collections import Counter, defaultdict


def open_text(path):
    path = str(path)
    if path.endswith('.gz'):
        return gzip.open(path, 'rt', encoding='utf-8', errors='replace')
    return open(path, 'r', encoding='utf-8', errors='replace')


def classify_variant(ref, alt):
    alt_list = [a for a in alt.split(',') if a and a != '.']
    if not alt_list:
        return 'no_alt'
    types = set()
    for a in alt_list:
        if len(ref) == 1 and len(a) == 1:
            types.add('SNP')
        elif len(ref) != len(a):
            types.add('INDEL')
        else:
            types.add('MNP_or_complex_substitution')
    return list(types)[0] if len(types) == 1 else 'mixed'


def parse_gt(sample_field, fmt_keys):
    if not sample_field or sample_field == '.':
        return 'missing'
    values = sample_field.split(':')
    fmt = dict(zip(fmt_keys, values))
    gt = fmt.get('GT', '')
    if not gt:
        return 'no_GT_field'
    if gt in ['.', './.', '.|.']:
        return 'missing'
    alleles = gt.replace('|', '/').split('/')
    if all(a == '0' for a in alleles):
        return 'hom_ref'
    if all(a != '0' and a != '.' for a in alleles):
        return 'hom_alt'
    if any(a != '0' and a != '.' for a in alleles):
        return 'het_or_mixed'
    return 'unknown'


def main():
    parser = argparse.ArgumentParser(description='Summarize VCF files.')
    parser.add_argument('vcf', help='Input .vcf or .vcf.gz')
    parser.add_argument('-o', '--out-dir', default='vcf_summary_out')
    parser.add_argument('--write-variants', action='store_true', help='Write variants_table.tsv')
    args = parser.parse_args()

    vcf_path = Path(args.vcf)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    total = 0
    multiallelic = 0
    chrom_counter = Counter()
    filter_counter = Counter()
    type_counter = Counter()
    qual_values = []
    sample_names = []
    sample_counter = defaultdict(Counter)

    variants_writer = None
    variants_handle = None
    if args.write_variants:
        variants_handle = open(out_dir / 'variants_table.tsv', 'w', encoding='utf-8', newline='')
        variants_writer = csv.writer(variants_handle, delimiter='\t')
        variants_writer.writerow(['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'variant_type', 'num_alt_alleles', 'INFO'])

    header_seen = False
    with open_text(vcf_path) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                header_seen = True
                parts = line.split('\t')
                if len(parts) > 9:
                    sample_names = parts[9:]
                continue
            if line.startswith('#'):
                continue

            parts = line.split('\t')
            if len(parts) < 8:
                continue

            chrom, pos, vid, ref, alt, qual, filt, info = parts[:8]
            fmt = parts[8] if len(parts) > 8 else ''
            sample_fields = parts[9:] if len(parts) > 9 else []

            total += 1
            chrom_counter[chrom] += 1
            filter_counter[filt] += 1
            if ',' in alt:
                multiallelic += 1

            vtype = classify_variant(ref, alt)
            type_counter[vtype] += 1

            if qual not in ['.', '']:
                try:
                    qual_values.append(float(qual))
                except ValueError:
                    pass

            if sample_names and sample_fields:
                fmt_keys = fmt.split(':') if fmt else []
                for sample, sample_field in zip(sample_names, sample_fields):
                    gt_class = parse_gt(sample_field, fmt_keys)
                    sample_counter[sample][gt_class] += 1

            if variants_writer:
                alts = [a for a in alt.split(',') if a and a != '.']
                variants_writer.writerow([chrom, pos, vid, ref, alt, qual, filt, vtype, len(alts), info])

    if variants_handle:
        variants_handle.close()

    if not header_seen:
        raise SystemExit('No #CHROM header found. Is this a valid VCF file?')

    with open(out_dir / 'vcf_overview.tsv', 'w', encoding='utf-8', newline='') as out:
        w = csv.writer(out, delimiter='\t')
        w.writerow(['metric', 'value'])
        w.writerow(['vcf_file', str(vcf_path)])
        w.writerow(['total_variant_records', total])
        w.writerow(['total_SNP_records', type_counter.get('SNP', 0)])
        w.writerow(['total_INDEL_records', type_counter.get('INDEL', 0)])
        w.writerow(['total_multiallelic_records', multiallelic])
        w.writerow(['num_reference_sequences_with_variants', len(chrom_counter)])
        w.writerow(['num_samples_in_vcf', len(sample_names)])
        w.writerow(['mean_QUAL', sum(qual_values) / len(qual_values) if qual_values else ''])
        w.writerow(['min_QUAL', min(qual_values) if qual_values else ''])
        w.writerow(['max_QUAL', max(qual_values) if qual_values else ''])

    with open(out_dir / 'variants_by_chrom.tsv', 'w', encoding='utf-8', newline='') as out:
        w = csv.writer(out, delimiter='\t')
        w.writerow(['CHROM', 'variant_records'])
        for chrom, count in chrom_counter.most_common():
            w.writerow([chrom, count])

    with open(out_dir / 'filter_summary.tsv', 'w', encoding='utf-8', newline='') as out:
        w = csv.writer(out, delimiter='\t')
        w.writerow(['FILTER', 'variant_records'])
        for filt, count in filter_counter.most_common():
            w.writerow([filt, count])

    with open(out_dir / 'variant_type_summary.tsv', 'w', encoding='utf-8', newline='') as out:
        w = csv.writer(out, delimiter='\t')
        w.writerow(['variant_type', 'variant_records'])
        for vtype, count in type_counter.most_common():
            w.writerow([vtype, count])

    with open(out_dir / 'sample_genotype_summary.tsv', 'w', encoding='utf-8', newline='') as out:
        w = csv.writer(out, delimiter='\t')
        w.writerow(['sample', 'hom_ref', 'het_or_mixed', 'hom_alt', 'missing', 'no_GT_field', 'unknown', 'total_records_seen'])
        for sample in sample_names:
            c = sample_counter[sample]
            w.writerow([sample, c.get('hom_ref', 0), c.get('het_or_mixed', 0), c.get('hom_alt', 0), c.get('missing', 0), c.get('no_GT_field', 0), c.get('unknown', 0), sum(c.values())])

    print('Done.')
    print('VCF:', vcf_path)
    print('Total variant records:', total)
    print('SNP records:', type_counter.get('SNP', 0))
    print('INDEL records:', type_counter.get('INDEL', 0))
    print('Reference sequences with variants:', len(chrom_counter))
    print('Samples in VCF:', len(sample_names))
    print('Output directory:', out_dir)


if __name__ == '__main__':
    main()
