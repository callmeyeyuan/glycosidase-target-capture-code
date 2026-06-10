# Code for targeted-capture sequencing analysis of human gut microbial glycoside hydrolase genes

This repository contains scripts used to generate processed data tables for a targeted-capture metagenomic sequencing dataset of human gut microbial glycoside hydrolase genes.

## Scripts

- `scripts/fastq_qc_summary.py`: summarizes FASTQ quality metrics.
- `scripts/sanitize_reference.py`: formats and sanitizes target and spike-in reference FASTA files.
- `scripts/run_target_alignment_qc.sh`: runs target-reference alignment and mapping QC.
- `scripts/summarize_target_alignment.py`: summarizes target-reference mapping results and generates sample-level and reference-level mapped-read tables.
- `scripts/generate_spikein_tables.py`: generates spike-in read-count and expected/detected spike-in check tables.
- `scripts/vcf_summary.py`: summarizes sample-level VCF files, including record counts, filter categories and reference-level variant counts.
- `scripts/generate_nonredundant_sites.py`: generates the cohort-level non-redundant variant table from sample-level VCF files.
- `scripts/run_snpeff_64.py`: runs SnpEff annotation for 64 sample-level VCF files.
- `scripts/summarize_snpeff_ann.py`: summarizes SnpEff annotation results.

## Input data

Raw FASTQ files are archived separately in GSA. Processed data files are archived separately in OMIX. This code release is intended to document the analysis steps used to generate quality-control, mapping, spike-in, variant and SnpEff summary tables.

## Notes

Sample metadata and spike-in addition design files were manually curated from de-identified sample information and spike-in addition records. SnpEff annotation was performed using a custom Ref2337 target-reference database.
