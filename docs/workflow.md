# Analysis workflow

This document summarizes the analysis workflow used to generate the processed data files for the targeted-capture metagenomic sequencing dataset of human gut microbial glycoside hydrolase genes.

## Step 1. FASTQ quality summary

Input:
- Paired-end FASTQ files deposited in GSA.

Script:
- `scripts/fastq_qc_summary.py`

Output:
- `target_capture_qc/fastq_qc_summary.full.tsv`

Description:
This step summarizes basic FASTQ-level quality metrics, including read counts, base counts and sequence quality statistics.

## Step 2. Reference preparation

Input:
- `reference_sequences/glycoside_hydrolase_targets_2337.fasta`
- `reference_sequences/synthetic_spikeins_20.fasta`

Script:
- `scripts/sanitize_reference.py`

Output:
- Sanitized target/spike-in reference files used for downstream alignment and annotation.

Description:
This step formats and sanitizes the target reference sequences and synthetic spike-in sequences before alignment and downstream analyses.

## Step 3. Target-reference alignment and mapping QC

Input:
- Paired-end FASTQ files.
- Combined target reference set containing 2,337 glycoside hydrolase target sequences and 20 synthetic spike-in sequences.

Scripts:
- `scripts/run_target_alignment_qc.sh`
- `scripts/summarize_target_alignment.py`

Outputs:
- `target_capture_qc/sample_mapping_summary.tsv`
- `target_capture_qc/reference_mapped_reads_long.tsv`

Description:
Reads were aligned to the combined target reference set. Mapping statistics were summarized at both sample level and reference-sequence level.

## Step 4. Spike-in detection summary

Input:
- `target_capture_qc/reference_mapped_reads_long.tsv`
- `metadata/spikein_addition_design.tsv`

Script:
- `scripts/generate_spikein_tables.py`

Outputs:
- `target_capture_qc/spikein_reads_matrix.tsv`
- `target_capture_qc/spikein_expected_check.tsv`

Description:
This step summarizes read counts mapped to synthetic spike-in sequences and compares expected versus detected spike-in sequences using the final spike-in identifier order.

## Step 5. VCF summarization

Input:
- Sample-level merged VCF files in `merged_vcf/`.

Script:
- `scripts/vcf_summary.py`

Outputs:
- `variant_summary/merge_counts.tsv`
- `variant_summary/merge_filter_counts.tsv`
- `variant_summary/merge_ref_counts.tsv`

Description:
This step summarizes merged VCF records at sample level, filter-category level and reference-sequence level.

## Step 6. Non-redundant variant generation

Input:
- Sample-level merged VCF files in `merged_vcf/`.

Script:
- `scripts/generate_nonredundant_sites.py`

Output:
- `nonredundant_variants/nonredundant_sites.tsv`

Description:
This step combines sample-level variant calls and generates a cohort-level non-redundant variant table.

## Step 7. SnpEff annotation

Input:
- Sample-level merged VCF files.
- Custom Ref2337 SnpEff target-reference database.

Script:
- `scripts/run_snpeff_64.py`

Output:
- `snpeff_annotation/*.ann.vcf.gz`

Description:
This step performs SnpEff annotation for each sample-level VCF file using a custom target-reference database.

## Step 8. SnpEff annotation summary

Input:
- SnpEff-annotated VCF files in `snpeff_annotation/`.

Script:
- `scripts/summarize_snpeff_ann.py`

Outputs:
- `snpeff_annotation/snpeff_sample_overview.tsv`
- `snpeff_annotation/snpeff_primary_category_total.tsv`
- `snpeff_annotation/snpeff_primary_effect_total.tsv`
- `snpeff_annotation/snpeff_primary_impact_total.tsv`
- `snpeff_annotation/snpeff_all_effect_total.tsv`
- `snpeff_annotation/snpeff_gene_total.tsv`
- `snpeff_annotation/snpeff_warning_total.tsv`

Description:
This step summarizes SnpEff annotation results across samples, including primary effects, predicted impacts, gene-level summaries and warning categories.
