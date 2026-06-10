# Example commands

The commands below illustrate how the scripts in this code release can be used. Paths should be modified according to the local directory structure.

## 1. FASTQ quality summary

```bash
python scripts/fastq_qc_summary.py \
  --fastq_dir /path/to/raw_fastq \
  --output target_capture_qc/fastq_qc_summary.full.tsv

##2. Sanitize reference sequences

```bash
python scripts/sanitize_reference.py \
  --target_fasta reference_sequences/glycoside_hydrolase_targets_2337.fasta \
  --spikein_fasta reference_sequences/synthetic_spikeins_20.fasta \
  --output_dir reference_sequences/sanitized_reference

##3. Target-reference alignment and mapping QC

```bash
bash scripts/run_target_alignment_qc.sh

Then summarize alignment results:

```bash
python scripts/summarize_target_alignment.py \
  --input_dir /path/to/alignment_outputs \
  --output_sample target_capture_qc/sample_mapping_summary.tsv \
  --output_reference target_capture_qc/reference_mapped_reads_long.tsv

##4. Generate spike-in tables

```bash
python scripts/generate_spikein_tables.py \
  --mapped_reads target_capture_qc/reference_mapped_reads_long.tsv \
  --spikein_design metadata/spikein_addition_design.tsv \
  --output_matrix target_capture_qc/spikein_reads_matrix.tsv \
  --output_check target_capture_qc/spikein_expected_check.tsv


##5. Summarize VCF files

```bash
python scripts/vcf_summary.py \
  --vcf_dir merged_vcf \
  --output_dir variant_summary


##6. Generate non-redundant variant table

```bash
python scripts/generate_nonredundant_sites.py \
  --vcf_dir merged_vcf \
  --output nonredundant_variants/nonredundant_sites.tsv

##7. Run SnpEff annotation

```bash
python scripts/run_snpeff_64.py \
  --vcf_dir merged_vcf \
  --output_dir snpeff_annotation \
  --database Ref2337

##8. Summarize SnpEff annotations

```bash
python scripts/summarize_snpeff_ann.py \
  --ann_vcf_dir snpeff_annotation \
  --output_dir snpeff_annotation

Notes

These commands are examples. The exact argument names may need to be adjusted according to the script implementation. Sample metadata and spike-in addition design files were manually curated from de-identified sample records.

