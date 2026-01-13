# ClinDnaVar_NF（ClinSeq DNAseq Variant Calling Pipeline）

ClinDnaVar_NF is a clinical-grade DNA sequencing variant calling workflow for WES, targeted panels, and WGS, built with Nextflow DSL2. It is also suitable for high-depth panel analysis of clonal hematopoiesis of indeterminate potential (CHIP) mutations. The pipeline covers end-to-end processing from raw FASTQ to annotated variant outputs: QC, read filtering, alignment, duplicate handling, BQSR, germline calling, multiple somatic callers, annotation, and MAF aggregation. Caller execution is configurable via parameters so you can tailor the workflow to the assay and use case.

## Features

- FastQC + fastp for QC and read filtering
- BWA MEM alignment + samtools sort/index
- GATK MarkDuplicates + BQSR
- GATK HaplotypeCaller (GVCF) + GenotypeGVCFs + VQSR post-processing
- Somatic callers: Mutect2, LoFreq, VarDict, Pisces
- Annotation (dbSNP + ANNOVAR) and merged MAF output

## Requirements

This pipeline expects the following tools to be available in PATH (or in your container/module environment):

- nextflow
- fastqc
- fastp
- bwa
- samtools
- gatk
- bcftools
- lofreq
- vardict-java
- var2vcf_valid.pl
- pisces
- annovar (table_annovar.pl)
- perl

Optional (currently commented in workflow):

- gencore
- bedtools
- pindel
- pindel2vcf
- deepvariant

Custom scripts in `bin/` are used by the workflow and should remain in PATH (Nextflow automatically adds `bin/` to PATH):

- `get_germline.pl`
- `lofreq_reformat.pl`
- `annovar2maf_multitype.pl`
- `maf_sort_by_pos.pl`
- `combind_maf_v2.pl`

## Input

### Samplesheet

CSV with header and at least the following columns:

- `ID`: sample ID
- `R1`: path to read1 fastq.gz
- `R2`: path to read2 fastq.gz

Example:

```csv
ID,R1,R2
Sample1,/path/to/Sample1_R1.fastq.gz,/path/to/Sample1_R2.fastq.gz
```

Default input is set in `main.nf`:

- `params.input = "$projectDir/bin/samplesheet_20240706.csv"`

You can override it with `--input`.

## Running

Basic run:

```bash
nextflow run main.nf --input samplesheet.test.csv
```

Select callers with `--callers` (comma-separated):

```bash
nextflow run main.nf --input samplesheet.test.csv --callers mutect2,lofreq
nextflow run main.nf --input samplesheet.test.csv --callers germline
nextflow run main.nf --input samplesheet.test.csv --callers mutect2,vardict,pisces,germline
```

Supported values:

- `mutect2`
- `lofreq`
- `vardict`
- `pisces`
- `germline` (or `hapcaller`)

If no callers are selected, annotation and MAF combine steps are skipped.

Enable UMI-based duplicate handling (front/back 4bp UMI) with `--use_umi`:

```bash
nextflow run main.nf --input samplesheet.test.csv --use_umi true
```

## Parameters

Key parameters in `main.nf`:

- `--input`: samplesheet CSV
- `--outdir`: output directory (default: `results`)
- `--reference`: reference fasta
- `--ref_dict`: reference dict
- `--tmp`: temp dir
- `--anno`: GATK bundle path
- `--snpdb`: dbSNP VCF
- `--bed`: capture BED
- `--intervals`: intervals list
- `--callers`: enabled callers (default: `mutect2,lofreq,vardict,pisces,germline`)
- `--use_umi`: enable UMI path (default: `false`)

## Outputs

Main output folders (under `--outdir`):

- `report/`: FastQC, fastp, GATK metrics and BQSR reports, VQSR PDFs
- `bam/`: BQSR BAM/BAI
- `gvcf/`: GVCF files from HaplotypeCaller
- `vcf/`: normalized VCFs for each caller
- `maf/`: annotated MAFs and combined MAF

## Workflow Overview

1. **QC**: FastQC
2. **Filtering**: fastp
3. **Alignment**: BWA MEM + samtools sort/index
4. **Pre-processing**: GATK MarkDuplicates + BQSR
5. **Germline** (optional): HaplotypeCaller (GVCF) + GenotypeGVCFs + VQSR
6. **Somatic callers** (optional): Mutect2, LoFreq, VarDict, Pisces
7. **Annotation**: dbSNP + ANNOVAR
8. **MAF merge**: combine multi-caller MAFs

## Project Layout

```
.
├── assets/               # reference bed/intervals and dbSNP VCF
├── bin/                  # helper scripts used by pipeline
├── modules/              # Nextflow DSL2 modules
├── main.nf               # main workflow
├── samplesheet.test.csv  # example samplesheet
└── README.md
```

## Notes

- The pipeline uses fixed paths for some reference resources in `main.nf` and `bin/get_germline.pl`. Update these paths to match your environment.
- `get_germline.pl` generates a SLURM script and executes it. Ensure SLURM is available or replace with your scheduler/job runner.

## Module Files

- `modules/qc.nf`
- `modules/align.nf`
- `modules/gatk.nf`
- `modules/callers.nf`
- `modules/annotation.nf`
