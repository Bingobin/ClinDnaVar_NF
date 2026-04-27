# ClinDnaVar_NF（ClinSeq DNAseq Variant Calling Pipeline）

ClinDnaVar_NF is a clinical-grade DNA sequencing variant calling workflow for WES, targeted panels, and WGS, built with Nextflow DSL2. It is also suitable for high-depth panel analysis of clonal hematopoiesis of indeterminate potential (CHIP) mutations. The pipeline covers end-to-end processing from raw FASTQ to annotated variant outputs: QC, read filtering, alignment, duplicate handling, BQSR, germline calling, multiple somatic callers, annotation, and MAF aggregation. Caller execution is configurable via parameters so you can tailor the workflow to the assay and use case.

## Related Publication

- Please cite this article as: Liu Yb, Xu Yy, Yang Sz, Song H, Jiao B, Tan Y, Strategies and challenges in the detection of clonal hematopoiesis: current advances and future perspectives, *LabMed Discovery*, https://doi.org/10.1016/j.lmd.2026.100120
  
<p align="center">
  <img src="assets/CH.png" alt="Clonal hematopoiesis review figure" width="700" />
</p>

## Features

- FastQC + fastp for QC and read filtering
- BWA MEM alignment + samtools sort/index
- GATK MarkDuplicates + BQSR
- GATK HaplotypeCaller (GVCF) + GenotypeGVCFs + VQSR post-processing
- Somatic callers: Mutect2, LoFreq, VarDict, Pisces
- Annotation (dbSNP + ANNOVAR) and merged MAF output
- Coverage and variant statistics (mosdepth_d4, bcftools stats)

## Requirements

This pipeline expects the following tools to be available in PATH (or in your container/module environment):

- nextflow
- fastqc
- fastp
- bwa
- samtools
- gatk
- bcftools
- mosdepth_d4
- lofreq
- vardict-java
- var2vcf_valid.pl
- pisces
- annovar (table_annovar.pl)
- perl
- cnvkit.py

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
- `annotate_maf_with_vcf_info.py` (append raw INFO from a reference/normal VCF onto a MAF by variant match, including ANNOVAR/MAF-style indels represented with `-` and shifted coordinates)

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

For BAM input (`--input_type bam`, treated as BQSR-ready), use:

- `ID`: sample ID
- `BAM`: path to aligned BAM
- `BAI`: path to BAM index

Example:

```csv
ID,BAM,BAI
Sample1,/path/to/Sample1.bam,/path/to/Sample1.bam.bai
```

Default input is set in `nextflow.config`:

- `params.input = "$projectDir/bin/samplesheet_20240706.csv"`

You can override it with `--input`.

## Running

Basic run:

```bash
nextflow run main.nf --input samplesheet.test.csv
```

Select assay type with `--assay_mode`:

```bash
nextflow run main.nf --input samplesheet.test.csv --assay_mode wes
nextflow run main.nf --input samplesheet.test.csv --assay_mode wgs --bed '' --intervals ''
nextflow run main.nf --input samplesheet.test.csv --assay_mode wgs --bed primary_contigs.bed --intervals primary_contigs.interval_list
nextflow run main.nf --input samplesheet.test.csv --assay_mode wgs --bed '' --intervals '' --cnvkit true --cnv_access /path/to/access-5kb-mappable.hg38.bed
nextflow run main.nf --input samplesheet.test.csv --assay_mode panel_no_umi
nextflow run main.nf --input samplesheet.test.csv --assay_mode panel_umi
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

Legacy switches are still supported:

```bash
nextflow run main.nf --input samplesheet.test.csv --use_umi true
nextflow run main.nf --input samplesheet.test.csv --no_umi_panel_call true
```

## Parameters

Key parameters in `nextflow.config`:

- `--input`: samplesheet CSV
- `--outdir`: output directory (default: `results`)
- `--reference`: reference fasta
- `--ref_dict`: reference dict
- `--tmp`: temp dir
- `--anno`: GATK bundle path
- `--snpdb`: dbSNP VCF
- `--bed`: capture BED
- `--intervals`: intervals list
- `--assay_mode`: `wes`, `wgs`, `panel_umi`, or `panel_no_umi` (default: `wes`)
- `--callers`: enabled callers (default: `mutect2,lofreq,vardict,pisces,germline`)
- `--use_umi`: legacy UMI switch; `assay_mode=panel_umi` is preferred
- `--no_umi_panel_call`: legacy non-UMI panel switch; `assay_mode=panel_no_umi` is preferred
- `--gatk_markdup_max_records`: Picard/GATK MarkDuplicates `--MAX_RECORDS_IN_RAM`; default is automatic (`1000000` for WGS, `task.memory.toGiga() * 400000` otherwise)
- `--depth_thresholds`: optional mosdepth thresholds, for example `100,200,500,1000`
- `--input_type`: `fastq` or `bam` (default: `fastq`)
- `--maf_retain_mode`: `exonic` or `all` for multicaller MAF merge output (default: `exonic`)
- `--cnvkit`: enable CNVKit batch run (default: `false`)
- `--cnv_method`: CNVKit method override; default is automatic (`wgs` for WGS, `hybrid` otherwise)
- `--cnv_targets`: target BED for CNVKit; defaults to `--bed` for non-WGS modes when unset
- `--cnv_access`: accessible regions BED for CNVKit, useful for WGS
- `--cnv_annotate`: refFlat annotation for CNVKit
- `--cnv_normal`: optional normal BAM(s) for CNVKit; use `true` or `self` to pass `--normal` with no file (self-normal)
- `--cnvkit_conda`: conda env name/path for CNVKit (default: `cnvkit`)
- `--only_cnv`: run only CNVKit after BQSR (default: `false`)
- `--only_depth`: run only depth after BQSR (default: `false`)

CNVKit runs with `--processes $task.cpus` inside the CNV module. In WGS mode, CNVKit is run with `--method wgs`, target BED is not passed, and `--cnv_access` is used when provided. If `--cnv_access` is unset, WGS mode uses `--bed` as the CNVKit `--access` regions file.

Default process resources (can be overridden per process):

- `process.cpus = 4`
- `process.memory = 4.GB * task.cpus`; GATK MarkDuplicates uses 80% of this value as Java heap and leaves the rest for JVM/native overhead. In WGS mode, VarDict uses `assets/hg38_5k_150bpOL_seg.txt.gz` as its genome tiling BED, `-f 0.01`, and 95% of the task memory as Java heap.

## Outputs

Main output folders (under `--outdir`):

- `report/`: FastQC, fastp, GATK metrics and BQSR reports, VQSR PDFs
- `bam/`: BQSR BAM/BAI
- `gvcf/`: GVCF files from HaplotypeCaller
- `vcf/`: normalized VCFs for each caller and per-caller VCF stats
- `maf/`: annotated MAFs and combined MAF
- `depth/`: mosdepth_d4 coverage outputs
- `cnvkit/`: CNVKit batch outputs (CNN reference + per-sample output dir, including `.cnr/.cns/.call.cns` and `*.cnv.{scatter,diagram}.pdf`)

## Workflow Overview

1. **QC**: FastQC
2. **Filtering**: fastp
3. **Alignment**: BWA MEM + samtools sort/index
4. **Pre-processing**: GATK MarkDuplicates + BQSR
5. **Germline** (optional): HaplotypeCaller (GVCF) + GenotypeGVCFs + VQSR
6. **Somatic callers** (optional): Mutect2, LoFreq, VarDict, Pisces
7. **Annotation**: dbSNP + ANNOVAR
8. **MAF merge**: combine multi-caller MAFs
9. **CNV** (optional): CNVKit batch on BQSR-ready BAM (`--cnvkit true`), plus explicit `cnvkit.py scatter/diagram` outputs.

## Project Layout

```
.
├── assets/               # reference bed/intervals and dbSNP VCF
├── bin/                  # helper scripts used by pipeline
├── modules/              # Nextflow DSL2 modules
├── main.nf               # main workflow
├── nextflow.config       # default params and resources
├── samplesheet.test.csv  # example samplesheet
└── README.md
```

## Notes

- The pipeline uses fixed paths for some reference resources in `nextflow.config` and `bin/get_germline.pl`. Update these paths to match your environment.
- When `--no_umi_panel_call true` is set, `MarkDuplicates` keeps duplicate reads in the BAM and only marks them, which is more suitable for non-UMI targeted panel somatic calling.
- Default mosdepth thresholds now follow `assay_mode`: `wes=10,20,50,100`, `wgs=10,20,30`, `panel_umi/panel_no_umi=100,200,500,1000`, unless `--depth_thresholds` is set explicitly.
- `--maf_retain_mode exonic` keeps coding/splice/frame-shift style records in the merged MAF; `--maf_retain_mode all` keeps all annotated records.
- For `--assay_mode wgs`, pass empty `--bed '' --intervals ''` for unrestricted WGS, or pass whole-genome-appropriate BED/interval-list files to keep primary contigs and exclude decoys, HLA alt contigs, or other scattered contigs. VarDict always uses the bundled WGS tiling BED `assets/hg38_5k_150bpOL_seg.txt.gz` in WGS mode.
- Conda is enabled in `nextflow.config`; ensure `conda` is available in PATH (or run with `-with-conda`).
- `get_germline.pl` generates a SLURM script and executes it. Ensure SLURM is available or replace with your scheduler/job runner.

## Module Files

- `modules/qc.nf`
- `modules/align.nf`
- `modules/gatk.nf`
- `modules/callers.nf`
- `modules/annotation.nf`
- `modules/stats.nf`
