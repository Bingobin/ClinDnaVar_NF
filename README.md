# ClinDnaVar_NF（ClinSeq DNAseq Variant Calling Pipeline）

ClinDnaVar_NF is a Nextflow-based DNA sequencing variant calling workflow for WES, targeted panels, and WGS. The workflow is designed for research and analysis of DNA variants, including variants with potential clinical relevance. It is also suitable for high-depth panel analysis of clonal hematopoiesis of indeterminate potential (CHIP) mutations. The pipeline covers end-to-end processing from raw FASTQ to annotated variant outputs: QC, read filtering, alignment, duplicate handling, BQSR, germline calling, multiple somatic callers, annotation, and MAF aggregation. Caller execution is configurable via parameters so you can tailor the workflow to the assay and use case.

## Features

- FastQC + fastp for QC and read filtering
- BWA MEM alignment + samtools sort/index
- GATK MarkDuplicates + BQSR
- GATK HaplotypeCaller (GVCF) + GenotypeGVCFs + VQSR post-processing
- Somatic callers: Mutect2, LoFreq, VarDict, Pisces
- Annotation (dbSNP + ANNOVAR) and merged MAF output
- Coverage and variant statistics (mosdepth_d4, bcftools stats)

## Requirements

### System requirements

ClinDnaVar_NF is a command-line Nextflow DSL2 workflow intended for Linux-based servers, high-performance computing clusters, or cloud analysis environments. The workflow was developed with Nextflow v26.04.1 build 12112 using DSL2 and standard Linux command-line tools.

No non-standard hardware is required for the workflow itself. For real WES, WGS, or high-depth panel datasets, use a Linux server or HPC environment with sufficient CPU cores, memory, and storage for BAM/CRAM, VCF, MAF, and annotation files. Default process resources are listed below and can be overridden in `nextflow.config`.

### Software dependencies

The workflow expects the following tools to be available in `PATH` or in the active container/module/Conda environment:

| Tool/resource | Tested version | Notes |
| --- | --- | --- |
| Nextflow | v26.04.1 build 12112 | DSL2 workflow engine |
| Java | Required by Nextflow, GATK, and VarDictJava | Use a Java version compatible with the installed Nextflow and GATK releases |
| FastQC | Required | FASTQ quality control |
| fastp | Required | Read filtering and UMI extraction in `panel_umi` mode |
| BWA | Required | Read alignment from FASTQ input |
| SAMtools | Required | BAM/CRAM processing and indexing |
| GATK BaseRecalibrator/ApplyBQSR | v4.3.0.0 | BQSR |
| GATK Mutect2 | v4.1.0.0 | Somatic variant calling in blood-only mode |
| BCFtools | v1.23 | VCF normalization, filtering, annotation, compression, and indexing |
| LoFreq | v2.1.5 | Somatic variant calling |
| VarDictJava | v1.8.3 | Somatic variant calling |
| Pisces | v5.2.10.49 | Somatic variant calling |
| ANNOVAR | Required | Functional annotation with `table_annovar.pl` |
| dbSNP | v155 | dbSNP annotation resource |
| COSMIC | v87 | Cancer mutation annotation resource |
| ClinVar | v20221231 | Clinical variant annotation resource |
| gnomAD | v3.1.2 | Population allele-frequency annotation resource |
| mosdepth_d4 | Required | Depth and coverage statistics |
| CNVKit | Required when `--cnvkit true` | CNV analysis |
| Perl | Required | Helper scripts in `bin/` |

Version numbers above list software or resource versions only.

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

## Installation

Clone the repository and enter the workflow directory:

```bash
git clone https://github.com/Bingobin/ClinDnaVar_NF.git
cd ClinDnaVar_NF
```

Install or load the required third-party tools listed above using your local module system, Conda environments, containers, or manually installed software. Nextflow automatically adds the repository `bin/` directory to `PATH` during workflow execution, so the bundled helper scripts do not require a separate installation step.

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

- `params.input = "$projectDir/assets/samplesheet.test.csv"`

This file is a samplesheet format template only; the repository does not include sequencing data. Replace its placeholder paths with paths to your own FASTQ files, or provide another samplesheet with `--input`.

## Running

Basic run:

```bash
nextflow run main.nf --input assets/samplesheet.test.csv
```

Select assay type with `--assay_mode`:

```bash
nextflow run main.nf --input assets/samplesheet.test.csv --assay_mode wes
nextflow run main.nf --input assets/samplesheet.test.csv --assay_mode wgs --bed '' --intervals ''
nextflow run main.nf --input assets/samplesheet.test.csv --assay_mode wgs --bed primary_contigs.bed --intervals primary_contigs.interval_list
nextflow run main.nf --input assets/samplesheet.test.csv --assay_mode wgs --bed '' --intervals '' --cnvkit true --cnv_access /path/to/access-5kb-mappable.hg38.bed
nextflow run main.nf --input assets/samplesheet.test.csv --assay_mode panel_no_umi
nextflow run main.nf --input assets/samplesheet.test.csv --assay_mode panel_umi
```

Select callers with `--callers` (comma-separated):

```bash
nextflow run main.nf --input assets/samplesheet.test.csv --callers mutect2,lofreq
nextflow run main.nf --input assets/samplesheet.test.csv --callers germline
nextflow run main.nf --input assets/samplesheet.test.csv --callers mutect2,vardict,pisces,germline
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
nextflow run main.nf --input assets/samplesheet.test.csv --use_umi true
nextflow run main.nf --input assets/samplesheet.test.csv --no_umi_panel_call true
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
- UMI mode assumes a dual 4 bp UMI structure, with one 4 bp UMI segment on each read (`fastp --umi_loc per_read --umi_len 4`)
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
│   └── samplesheet.test.csv
├── bin/                  # helper scripts used by pipeline
├── modules/              # Nextflow DSL2 modules
├── main.nf               # main workflow
├── nextflow.config       # default params and resources
├── LICENSE               # GPL-3.0 license
└── README.md
```

## Related Publication

- Please cite this article as: Liu Yb, Xu Yy, Yang Sz, Song H, Jiao B, Tan Y, Strategies and challenges in the detection of clonal hematopoiesis: current advances and future perspectives, *LabMed Discovery*, https://doi.org/10.1016/j.lmd.2026.100120

<p align="center">
  <img src="assets/CH.png" alt="Clonal hematopoiesis review figure" width="700" />
</p>

## License

ClinDnaVar_NF is released under the GNU General Public License v3.0 (GPL-3.0). See `LICENSE` for the full license text.

## Implementation Details

The workflow logic includes optional region-restricted input processing, BQSR, blood-only multi-caller somatic variant detection, VCF normalization, annotation, MAF conversion, multi-caller merge logic, and CHIP-specific filtering/curation. The executable workflow steps are summarized in the `Workflow Overview` section above and implemented in `main.nf`, `modules/`, and `bin/`.

## Notes

- The pipeline uses placeholder paths for reference resources in `nextflow.config` and helper scripts. Update these paths to match your environment before running on real data.
- In `panel_umi` mode, or when `--use_umi true` is used, fastp extracts UMIs as two 4 bp segments by default, one from each read, and writes them to read names with the `UMI` prefix for downstream UMI duplicate removal.
- When `--no_umi_panel_call true` is set, `MarkDuplicates` keeps duplicate reads in the BAM and only marks them, which is more suitable for non-UMI targeted panel somatic calling.
- Default mosdepth thresholds now follow `assay_mode`: `wes=10,20,50,100`, `wgs=10,20,30`, `panel_umi/panel_no_umi=100,200,500,1000`, unless `--depth_thresholds` is set explicitly.
- Pisces defaults follow `assay_mode`: `wes` uses `--minvf 0.02 --mindpfilter 20`, `wgs` uses `--minvf 0.05 --mindpfilter 10`, and both `panel_umi` and `panel_no_umi` use `--minvf 0.0005 --mindpfilter 500`.
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
