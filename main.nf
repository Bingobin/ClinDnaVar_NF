#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input = "$projectDir/bin/samplesheet_20240706.csv"
params.outdir = "results"
params.reference = "/lustre/home/acct-medkkw/medlyb/database/annotation/gatk_ann/hg38/bwaindex2/Homo_sapiens_assembly38.fasta"
params.ref_dict = "/lustre/home/acct-medkkw/medlyb/database/annotation/gatk_ann/hg38/bwaindex2/Homo_sapiens_assembly38.dict"
params.tmp = "/lustre/home/acct-medkkw/medlyb/tmp"
params.anno = "/lustre/home/acct-medkkw/medlyb/database/annotation/gatk_ann/hg38/hg38bundle"
params.snpdb = "$projectDir/assets/wes.target.hg38.f.dbsnp_155.vcf.gz"
params.bed = "$projectDir/assets/wes.target.hg38.f.bed"
params.intervals = "$projectDir/assets/wes.target.hg38.f.intervals"
params.callers = "mutect2,lofreq,vardict,pisces,germline"
params.use_umi = false

log.info """\
    C H I P _ P I P E - N F    P I P E L I N E
    ==========================================
    Sample Info :  ${params.input}
    Project Dir :  $projectDir
    Reference   :  ${params.reference}
    Capture BED :  ${params.bed}
    Out Dir     :  ${params.outdir}
    """
    .stripIndent(true)

include { FASTQC; FastpFilter } from './modules/qc'
include { ALING_REF; RM_UMI_DUP } from './modules/align'
include { GATK_rmdup; GATK_BQSR; GENOTYPE; GATK_CALL_GERM } from './modules/gatk'
include {
    DeepVariant_CALL
    Mutect2_Call
    Pindel_Split_BED
    Pindel_Call
    Pindel_Merge
    LoFreq_Call
    VarDict_Call
    Pisces_Call
} from './modules/callers'
include { MUT_ANNOTATE; MAF_COMBIND } from './modules/annotation'
include { BAM_DEPTH; VCF_STATS } from './modules/stats'

workflow {
    Channel.fromPath(params.input)
        .splitCsv(header: true)
        .map{["${it.ID}" ,["${it.R1}", "${it.R2}"]]}
        .set {ch_rawfastq}
//    ch_rawfastq.view()
    FASTQC(ch_rawfastq)
    ch_cleanfastq = FastpFilter{ch_rawfastq}
    ch_aligned_bam = ALING_REF(ch_cleanfastq.fastq)
    if (params.use_umi) {
        ch_umi_bam = RM_UMI_DUP(ch_aligned_bam.bam)
        ch_bqsr = GATK_BQSR(ch_umi_bam.bam)
    } else {
        ch_rmdup = GATK_rmdup(ch_aligned_bam.bam)
        ch_bqsr = GATK_BQSR(ch_rmdup.bam)
    }
    BAM_DEPTH(ch_bqsr.bam)

    def caller_list = params.callers.split(',').collect { it.trim().toLowerCase() }.findAll { it }
    def anno_inputs = []

    if (caller_list.contains('germline') || caller_list.contains('hapcaller')) {
        ch_gvcf = GENOTYPE(ch_bqsr.bam)
        ch_germline = GATK_CALL_GERM(ch_gvcf)
        anno_inputs << ch_germline.vcf
        VCF_STATS(ch_germline.vcf)
    }

    if (caller_list.contains('mutect2')) {
        ch_mutect2 = Mutect2_Call(ch_bqsr.bam)
        anno_inputs << ch_mutect2
        VCF_STATS(ch_mutect2)
    }
//    ch_pindel_s = Pindel_Split_BED(ch_bqsr.bam, params.bed).transpose() | Pindel_Call
//    ch_pindel = Pindel_Merge(ch_pindel_s.groupTuple())

    if (caller_list.contains('lofreq')) {
        ch_lofreq = LoFreq_Call(ch_bqsr.bam)
        anno_inputs << ch_lofreq
        VCF_STATS(ch_lofreq)
    }

    if (caller_list.contains('vardict')) {
        ch_vardict = VarDict_Call(ch_bqsr.bam)
        anno_inputs << ch_vardict
        VCF_STATS(ch_vardict)
    }

    if (caller_list.contains('pisces')) {
        ch_pisces = Pisces_Call(ch_bqsr.bam)
        anno_inputs << ch_pisces
        VCF_STATS(ch_pisces)
    }
//    ch_anno = MUT_ANNOTATE(ch_mutect2.mix(ch_lofreq, ch_vardict, ch_pisces, ch_pindel))
    if (anno_inputs) {
        ch_anno_in = anno_inputs.size() == 1 ? anno_inputs[0] : Channel.mix(anno_inputs as Channel[])
        ch_anno = MUT_ANNOTATE(ch_anno_in)
        MAF_COMBIND(ch_anno.maf.groupTuple())
    } else {
        log.warn "No callers selected in params.callers; skipping MUT_ANNOTATE and MAF_COMBIND"
    }
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! See results --> $params.outdir\n" : "Oops.. someting went wrong" )
}
