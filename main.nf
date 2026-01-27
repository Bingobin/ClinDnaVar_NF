#!/usr/bin/env nextflow

nextflow.enable.dsl=2

log.info """\
    C l i n D n a V a r _ N F   P I P E L I N E
    ============================================
    Sample Info :  ${params.input}
    Project Dir :  $projectDir
    Reference   :  ${params.reference}
    Capture BED :  ${params.bed}
    Input Type  :  ${params.input_type}
    CNVKit Run  :  ${params.cnvkit}
    Only CNV    :  ${params.only_cnv}
    Only Depth  :  ${params.only_depth}
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
include { CNVKIT_BATCH } from './modules/cnv'

workflow {
    def input_header = new File(params.input.toString()).withReader { it.readLine() }
    def input_sep = (input_header != null && input_header.contains('\t')) ? '\t' : ','
    def input_type = params.input_type ? params.input_type.toString().toLowerCase() : 'fastq'
    if (input_type == 'bam') {
        Channel.fromPath(params.input)
            .splitCsv(header: true, sep: input_sep)
            .map { tuple(it.ID.toString(), file(it.BAM), file(it.BAI)) }
            .set { ch_input_bam }
    } else {
        Channel.fromPath(params.input)
            .splitCsv(header: true, sep: input_sep)
            .map{["${it.ID}" ,["${it.R1}", "${it.R2}"]]}
            .set {ch_rawfastq}
    //    ch_rawfastq.view()
        FASTQC(ch_rawfastq)
        ch_cleanfastq = FastpFilter{ch_rawfastq}
        ch_input_bam = ALING_REF(ch_cleanfastq.fastq).bam
    }
    if (input_type == 'bam') {
        ch_bqsr = ch_input_bam
    } else if (params.use_umi) {
        ch_umi_bam = RM_UMI_DUP(ch_input_bam)
        ch_bqsr = GATK_BQSR(ch_umi_bam.bam)
    } else {
        ch_rmdup = GATK_rmdup(ch_input_bam)
        ch_bqsr = GATK_BQSR(ch_rmdup.bam)
    }
    if (params.only_cnv) {
        CNVKIT_BATCH(ch_bqsr.bam)
        return
    }
    if (params.only_depth) {
        BAM_DEPTH(ch_bqsr.bam)
        return
    }
    BAM_DEPTH(ch_bqsr.bam)
    if (params.cnvkit) {
        CNVKIT_BATCH(ch_bqsr.bam)
    }

    def caller_list = params.callers.split(',').collect { it.trim().toLowerCase() }.findAll { it }
    def anno_inputs = []
    def vcf_stats_inputs = []

    if (caller_list.contains('germline') || caller_list.contains('hapcaller')) {
        ch_gvcf = GENOTYPE(ch_bqsr.bam)
        ch_germline = GATK_CALL_GERM(ch_gvcf)
        anno_inputs << ch_germline.vcf
        vcf_stats_inputs << ch_germline.vcf
    }

    if (caller_list.contains('mutect2')) {
        ch_mutect2 = Mutect2_Call(ch_bqsr.bam)
        anno_inputs << ch_mutect2
        vcf_stats_inputs << ch_mutect2
    }
//    ch_pindel_s = Pindel_Split_BED(ch_bqsr.bam, params.bed).transpose() | Pindel_Call
//    ch_pindel = Pindel_Merge(ch_pindel_s.groupTuple())

    if (caller_list.contains('lofreq')) {
        ch_lofreq = LoFreq_Call(ch_bqsr.bam)
        anno_inputs << ch_lofreq
        vcf_stats_inputs << ch_lofreq
    }

    if (caller_list.contains('vardict')) {
        ch_vardict = VarDict_Call(ch_bqsr.bam)
        anno_inputs << ch_vardict
        vcf_stats_inputs << ch_vardict
    }

    if (caller_list.contains('pisces')) {
        ch_pisces = Pisces_Call(ch_bqsr.bam)
        anno_inputs << ch_pisces
        vcf_stats_inputs << ch_pisces
    }
    if (vcf_stats_inputs) {
        ch_stats_in = vcf_stats_inputs.size() == 1 \
            ? vcf_stats_inputs[0] \
            : vcf_stats_inputs[0].mix(*vcf_stats_inputs[1..-1])
        VCF_STATS(ch_stats_in)
    }
//    ch_anno = MUT_ANNOTATE(ch_mutect2.mix(ch_lofreq, ch_vardict, ch_pisces, ch_pindel))
    if (anno_inputs) {
        ch_anno_in = anno_inputs.size() == 1 \
            ? anno_inputs[0] \
            : anno_inputs[0].mix(*anno_inputs[1..-1])
        ch_anno = MUT_ANNOTATE(ch_anno_in)
        MAF_COMBIND(ch_anno.maf.groupTuple())
    } else {
        log.warn "No callers selected in params.callers; skipping MUT_ANNOTATE and MAF_COMBIND"
    }
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! See results --> $params.outdir\n" : "Oops.. someting went wrong" )
}
