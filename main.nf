#!/usr/bin/env nextflow

nextflow.enable.dsl=2

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
    log.info("""\
        C l i n D n a V a r _ N F   P I P E L I N E
        ============================================
        Sample Info :  ${params.input}
        Project Dir :  $projectDir
        Reference   :  ${params.reference}
        Assay Mode  :  ${params.assay_mode}
        Capture BED :  ${params.bed}
        Input Type  :  ${params.input_type}
        No UMI Call :  ${params.no_umi_panel_call}
        CNVKit Run  :  ${params.cnvkit}
        Only CNV    :  ${params.only_cnv}
        Only Depth  :  ${params.only_depth}
        Out Dir     :  ${params.outdir}
        """
        .stripIndent(true))

    def input_header = new File(params.input.toString()).withReader { it.readLine() }
    def input_sep = (input_header != null && input_header.contains('\t')) ? '\t' : ','
    def input_type = params.input_type ? params.input_type.toString().toLowerCase() : 'fastq'
    def assay_mode = params.assay_mode ? params.assay_mode.toString().toLowerCase() : 'wes'
    def valid_assay_modes = ['wes', 'wgs', 'panel_umi', 'panel_no_umi'] as Set
    if (!valid_assay_modes.contains(assay_mode)) {
        error "Unsupported --assay_mode '${params.assay_mode}'. Supported values: wes, wgs, panel_umi, panel_no_umi"
    }
    def use_umi_mode = assay_mode == 'panel_umi' || params.use_umi

    def ch_input_bam
    if (input_type == 'bam') {
        ch_input_bam = Channel.fromPath(params.input)
            .splitCsv(header: true, sep: input_sep)
            .map {["${it.ID}" ,"${it.BAM}", "${it.BAI}"]}
    } else {
        def ch_rawfastq = Channel.fromPath(params.input)
            .splitCsv(header: true, sep: input_sep)
            .map{["${it.ID}" ,["${it.R1}", "${it.R2}"]]}
    //    ch_rawfastq.view()
        FASTQC(ch_rawfastq)
        def ch_cleanfastq = FastpFilter(ch_rawfastq)
        ch_input_bam = ALING_REF(ch_cleanfastq.fastq).bam
    }
    def ch_bqsr
    if (input_type == 'bam') {
        ch_bqsr = [bam: ch_input_bam]
    } else if (use_umi_mode) {
        def ch_umi_bam = RM_UMI_DUP(ch_input_bam)
        ch_bqsr = GATK_BQSR(ch_umi_bam.bam)
    } else {
        def ch_rmdup = GATK_rmdup(ch_input_bam)
        ch_bqsr = GATK_BQSR(ch_rmdup.bam)
    }

    if (params.only_cnv) {
        CNVKIT_BATCH(ch_bqsr.bam)
    } else if (params.only_depth) {
        BAM_DEPTH(ch_bqsr.bam)
    } else {
        BAM_DEPTH(ch_bqsr.bam)
        if (params.cnvkit) {
            CNVKIT_BATCH(ch_bqsr.bam)
        }

        def caller_list = params.callers.split(',').collect { it.trim().toLowerCase() }.findAll { it }
        def anno_inputs = []
        def vcf_stats_inputs = []

        if (caller_list.contains('germline') || caller_list.contains('hapcaller')) {
            def ch_gvcf = GENOTYPE(ch_bqsr.bam)
            def ch_germline = GATK_CALL_GERM(ch_gvcf)
            anno_inputs << ch_germline.vcf
            vcf_stats_inputs << ch_germline.vcf
        }

        if (caller_list.contains('mutect2')) {
            def ch_mutect2 = Mutect2_Call(ch_bqsr.bam)
            anno_inputs << ch_mutect2
            vcf_stats_inputs << ch_mutect2
        }
//        ch_pindel_s = Pindel_Split_BED(ch_bqsr.bam, params.bed).transpose() | Pindel_Call
//        ch_pindel = Pindel_Merge(ch_pindel_s.groupTuple())

        if (caller_list.contains('lofreq')) {
            def ch_lofreq = LoFreq_Call(ch_bqsr.bam)
            anno_inputs << ch_lofreq
            vcf_stats_inputs << ch_lofreq
        }

        if (caller_list.contains('vardict')) {
            def ch_vardict = VarDict_Call(ch_bqsr.bam)
            anno_inputs << ch_vardict
            vcf_stats_inputs << ch_vardict
        }

        if (caller_list.contains('pisces')) {
            def ch_pisces = Pisces_Call(ch_bqsr.bam)
            anno_inputs << ch_pisces
            vcf_stats_inputs << ch_pisces
        }
        if (vcf_stats_inputs) {
            def ch_stats_in = vcf_stats_inputs.size() == 1 \
                ? vcf_stats_inputs[0] \
                : vcf_stats_inputs.inject(null) { acc, ch -> acc == null ? ch : acc.mix(ch) }
            VCF_STATS(ch_stats_in)
        }
//        ch_anno = MUT_ANNOTATE(ch_mutect2.mix(ch_lofreq, ch_vardict, ch_pisces, ch_pindel))
        if (anno_inputs) {
            def ch_anno_in = anno_inputs.size() == 1 \
                ? anno_inputs[0] \
                : anno_inputs.inject(null) { acc, ch -> acc == null ? ch : acc.mix(ch) }
            def ch_anno = MUT_ANNOTATE(ch_anno_in)
            MAF_COMBIND(ch_anno.maf.groupTuple())
        } else {
            log.warn "No callers selected in params.callers; skipping MUT_ANNOTATE and MAF_COMBIND"
        }
    }
}
