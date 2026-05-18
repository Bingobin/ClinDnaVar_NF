def paramValue(value) {
    if (value == null) {
        return ""
    }

    def text = value.toString().trim()
    if (!text || text.equalsIgnoreCase("true") || text.equalsIgnoreCase("false") || text.equalsIgnoreCase("null")) {
        return ""
    }

    return text
}

def gatkIntervalsArg(intervals) {
    def value = paramValue(intervals)
    return value ? "-L ${value}" : ""
}

def gatkJavaOptions(task, extra = "") {
    def heapGb = Math.max(1, Math.floor(task.memory.toGiga() * 0.8) as int)
    def options = "-Xmx${heapGb}g"
    return extra ? "${options} ${extra}" : options
}

process GATK_rmdup {

    tag { "MarkDuplicates on $sample_id" }
    publishDir "$params.outdir/report", pattern: "*.metrics", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.rmdup.bam"), path("${sample_id}.rmdup.bam.bai"), emit: bam
    tuple val(sample_id), path("${sample_id}.rmdup.metrics"), emit: report

    script:
    def assayMode = params.assay_mode ? params.assay_mode.toString().toLowerCase() : "wes"
    def removeDup = (params.no_umi_panel_call || assayMode == "panel_no_umi") ? "false" : "true"
    def memoryGb = Math.max(1, Math.floor(task.memory.toGiga()) as int)
    def defaultMaxRecords = assayMode == "wgs" ? 1000000 : memoryGb * 400000
    def maxRecords = paramValue(params.gatk_markdup_max_records) ? params.gatk_markdup_max_records.toString().toInteger() : defaultMaxRecords
    def javaOptions = gatkJavaOptions(task, "-XX:+UseParallelGC")
    """
    gatk --java-options "${javaOptions}" MarkDuplicates -I ${bam[0]} -O ${sample_id}.rmdup.bam -M ${sample_id}.rmdup.metrics --REMOVE_SEQUENCING_DUPLICATES ${removeDup} --ASSUME_SORT_ORDER coordinate --VALIDATION_STRINGENCY LENIENT --MAX_RECORDS_IN_RAM ${maxRecords} --TMP_DIR ${params.tmp}
    samtools index ${sample_id}.rmdup.bam
    """
}

process GATK_BQSR {
    tag { "BQSR on $sample_id" }

    publishDir "$params.outdir/bam", mode:'copy', pattern: "*.bqsr.bam"
    publishDir "$params.outdir/bam", mode:'copy', pattern: "*.bqsr.bai"
    publishDir "$params.outdir/report", mode:'copy', pattern: "*.recal_data.table"
    publishDir "$params.outdir/report", mode:'copy', pattern: "*.bqsr.bam.stats"

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path( "${sample_id}.bqsr.bam"), path("${sample_id}.bqsr.bai"), emit: 'bam'
    tuple val(sample_id), path( "${sample_id}.recal_data.table"), path( "${sample_id}.bqsr.bam.stats") , emit: 'report'

    script:
    def intervalsArg = gatkIntervalsArg(params.intervals)
    def javaOptions = gatkJavaOptions(task)
    def applyJavaOptions = gatkJavaOptions(task, "-Dsamjdk.compression_level=6")
    """
    gatk --java-options "${javaOptions}" BaseRecalibrator -I ${bam[0]} -R ${params.reference} ${intervalsArg} \
    --known-sites $params.anno/Mills_and_1000G_gold_standard.indels.hg38.vcf \
    --known-sites $params.anno/dbsnp_146.hg38.vcf \
    --known-sites $params.anno/Homo_sapiens_assembly38.known_indels.vcf \
    -O ${sample_id}.recal_data.table
    gatk --java-options "${applyJavaOptions}" ApplyBQSR -I ${bam[0]} -R ${params.reference} --bqsr-recal-file ${sample_id}.recal_data.table -O ${sample_id}.bqsr.bam
    samtools stats -@ $task.cpus ${sample_id}.bqsr.bam > ${sample_id}.bqsr.bam.stats
    """
}

process GENOTYPE {
    tag { "GenotypeGVCFs on $sample_id" }
    publishDir "$params.outdir/gvcf", mode:'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.g.vcf.*")

    script:
    def intervalsArg = gatkIntervalsArg(params.intervals)
    def javaOptions = gatkJavaOptions(task)
    """
    gatk --java-options "${javaOptions}" HaplotypeCaller -R ${params.reference} -I ${bam[0]}  -ERC GVCF -O ${sample_id}.g.vcf.gz ${intervalsArg} -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation
    """
}

process GATK_GENOTYPE_GVCF {
    tag { "GenotypeGVCFs on $sample_id" }

    input:
    tuple val(sample_id), path(gvcf)

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), path("${sample_id}.vcf.gz.tbi"), emit: vcf

    script:
    def memory = gatkJavaOptions(task)
    def intervalsArg = gatkIntervalsArg(params.intervals)
    def mergeIntervalsArg = intervalsArg ? "--merge-input-intervals" : ""
    def genotypeInput = intervalsArg ? "gendb://${sample_id}_DB" : "${gvcf[0]}"
    def genomicsDbImport = intervalsArg ? """
    gatk --java-options "${memory}" GenomicsDBImport \\
        --genomicsdb-workspace-path ${sample_id}_DB \\
        --tmp-dir ${params.tmp} \\
        ${intervalsArg} \\
        -V ${gvcf[0]} \\
        ${mergeIntervalsArg}
    """ : """
    echo "Skipping GenomicsDBImport for ${sample_id}; no GATK intervals were provided"
    """
    def cleanupGenomicsDb = intervalsArg ? "rm -rf ${sample_id}_DB" : ""
    """
    ${genomicsDbImport}

    gatk --java-options "${memory}" GenotypeGVCFs \\
        -R ${params.reference} \\
        -V ${genotypeInput} \\
        -O ${sample_id}.vcf.gz \\
        --tmp-dir ${params.tmp} \\
        --dbsnp ${params.snpdb}

    ${cleanupGenomicsDb}
    """
}

process GATK_SNP_VQSR {
    tag { "SNP VQSR on $sample_id" }
    publishDir "$params.outdir/report", mode:'copy', pattern: "*.pdf"

    input:
    tuple val(sample_id), path(vcf), path(tbi)

    output:
    tuple val(sample_id), path("${sample_id}.VQSR.snp.vcf.gz"), path("${sample_id}.VQSR.snp.vcf.gz.tbi"), emit: vcf
    tuple val(sample_id), path("*.pdf"), emit: report

    script:
    def memory = gatkJavaOptions(task)
    """
    gatk --java-options "${memory}" SelectVariants \\
        -R ${params.reference} \\
        -V ${vcf} \\
        --tmp-dir ${params.tmp} \\
        -O ${sample_id}.raw.snp.vcf.gz \\
        --exclude-non-variants \\
        --select-type-to-include SNP

    gatk --java-options "${memory}" VariantRecalibrator \\
        -R ${params.reference} \\
        -V ${sample_id}.raw.snp.vcf.gz \\
        --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${params.anno}/hapmap_3.3.hg38.vcf.gz \\
        --resource:omni,known=false,training=true,truth=false,prior=12.0 ${params.anno}/1000G_omni2.5.hg38.vcf.gz \\
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${params.anno}/1000G_phase1.snps.high_confidence.hg38.vcf \\
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${params.anno}/dbsnp_146.hg38.vcf \\
        -an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum \\
        -mode SNP \\
        --output ${sample_id}.recalibrate_SNP.recal \\
        --tranches-file ${sample_id}.recalibrate_SNP.tranches \\
        --rscript-file ${sample_id}.recalibrate_SNP_plots.R

    gatk --java-options "${memory}" ApplyVQSR \\
        -R ${params.reference} \\
        -V ${sample_id}.raw.snp.vcf.gz \\
        --tranches-file ${sample_id}.recalibrate_SNP.tranches \\
        --recal-file ${sample_id}.recalibrate_SNP.recal \\
        -O ${sample_id}.VQSR.snp.vcf.gz \\
        -mode SNP

    rm -f ${sample_id}.raw.snp.vcf.gz ${sample_id}.raw.snp.vcf.gz.tbi
    """
}

process GATK_INDEL_VQSR {
    tag { "INDEL VQSR on $sample_id" }
    publishDir "$params.outdir/report", mode:'copy', pattern: "*.pdf"

    input:
    tuple val(sample_id), path(vcf), path(tbi)

    output:
    tuple val(sample_id), path("${sample_id}.VQSR.indel.vcf.gz"), path("${sample_id}.VQSR.indel.vcf.gz.tbi"), emit: vcf
    tuple val(sample_id), path("*.pdf"), emit: report

    script:
    def memory = gatkJavaOptions(task)
    """
    gatk --java-options "${memory}" SelectVariants \\
        -R ${params.reference} \\
        -V ${vcf} \\
        --tmp-dir ${params.tmp} \\
        -O ${sample_id}.raw.indel.vcf.gz \\
        --exclude-non-variants \\
        --select-type-to-include INDEL

    gatk --java-options "${memory}" VariantRecalibrator \\
        -R ${params.reference} \\
        -V ${sample_id}.raw.indel.vcf.gz \\
        --resource:mills,known=true,training=true,truth=true,prior=12.0 ${params.anno}/Mills_and_1000G_gold_standard.indels.hg38.vcf \\
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${params.anno}/dbsnp_146.hg38.vcf \\
        -an SOR -an MQ -an MQRankSum -an ReadPosRankSum \\
        --max-gaussians 4 \\
        -mode INDEL \\
        -O ${sample_id}.recalibrate_INDEL.recal \\
        --tranches-file ${sample_id}.recalibrate_INDEL.tranches \\
        --rscript-file ${sample_id}.recalibrate_INDEL_plots.R

    gatk --java-options "${memory}" ApplyVQSR \\
        -R ${params.reference} \\
        -V ${sample_id}.raw.indel.vcf.gz \\
        --tranches-file ${sample_id}.recalibrate_INDEL.tranches \\
        --recal-file ${sample_id}.recalibrate_INDEL.recal \\
        -O ${sample_id}.VQSR.indel.vcf.gz \\
        -mode INDEL

    rm -f ${sample_id}.raw.indel.vcf.gz ${sample_id}.raw.indel.vcf.gz.tbi
    """
}

process GATK_MERGE_VQSR {
    tag { "Merge VQSR on $sample_id" }

    input:
    tuple val(sample_id), path(snp_vcf), path(snp_tbi), path(indel_vcf), path(indel_tbi)

    output:
    tuple val(sample_id), path("${sample_id}.VQSR.sort.vcf.gz"), path("${sample_id}.VQSR.sort.vcf.gz.tbi"), emit: vcf

    script:
    def memory = gatkJavaOptions(task)
    """
    gatk --java-options "${memory}" SortVcf -I ${snp_vcf} -O ${sample_id}.VQSR.snp.sort.vcf.gz
    gatk --java-options "${memory}" SortVcf -I ${indel_vcf} -O ${sample_id}.VQSR.indel.sort.vcf.gz
    gatk --java-options "${memory}" MergeVcfs -I ${sample_id}.VQSR.snp.sort.vcf.gz -I ${sample_id}.VQSR.indel.sort.vcf.gz -O ${sample_id}.VQSR.sort.vcf.gz
    gatk --java-options "${memory}" IndexFeatureFile -I ${sample_id}.VQSR.sort.vcf.gz

    rm -f ${sample_id}.VQSR.snp.sort.vcf.gz ${sample_id}.VQSR.snp.sort.vcf.gz.tbi
    rm -f ${sample_id}.VQSR.indel.sort.vcf.gz ${sample_id}.VQSR.indel.sort.vcf.gz.tbi
    """
}

process GATK_NORM_GERMLINE {
    tag { "Normalize HapCaller on $sample_id" }
    publishDir "$params.outdir/vcf", mode:'copy', pattern: "*.HapCaller.norm.vcf.*"

    input:
    tuple val(sample_id), path(vcf), path(tbi)

    output:
    tuple val(sample_id), val("HapCaller"), path("${sample_id}.HapCaller.norm.vcf.*"), emit: vcf

    script:
    """
    bcftools norm --threads $task.cpus --check-ref w --atomize --multiallelics -any -f ${params.reference} -Oz -o ${sample_id}.HapCaller.norm.vcf.gz ${vcf}
    bcftools index --threads $task.cpus -t ${sample_id}.HapCaller.norm.vcf.gz
    """
}

workflow GATK_CALL_GERM {
    take:
    gvcf

    main:
    def ch_genotyped = GATK_GENOTYPE_GVCF(gvcf)
    def ch_snp = GATK_SNP_VQSR(ch_genotyped.vcf)
    def ch_indel = GATK_INDEL_VQSR(ch_genotyped.vcf)
    def ch_merged = GATK_MERGE_VQSR(ch_snp.vcf.join(ch_indel.vcf))
    def ch_norm = GATK_NORM_GERMLINE(ch_merged.vcf)

    emit:
    vcf = ch_norm.vcf
    report = ch_snp.report.mix(ch_indel.report)
}
