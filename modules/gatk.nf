process GATK_rmdup {

    tag "MarkDuplicates on $sample_id"
    publishDir "$params.outdir/report", pattern: "*.metrics", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.rmdup.bam"), path("${sample_id}.rmdup.bam.bai"), emit: bam
    tuple val(sample_id), path("${sample_id}.rmdup.metrics"), emit: report

    script:
    """
    gatk --java-options "-Xmx${task.cpus * 4}g -XX:+UseParallelGC" MarkDuplicates -I ${bam[0]} -O ${sample_id}.rmdup.bam -M ${sample_id}.rmdup.metrics --REMOVE_SEQUENCING_DUPLICATES true --ASSUME_SORT_ORDER coordinate --VALIDATION_STRINGENCY LENIENT --MAX_RECORDS_IN_RAM  ${task.cpus * 4 * 500000}
    samtools index ${sample_id}.rmdup.bam
    """
}

process GATK_BQSR {
    tag "BQSR on $sample_id"

    publishDir "$params.outdir/bam", mode:'copy', pattern: "${sample_id}.bqsr.bam",  
    publishDir "$params.outdir/bam", mode:'copy', pattern: "${sample_id}.bqsr.bai"
    publishDir "$params.outdir/report", mode:'copy', pattern: "${sample_id}.recal_data.table"
    publishDir "$params.outdir/report", mode:'copy', pattern: "${sample_id}.bqsr.bam.stats"

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path( "${sample_id}.bqsr.bam"), path("${sample_id}.bqsr.bai"), emit: 'bam'
    tuple val(sample_id), path( "${sample_id}.recal_data.table"), path( "${sample_id}.bqsr.bam.stats") , emit: 'report'

    script:
    """
    gatk --java-options "-Xmx${task.cpus * 4}g" BaseRecalibrator -I ${bam[0]} -R ${params.reference} -L ${params.intervals} \
    --known-sites $params.anno/Mills_and_1000G_gold_standard.indels.hg38.vcf \
    --known-sites $params.anno/dbsnp_146.hg38.vcf \
    --known-sites $params.anno/Homo_sapiens_assembly38.known_indels.vcf \
    -O ${sample_id}.recal_data.table
    gatk --java-options "-Xmx${task.cpus * 4}g -Dsamjdk.compression_level=6" ApplyBQSR -I ${bam[0]} -R ${params.reference} --bqsr-recal-file ${sample_id}.recal_data.table -O ${sample_id}.bqsr.bam
    samtools stats -@ $task.cpus ${sample_id}.bqsr.bam > ${sample_id}.bqsr.bam.stats
    """
}

process GENOTYPE {
    tag "GenotypeGVCFs on $sample_id"
    publishDir "$params.outdir/gvcf", mode:'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.g.vcf.*")

    script:
    """
    gatk --java-options "-Xmx${task.cpus * 4}g" HaplotypeCaller -R ${params.reference} -I ${bam[0]}  -ERC GVCF -O ${sample_id}.g.vcf.gz -L ${params.intervals} -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation
    """
}

process GATK_CALL_GERM {
    tag "GATK call germline on $sample_id"
    publishDir "$params.outdir/vcf", mode:'copy', pattern: "${sample_id}.HapCaller.norm.vcf.*"
    publishDir "$params.outdir/report", mode:'copy', pattern: "*.pdf"

    input:
    tuple val(sample_id), path(gvcf)

    output:
    tuple val(sample_id), val("HapCaller"), path("${sample_id}.HapCaller.norm.vcf.*"), emit: 'vcf'
    tuple val(sample_id), path("*.pdf"), emit: 'report'

    script:
    """
    gatk --java-options "-Xmx${task.cpus * 4}g -Xms${task.cpus * 4}g" GenomicsDBImport --genomicsdb-workspace-path ${sample_id}_DB --tmp-dir ${params.tmp}  -L ${params.intervals} -V ${gvcf[0]} --merge-input-intervals
    gatk --java-options "-Xmx${task.cpus * 4}g" GenotypeGVCFs -R ${params.reference} -V gendb://${sample_id}_DB -O ${sample_id}.vcf.gz --tmp-dir ${params.tmp} --dbsnp ${params.snpdb}
    get_germline.pl ${sample_id} $task.cpus > ${sample_id}.slurm
    sh ${sample_id}.slurm
    bcftools norm --threads $task.cpus  --check-ref w  --multiallelics -any -f ${params.reference} -Oz -o ${sample_id}.HapCaller.norm.vcf.gz  ${sample_id}.VQSR.sort.vcf.gz
    bcftools index --threads $task.cpus -t ${sample_id}.HapCaller.norm.vcf.gz
    """
}
