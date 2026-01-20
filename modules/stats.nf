process BAM_DEPTH {

    tag "mosdepth_d4 on $sample_id"
    publishDir "$params.outdir/depth", mode:'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.depth.*"), emit: depth

    script:
    """
    mosdepth_d4 -t $task.cpus -b ${params.bed} -f ${params.reference} -T 10,20,50,100 ${sample_id}.depth ${bam[0]}
    """
}

process VCF_STATS {
    tag "bcftools stats on $sample_id $type"
    publishDir "$params.outdir/vcf", mode:'copy', pattern: "*.bcftools.stats"

    input:
    tuple val(sample_id), val(type), path(vcf)

    output:
    tuple val(sample_id), val(type), path("${sample_id}.${type}.bcftools.stats")

    script:
    """
    bcftools stats ${vcf[0]} -F ${params.reference} -R ${params.bed} --threads $task.cpus > ${sample_id}.${type}.bcftools.stats
    """
}
