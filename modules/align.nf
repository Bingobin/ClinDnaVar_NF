process ALING_REF {

    tag { "bwa on $sample_id" }

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), emit: bam

    script:
    """
    set -euo pipefail

    bwa mem -t $task.cpus -M -R "@RG\\tID:${sample_id}\\tLB:LIB1\\tSM:${sample_id}\\tPL:ILLUMINA" -Y ${params.reference}  ${reads[0]} ${reads[1]} \
    | samtools view -@ $task.cpus -b -o ${sample_id}.bam -

    samtools quickcheck -v ${sample_id}.bam
    samtools sort -@ $task.cpus -T ${sample_id}.samtools_sort -o ${sample_id}.sorted.bam ${sample_id}.bam
    samtools quickcheck -v ${sample_id}.sorted.bam
    samtools index -@ $task.cpus ${sample_id}.sorted.bam
    rm ${sample_id}.bam
    """
}

process RM_UMI_DUP {

    tag { "gencore on $sample_id" }
    publishDir "$params.outdir/report", mode:'copy', pattern: "*report*"

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.umi.bam"), path("${sample_id}.sorted.umi.bam.bai"), emit: bam
    tuple val(sample_id), path("${sample_id}.umi_report.*"), emit: report

    script:
    """
    gencore -i ${bam[0]}  -o ${sample_id}.umi.bam -u UMI -b ${params.bed} \
    -r ${params.reference} \
    -j ${sample_id}.umi_report.json -h ${sample_id}.umi_report.html
    samtools view -@$task.cpus -h ${sample_id}.umi.bam |  awk -F "\t" '{if(\$10~"^[ATCGN]*\$"){print \$0}}' | samtools view -@$task.cpus -Sb - | samtools sort -@ $task.cpus  - -o ${sample_id}.sorted.umi.bam
    samtools index ${sample_id}.sorted.umi.bam
    rm ${sample_id}.umi.bam
    """
}
