process FASTQC {
    cpus 2
    tag { "fastqc in $sample_id" }
    publishDir "$params.outdir/report", mode:'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.{html,zip}")

    script:
    """
    if [[ "${reads[0]}" == *.gz ]]; then
        ext=".fastq.gz"
    else
        ext=".fastq"
    fi
    mv ${reads[0]} ${sample_id}_R1\${ext}
    mv ${reads[1]} ${sample_id}_R2\${ext}
    fastqc -t $task.cpus ${sample_id}_R1\${ext} ${sample_id}_R2\${ext}
    """
}

process FastpFilter {

    tag { "fastp on $sample_id" }
    publishDir "$params.outdir/report", mode:'copy', pattern: "*report*"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_clean_R*.fq.gz"), emit: fastq
    tuple val(sample_id), path("${sample_id}.filter_report.*"), emit: report

    script:
    def use_umi = params.assay_mode?.toString()?.toLowerCase() == 'panel_umi' || params.use_umi?.toString()?.toBoolean()
    def umi_opts = use_umi ? '-U --umi_loc per_read --umi_prefix UMI --umi_len 4' : ''
    """
    fastp -i ${reads[0]} -I ${reads[1]} -o ${sample_id}_clean_R1.fq.gz -O ${sample_id}_clean_R2.fq.gz ${umi_opts} -w $task.cpus -j ${sample_id}.filter_report.json -h ${sample_id}.filter_report.html
    """
}
