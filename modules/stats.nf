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

def optionArg(option, value) {
    def text = paramValue(value)
    return text ? "${option} ${text}" : ""
}

process BAM_DEPTH {

    tag { "mosdepth_d4 on $sample_id" }
    publishDir "$params.outdir/depth", mode:'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.depth.*"), emit: depth

    script:
    def assayMode = params.assay_mode ? params.assay_mode.toString().toLowerCase() : "wes"
    def bedArg = optionArg("-b", params.bed)
    def thresholds = params.depth_thresholds ? params.depth_thresholds : (assayMode.startsWith("panel") ? "100,200,500,1000" : (assayMode == "wgs" ? "10,20,30" : "10,20,50,100"))
    """
    mosdepth_d4 -t $task.cpus ${bedArg} -f ${params.reference} -T ${thresholds} ${sample_id}.depth ${bam[0]}
    """
}

process VCF_STATS {
    tag { "bcftools stats on $sample_id $type" }
    publishDir "$params.outdir/vcf", mode:'copy', pattern: "*.bcftools.stats"

    input:
    tuple val(sample_id), val(type), path(vcf)

    output:
    tuple val(sample_id), val(type), path("${sample_id}.${type}.bcftools.stats")

    script:
    def regionsArg = optionArg("-R", params.bed)
    """
    bcftools stats ${vcf[0]} -F ${params.reference} ${regionsArg} --threads $task.cpus > ${sample_id}.${type}.bcftools.stats
    """
}
