def cnvkitParamValue(value) {
    if (value == null) {
        return ""
    }

    def text = value.toString().trim()
    if (!text || text.equalsIgnoreCase("true") || text.equalsIgnoreCase("false") || text.equalsIgnoreCase("null")) {
        return ""
    }

    return text
}

process CNVKIT_BATCH {

    tag "cnvkit batch on $sample_id"
    publishDir "$params.outdir/cnvkit", mode:'copy'
    conda params.cnvkit_conda

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.cnvkit.out"), path("${sample_id}.cnvkit.cnn"), emit: cnv
    tuple val(sample_id), path("${sample_id}.cnvkit.out/${sample_id}*.cnr"), emit: cnr
    tuple val(sample_id), path("${sample_id}.cnvkit.out/${sample_id}*.cns"), emit: cns
    tuple val(sample_id), path("${sample_id}.cnvkit.out/${sample_id}*.call.cns"), emit: call
    tuple val(sample_id), path("${sample_id}.cnvkit.out/${sample_id}.cnv.scatter.pdf"), emit: scatter_pdf
    tuple val(sample_id), path("${sample_id}.cnvkit.out/${sample_id}.cnv.diagram.pdf"), emit: diagram_pdf

    script:
    def assayMode = params.assay_mode ? params.assay_mode.toString().toLowerCase() : "wes"
    def methodValue = cnvkitParamValue(params.cnv_method)
    def method = methodValue ? methodValue.toLowerCase() : (assayMode == "wgs" ? "wgs" : "hybrid")
    def method_arg = "--method ${method}"

    def normal_arg = ""
    if (params.cnv_normal) {
        def normalValue = params.cnv_normal.toString().trim()
        if (normalValue.equalsIgnoreCase("true") || normalValue.equalsIgnoreCase("self")) {
            normal_arg = "--normal"
        } else {
            normalValue = cnvkitParamValue(params.cnv_normal)
            normal_arg = normalValue ? "--normal ${normalValue}" : ""
        }
    }

    def targetsValue = cnvkitParamValue(params.cnv_targets)
    if (!targetsValue) {
        targetsValue = cnvkitParamValue(params.bed)
    }
    def targets_arg = (method == "wgs" || !targetsValue) ? "" : "--targets ${targetsValue}"
    def accessValue = cnvkitParamValue(params.cnv_access)
    if (!accessValue && method == "wgs") {
        accessValue = cnvkitParamValue(params.bed)
    }
    def access_arg = accessValue ? "--access ${accessValue}" : ""
    def annotateValue = cnvkitParamValue(params.cnv_annotate)
    def annotate_arg = annotateValue ? "--annotate ${annotateValue}" : ""
    """
    mv $bam ${sample_id}.bam
    mv $bai ${sample_id}.bai
    cnvkit.py batch ${sample_id}.bam \
      ${normal_arg} \
      ${method_arg} \
      --fasta ${params.reference} \
      ${targets_arg} \
      ${access_arg} \
      ${annotate_arg} \
      --output-reference ${sample_id}.cnvkit.cnn \
      --output-dir ${sample_id}.cnvkit.out \
      --processes $task.cpus
    cnvkit.py call ${sample_id}.cnvkit.out/${sample_id}.cns \
      -o ${sample_id}.cnvkit.out/${sample_id}.call.cns
    cnvkit.py scatter -s ${sample_id}.cnvkit.out/${sample_id}.call.cns \
      ${sample_id}.cnvkit.out/${sample_id}.cnr \
      -o ${sample_id}.cnvkit.out/${sample_id}.cnv.scatter.pdf
    cnvkit.py diagram -s ${sample_id}.cnvkit.out/${sample_id}.call.cns \
      ${sample_id}.cnvkit.out/${sample_id}.cnr \
      -o ${sample_id}.cnvkit.out/${sample_id}.cnv.diagram.pdf
    """
}
