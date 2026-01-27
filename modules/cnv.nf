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
    def normal_arg = ""
    if (params.cnv_normal) {
        def normal_val = params.cnv_normal.toString()
        normal_arg = (normal_val == "true" || normal_val == "self") ? "--normal" : "--normal ${params.cnv_normal}"
    }
    def targets_arg = params.cnv_targets ? "--targets ${params.bed}" : ""
    def annotate_arg = params.cnv_annotate ? "--annotate ${params.cnv_annotate}" : ""
    """
    mv $bam ${sample_id}.bam
    mv $bai ${sample_id}.bai
    cnvkit.py batch ${sample_id}.bam \
      ${normal_arg} \
      --fasta ${params.reference} \
      ${targets_arg} \
      ${annotate_arg} \
      --output-reference ${sample_id}.cnvkit.cnn \
      --output-dir ${sample_id}.cnvkit.out \
      --processes $task.cpus
    cnvkit.py scatter -s ${sample_id}.cnvkit.out/${sample_id}.cns \
      ${sample_id}.cnvkit.out/${sample_id}.cnr \
      -o ${sample_id}.cnvkit.out/${sample_id}.cnv.scatter.pdf
    cnvkit.py diagram -s ${sample_id}.cnvkit.out/${sample_id}.cns \
      ${sample_id}.cnvkit.out/${sample_id}.cnr \
      -o ${sample_id}.cnvkit.out/${sample_id}.cnv.diagram.pdf
    """
}
