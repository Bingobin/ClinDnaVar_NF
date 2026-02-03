process MUT_ANNOTATE {
    cpus 1
    tag "annovar on $sample_id"
    publishDir "$params.outdir/maf", mode:'copy', pattern: "*.snpdb155.hg38_multianno.maf.gz"
    publishDir "$params.outdir/maf", mode:'copy', pattern: "*.snpdb155.hg38_multianno.vcf.*"

    input:
    tuple val(sample_id), val(type), path(vcf)

    output:
    tuple val(sample_id), path("${sample_id}.${type}.snpdb155.hg38_multianno.maf.gz"), emit: maf
    tuple val(sample_id), path("${sample_id}.${type}.snpdb155.hg38_multianno.vcf.*"), emit: vcf

    script:
    """
    bcftools annotate --threads $task.cpus -a ${params.snpdb} -c INFO/RS,INFO/SSR,INFO/GENEINFO,INFO/PSEUDOGENEINFO,INFO/SAO,INFO/SSR,INFO/VC,INFO/PM,INFO/NSF,INFO/NSM,INFO/NSN,INFO/SYN,INFO/U3,INFO/U5,INFO/ASS,INFO/DSS,INFO/INT,INFO/R3,INFO/R5,INFO/GNO,INFO/PUB,INFO/FREQ,INFO/COMMON ${vcf[0]} -o   ${sample_id}.${type}.snpdb155.vcf
    perl ${params.annovar_path}/table_annovar.pl  ${sample_id}.${type}.snpdb155.vcf \
        ${params.annovar_db}/ \
        -buildver hg38 \
        -protocol ${params.annovar_protocol} \
        --gff3dbfile ${params.annovar_gff3db} \
        -operation ${params.annovar_operation} \
        -vcfinput --remove --polish \
        --outfile  ${sample_id}.${type}.snpdb155
    rm  ${sample_id}.${type}.snpdb155.vcf  ${sample_id}.${type}.snpdb155.avinput &&\
    annovar2maf_multitype.pl  ${sample_id}.${type}.snpdb155.hg38_multianno.txt $sample_id $type ${params.source_id} > ${sample_id}.${type}.snpdb155.hg38_multianno.maf
    bcftools view --threads $task.cpus -Oz -l 9  -o  ${sample_id}.${type}.snpdb155.hg38_multianno.vcf.gz  ${sample_id}.${type}.snpdb155.hg38_multianno.vcf && rm  ${sample_id}.${type}.snpdb155.hg38_multianno.vcf
    bcftools index --threads $task.cpus -t  ${sample_id}.${type}.snpdb155.hg38_multianno.vcf.gz
    gzip ${sample_id}.${type}.snpdb155.hg38_multianno.maf
    """
}

process MAF_COMBIND {
    tag "combind_maf on $sample_id"
    publishDir "$params.outdir/maf", mode:'copy'

    input:
    tuple val(sample_id), path(maf)

    output:
    tuple val(sample_id), path("${sample_id}.multicaller.combind.maf")

    script:
    """
    zcat ${maf[0]} | head -n 1 > maf.header
    zcat *.maf.gz | grep -v "Hugo_Symbol" > merge.tmp.file
    cat maf.header merge.tmp.file > ${sample_id}.multicaller.merge.maf && rm merge.tmp.file maf.header
    maf_sort_by_pos.pl ${sample_id}.multicaller.merge.maf > a && mv a ${sample_id}.multicaller.merge.maf
    combind_maf_v2.pl -maf ${sample_id}.multicaller.merge.maf > ${sample_id}.multicaller.combind.maf
    maf_sort_by_pos.pl ${sample_id}.multicaller.combind.maf > a && mv a ${sample_id}.multicaller.combind.maf
    """
}
