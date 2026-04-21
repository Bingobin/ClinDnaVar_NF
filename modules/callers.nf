process DeepVariant_CALL {
    tag "DeepVariant on $sample_id"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), emit: vcf
    tuple val(sample_id), path("${sample_id}.g.vcf.gz"), emit: gvcf

    script:
    """
    deepvariant --model_type=WES --ref=${params.reference} --reads=${bam[0]} --output_vcf=${sample_id}.vcf.gz --output_gvcf= ${sample_id}.g.vcf.gz --num_shards=$task.cpus --regions ${params.bed}
    """
}

process Mutect2_Call {
    tag "mutect2 on $sample_id"
    publishDir "$params.outdir/vcf", mode:'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), val("Mutect2"), path("${sample_id}.mutect2.norm.vcf.*")

    script:
    def intervalsArg = params.intervals ? "-L ${params.intervals}" : ""
    """
    gatk --java-options "-Xmx${task.cpus * 4}g" Mutect2 --native-pair-hmm-threads $task.cpus -R ${params.reference} -I $bam -O ${sample_id}.mutect2.vcf.gz ${intervalsArg} --read-index $bai --f1r2-tar-gz mutect2.f1r2.tar.gz --max-reads-per-alignment-start 0
    gatk --java-options "-Xmx${task.cpus * 4}g" LearnReadOrientationModel -I  mutect2.f1r2.tar.gz -O mutect2.atrifact_prior.tar.gz
    gatk --java-options "-Xmx${task.cpus * 4}g" FilterMutectCalls -R ${params.reference} -V ${sample_id}.mutect2.vcf.gz -ob-priors mutect2.atrifact_prior.tar.gz -O ${sample_id}.mutect2.filtered.vcf.gz
    bcftools norm --threads $task.cpus  --check-ref w --atomize --multiallelics -any -f ${params.reference} -Oz -o ${sample_id}.mutect2.norm.vcf.gz  ${sample_id}.mutect2.filtered.vcf.gz
    bcftools index --threads $task.cpus -t ${sample_id}.mutect2.norm.vcf.gz
    """
}

process Pindel_Split_BED {
    tag "split bed file for Pindel"

    input:
    tuple val(sample_id), path(bam), path(bai)
    path(bed)

    output:
    tuple val(sample_id), path(bam), path(bai), path("splitBed.*.bed")

    script:
    """
    bedtools makewindows -b $bed -w 20250 -s 20000 > ${bed}.windows.bed
    split -d --additional-suffix .bed -n l/${task.cpus}  ${bed}.windows.bed splitBed.
    """
}

process Pindel_Call {
    cpus 1
    tag "pindel on $sample_id  ${file(bed).getBaseName()}"

    input:
    tuple val(sample_id), path(bam), path(bai), path(bed)

    output:
    tuple val(sample_id), path("*.pindel_result.txt")

    script:
    """
    echo -e "${bam[0]}\t400\t$sample_id" >  pindel_config.txt
    pindel -T 4 -i pindel_config.txt -f ${params.reference} -o pindel_$sample_id -M 3 -j $bed
    cat *_D *_SI *_TD  *_INV  | grep "ChrID" > ${file(bed).getBaseName()}.pindel_result.txt
    """
}

process Pindel_Merge {
    tag "pindel merge on $sample_id"
    publishDir "$params.outdir/vcf", mode:'copy'

    input:
    tuple val(sample_id), path(result)

    output:
    tuple val(sample_id), val("Pindel"), path("${sample_id}.pindel.norm.vcf.*")

    script:
    def ref = file(params.reference)
    """
    cat *.pindel_result.txt > ${sample_id}.pindel_result.merge.txt
    pindel2vcf -p ${sample_id}.pindel_result.merge.txt -r ${params.reference} -R ${ref.getSimpleName()} -v ${sample_id}.pindel.vcf -d 20240628 -G
    gatk --java-options "-Xmx${task.cpus * 4}g" MergeVcfs -I ${sample_id}.pindel.vcf -O ${sample_id}.pindel.vcf.gz -D ${params.ref_dict}
    bcftools norm  --threads $task.cpus --check-ref w --atomize --multiallelics -any -f ${params.reference} -Oz -o ${sample_id}.pindel.norm.vcf.gz  ${sample_id}.pindel.vcf.gz
    bcftools index --threads $task.cpus -t ${sample_id}.pindel.norm.vcf.gz
    """
}

process LoFreq_Call {
    tag "lofreq on $sample_id"
    publishDir "$params.outdir/vcf", mode:'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), val("LoFreq"), path("${sample_id}.lofreq.norm.vcf.*")

    script:
    """
    lofreq indelqual --dindel -f ${params.reference} -o ${sample_id}.indel.bam ${bam}
    samtools index ${sample_id}.indel.bam
    lofreq call-parallel --pp-threads $task.cpus --call-indels -f ${params.reference} -o ${sample_id}.lofreq.vcf ${sample_id}.indel.bam --bed ${params.bed}
    lofreq_reformat.pl  ${sample_id}.lofreq.vcf ${sample_id} > ${sample_id}.lofreq.reformat.vcf
    gatk --java-options "-Xmx${task.cpus * 4}g" MergeVcfs -I ${sample_id}.lofreq.reformat.vcf -O ${sample_id}.lofreq.reformat.vcf.gz -D ${params.ref_dict}
    bcftools norm  --threads $task.cpus --check-ref w --atomize --multiallelics -any -f ${params.reference} -Oz  -o ${sample_id}.lofreq.norm.vcf.gz  ${sample_id}.lofreq.reformat.vcf.gz
    bcftools index --threads $task.cpus  -t ${sample_id}.lofreq.norm.vcf.gz
    """

}

process VarDict_Call {
    tag "vardict on $sample_id"
    publishDir "$params.outdir/vcf", mode:'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), val("VarDict"), path("${sample_id}.vardict.norm.vcf.*")

    script:
    """
    vardict-java -U -G ${params.reference} -f 0.0001 -N $sample_id  -b ${bam} -deldupvar -Q 10  -c 1 -S 2 -E 3 -g 4 -F 0x704 -th $task.cpus  -fisher ${params.bed} | var2vcf_valid.pl -N $sample_id -E -f 0.0001 > ${sample_id}.vardict.vcf
    gatk --java-options "-Xmx${task.cpus * 4}g" MergeVcfs -I ${sample_id}.vardict.vcf -O ${sample_id}.vardict.vcf.gz -D ${params.ref_dict}
    bcftools filter --threads $task.cpus -e "((FMT/AF[0] * FMT/DP < 6) && ((INFO/MQ < 55.0 && INFO/NM > 1.0) || (INFO/MQ < 60.0 && INFO/NM > 3.0) || (FMT/DP < 6500) || (INFO/QUAL < 27)))" -Oz  -o ${sample_id}.vardict.f.vcf.gz ${sample_id}.vardict.vcf.gz
    bcftools index ${sample_id}.vardict.f.vcf.gz
    bcftools norm  --threads $task.cpus --check-ref w --atomize --multiallelics -any -f ${params.reference} -Oz -o ${sample_id}.vardict.norm.vcf.gz  ${sample_id}.vardict.f.vcf.gz
    bcftools index --threads $task.cpus  -t ${sample_id}.vardict.norm.vcf.gz
    """
}

process Pisces_Call {
    tag "pisces on $sample_id"
    publishDir "$params.outdir/vcf", mode:'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), val("Pisces"), path("${sample_id}.pisces.norm.vcf.*")

    script:
    """
    mv $bam ${sample_id}.bam &&  mv $bai ${sample_id}.bam.bai
    pisces -o pisces_result -bam ${sample_id}.bam -g ${file(params.reference).getParent()} -i ${params.bed} -t $task.cpus --minvf 0.0005 --callmnvs false  --mindpfilter 500 --mindp 5 --threadbychr true --ssfilter false --minvq 0 --vqfilter 20  --minbq 30 --reportnocalls true --gvcf false
    bcftools filter -i 'FILTER=="PASS"' pisces_result/${sample_id}.vcf -o ${sample_id}.pass.vcf
    bcftools norm  --threads $task.cpus --check-ref w --atomize --multiallelics -any -f ${params.reference} -Oz -o ${sample_id}.pisces.norm.vcf.gz  ${sample_id}.pass.vcf
    bcftools index --threads $task.cpus  -t ${sample_id}.pisces.norm.vcf.gz
    """
}
