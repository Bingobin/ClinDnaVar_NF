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

def optionArg(option, value) {
    def text = paramValue(value)
    return text ? "${option} ${text}" : ""
}

def gatkJavaOptions(task) {
    def heapGb = Math.max(1, Math.floor(task.memory.toGiga() * 0.8) as int)
    return "-Xmx${heapGb}g"
}

def vardictJavaOptions(task) {
    def heapGb = Math.max(1, Math.floor(task.memory.toGiga() * 0.8) as int)
    return "-Xmx${heapGb}g"
}

process DeepVariant_CALL {
    tag "DeepVariant on $sample_id"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), emit: vcf
    tuple val(sample_id), path("${sample_id}.g.vcf.gz"), emit: gvcf

    script:
    def regionsArg = optionArg("--regions", params.bed)
    def modelType = params.assay_mode && params.assay_mode.toString().toLowerCase() == "wgs" ? "WGS" : "WES"
    """
    deepvariant --model_type=${modelType} --ref=${params.reference} --reads=${bam[0]} --output_vcf=${sample_id}.vcf.gz --output_gvcf=${sample_id}.g.vcf.gz --num_shards=$task.cpus ${regionsArg}
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
    def intervalsArg = gatkIntervalsArg(params.intervals)
    def javaOptions = gatkJavaOptions(task)
    """
    gatk --java-options "${javaOptions}" Mutect2 --native-pair-hmm-threads $task.cpus -R ${params.reference} -I $bam -O ${sample_id}.mutect2.vcf.gz ${intervalsArg} --read-index $bai --f1r2-tar-gz mutect2.f1r2.tar.gz --max-reads-per-alignment-start 0
    gatk --java-options "${javaOptions}" LearnReadOrientationModel -I  mutect2.f1r2.tar.gz -O mutect2.atrifact_prior.tar.gz
    gatk --java-options "${javaOptions}" FilterMutectCalls -R ${params.reference} -V ${sample_id}.mutect2.vcf.gz -ob-priors mutect2.atrifact_prior.tar.gz -O ${sample_id}.mutect2.filtered.vcf.gz
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
    def javaOptions = gatkJavaOptions(task)
    """
    cat *.pindel_result.txt > ${sample_id}.pindel_result.merge.txt
    pindel2vcf -p ${sample_id}.pindel_result.merge.txt -r ${params.reference} -R ${ref.getSimpleName()} -v ${sample_id}.pindel.vcf -d 20240628 -G
    gatk --java-options "${javaOptions}" MergeVcfs -I ${sample_id}.pindel.vcf -O ${sample_id}.pindel.vcf.gz -D ${params.ref_dict}
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
    def bedArg = optionArg("--bed", params.bed)
    def javaOptions = gatkJavaOptions(task)
    """
    lofreq indelqual --dindel -f ${params.reference} -o ${sample_id}.indel.bam ${bam}
    samtools index ${sample_id}.indel.bam
    lofreq call-parallel --pp-threads $task.cpus --call-indels -f ${params.reference} -o ${sample_id}.lofreq.vcf ${sample_id}.indel.bam ${bedArg}
    lofreq_reformat.pl  ${sample_id}.lofreq.vcf ${sample_id} > ${sample_id}.lofreq.reformat.vcf
    awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} (\$4 !~ /[MRWSYKVHDBNmrwsykvhdbn]/ && \$5 !~ /[MRWSYKVHDBNmrwsykvhdbn]/) {print}' ${sample_id}.lofreq.reformat.vcf > ${sample_id}.lofreq.reformat.filtered.vcf && mv ${sample_id}.lofreq.reformat.filtered.vcf ${sample_id}.lofreq.reformat.vcf
    gatk --java-options "${javaOptions}" MergeVcfs -I ${sample_id}.lofreq.reformat.vcf -O ${sample_id}.lofreq.reformat.vcf.gz -D ${params.ref_dict}
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
    def bedValue = paramValue(params.bed)
    def javaOptions = gatkJavaOptions(task)
    def vardictOptions = vardictJavaOptions(task)
    """
    JAVA_TOOL_OPTIONS="${vardictOptions}" vardict-java -U -G ${params.reference} -f 0.0001 -N $sample_id  -b ${bam} -deldupvar -Q 10  -c 1 -S 2 -E 3 -g 4 -F 0x704 -th $task.cpus  -fisher ${bedValue} | var2vcf_valid.pl -N $sample_id -E -f 0.0001 > ${sample_id}.vardict.vcf
    gatk --java-options "${javaOptions}" MergeVcfs -I ${sample_id}.vardict.vcf -O ${sample_id}.vardict.vcf.gz -D ${params.ref_dict}
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
    def intervalsArg = optionArg("-i", params.bed)
    """
    mv $bam ${sample_id}.bam &&  mv $bai ${sample_id}.bam.bai
    pisces -o pisces_result -bam ${sample_id}.bam -g ${file(params.reference).getParent()} ${intervalsArg} -t $task.cpus --minvf 0.0005 --callmnvs false  --mindpfilter 500 --mindp 5 --threadbychr true --ssfilter false --minvq 0 --vqfilter 20  --minbq 30 --reportnocalls true --gvcf false
    bcftools filter -i 'FILTER=="PASS"' pisces_result/${sample_id}.vcf -o ${sample_id}.pass.vcf
    bcftools norm  --threads $task.cpus --check-ref w --atomize --multiallelics -any -f ${params.reference} -Oz -o ${sample_id}.pisces.norm.vcf.gz  ${sample_id}.pass.vcf
    bcftools index --threads $task.cpus  -t ${sample_id}.pisces.norm.vcf.gz
    """
}
