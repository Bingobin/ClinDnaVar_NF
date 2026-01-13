#!/bin/bash

PREFIX=$1
GVCF=$2
INT=$3


GENOME=/lustre/home/acct-medkkw/medlyb/database/annotation/gatk_ann/hg38/bwaindex2/Homo_sapiens_assembly38.fasta
TMP=/lustre/home/acct-medkkw/medlyb/tmp
ANNO=/lustre/home/acct-medkkw/medlyb/database/annotation/gatk_ann/hg38/hg38bundle


#########################################
echo "Mission $PREFIX GenotypeGVCFs starts at  `date`" &&\
gatk --java-options "-Xmx14g -Xms14g" GenomicsDBImport --genomicsdb-workspace-path ${PREFIX}_DB --tmp-dir $TMP  -L $INT -V $GVCF --merge-input-intervals &&\
gatk --java-options "-Xmx14g" GenotypeGVCFs -R $GENOME -V gendb://${PREFIX}_DB -O ${PREFIX}.vcf.gz --tmp-dir $TMP --dbsnp $ANNO/../dbsnp_155.hg38.vcf.gz &&\
echo "Mission completes at `date`"
########################################
echo "Mission $PREFIX SNP_recal starts at  `date`" &&\
gatk SelectVariants -R $GENOME -V ${PREFIX}.vcf.gz --tmp-dir $TMP -O ${PREFIX}.raw.snp.vcf.gz --exclude-non-variants --select-type-to-include SNP &&\
gatk VariantRecalibrator -R $GENOME -V ${PREFIX}.raw.snp.vcf.gz \
	 --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $ANNO/hapmap_3.3.hg38.vcf.gz  \
	 --resource:omni,known=false,training=true,truth=false,prior=12.0 $ANNO/1000G_omni2.5.hg38.vcf.gz \
	 --resource:1000G,known=false,training=true,truth=false,prior=10.0 $ANNO/1000G_phase1.snps.high_confidence.hg38.vcf \
	 --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $ANNO/dbsnp_146.hg38.vcf \
	 -an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum \
	 -mode SNP \
	 --output ${PREFIX}.recalibrate_SNP.recal \
	 --tranches-file ${PREFIX}.recalibrate_SNP.tranches \
	 --rscript-file ${PREFIX}.recalibrate_SNP_plots.R &&\
gatk ApplyVQSR -R $GENOME -V ${PREFIX}.raw.snp.vcf.gz --tranches-file ${PREFIX}.recalibrate_SNP.tranches --recal-file ${PREFIX}.recalibrate_SNP.recal -O ${PREFIX}.VQSR.snp.vcf.gz -mode SNP &&\
rm ${PREFIX}.raw.snp.vcf.gz ${PREFIX}.raw.snp.vcf.gz.tbi &&\
echo "Mission completes at `date`"
########################################
echo "Mission $PREFIX INDEL_recal starts at  `date`" &&\
gatk SelectVariants -R $GENOME -V ${PREFIX}.vcf.gz --tmp-dir $TMP -O ${PREFIX}.raw.indel.vcf.gz --exclude-non-variants --select-type-to-include INDEL &&\
gatk VariantRecalibrator -R $GENOME -V ${PREFIX}.raw.indel.vcf.gz \
		 --resource:mills,known=true,training=true,truth=true,prior=12.0 $ANNO/Mills_and_1000G_gold_standard.indels.hg38.vcf  \
		 --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $ANNO/dbsnp_146.hg38.vcf \
		 -an SOR -an MQ -an MQRankSum -an ReadPosRankSum --max-gaussians 4 \
		 -mode INDEL \
		 -O ${PREFIX}.recalibrate_INDEL.recal \
		 --tranches-file ${PREFIX}.recalibrate_INDEL.tranches \
		 --rscript-file ${PREFIX}.recalibrate_INDEL_plots.R &&\
gatk ApplyVQSR -R $GENOME -V ${PREFIX}.raw.indel.vcf.gz --tranches-file ${PREFIX}.recalibrate_INDEL.tranches --recal-file ${PREFIX}.recalibrate_INDEL.recal -O ${PREFIX}.VQSR.indel.vcf.gz -mode INDEL &&\
#rm ${PREFIX}.raw.indel.vcf.gz  ${PREFIX}.raw.indel.vcf.gz.tbi &&\
echo "Mission completes at `date`"
#######################################
echo "Mission $PREFIX Merge starts at  `date`" &&\
gatk SortVcf -I ${PREFIX}.VQSR.indel.vcf.gz -O ${PREFIX}.VQSR.indel.sort.vcf.gz &&\
gatk SortVcf -I ${PREFIX}.VQSR.snp.vcf.gz -O ${PREFIX}.VQSR.snp.sort.vcf.gz &&\
gatk MergeVcfs -I ${PREFIX}.VQSR.snp.sort.vcf.gz -I ${PREFIX}.VQSR.indel.sort.vcf.gz -O ${PREFIX}.VQSR.sort.vcf.gz &&\
rm ${PREFIX}.VQSR.indel.vcf.gz ${PREFIX}.VQSR.indel.vcf.gz.tbi &&\
rm ${PREFIX}.VQSR.indel.sort.vcf.gz   ${PREFIX}.VQSR.indel.sort.vcf.gz.tbi &&\
rm ${PREFIX}.VQSR.snp.vcf.gz ${PREFIX}.VQSR.snp.vcf.gz.tbi &&\
rm ${PREFIX}.VQSR.snp.sort.vcf.gz ${PREFIX}.VQSR.snp.sort.vcf.gz.tbi &&\
rm -r ${PREFIX}_DB &&\
#rm ${PREFIX}.recalibrate_INDEL.recal ${PREFIX}.recalibrate_INDEL.recal.idx &&\
rm ${PREFIX}.recalibrate_SNP.recal ${PREFIX}.recalibrate_SNP.recal.idx &&\
rm ${PREFIX}.vcf.gz ${PREFIX}.vcf.gz.tbi &&\
echo "Mission completes at `date`"
######################################
#echo "Mission $PREFIX ANNOVAR starts at  `date`" &&\
#bcftools annotate --threads 2 -a $ANNO/../dbsnp_155.hg38.vcf.gz -c INFO/SSR,INFO/GENEINFO,INFO/PSEUDOGENEINFO,INFO/SAO,INFO/SSR,INFO/VC,INFO/PM,INFO/NSF,INFO/NSM,INFO/NSN,INFO/SYN,INFO/U3,INFO/U5,INFO/ASS,INFO/DSS,INFO/INT,INFO/R3,INFO/R5,INFO/GNO,INFO/PUB,INFO/FREQ,INFO/COMMON ${PREFIX}.VQSR.sort.vcf.gz -o  ${PREFIX}.VQSR.sort.snpdb155.vcf &&\
#		 rm ${PREFIX}.VQSR.sort.vcf.gz ${PREFIX}.VQSR.sort.vcf.gz.tbi &&\
#perl /lustre/home/acct-medkkw/medlyb/soft/annovar/table_annovar.pl ${PREFIX}.VQSR.sort.snpdb155.vcf \
#	 /lustre/home/acct-medkkw/medlyb/soft/annovar/humandb/ \
#	 -buildver hg38 \
#	 -protocol refGene,rmsk,gff3,cosmic87,mcap,revel,clinvar_20221231,exac03,dbnsfp35a,gnomad312_genome \
#	 --gff3dbfile hg38_genehancer.txt \
#	 -operation g,r,r,f,f,f,f,f,f,f, \
#	 -vcfinput --remove --polish \
#	 --outfile ${PREFIX}.VQSR.sort.snpdb155 &&\
#rm ${PREFIX}.VQSR.sort.snpdb155.vcf ${PREFIX}.VQSR.sort.snpdb155.avinput &&\
#perl ../annovar2maf_germline_sin.pl ${PREFIX}.VQSR.sort.snpdb155.hg38_multianno.txt  ${PREFIX} > ${PREFIX}.VQSR.sort.snpdb155.hg38_multianno.maf &&\
#	 gzip ${PREFIX}.VQSR.sort.snpdb155.hg38_multianno.txt &&\
#	 bcftools view -Oz -l 9  -o ${PREFIX}.VQSR.sort.snpdb155.hg38_multianno.vcf.gz ${PREFIX}.VQSR.sort.snpdb155.hg38_multianno.vcf && rm ${PREFIX}.VQSR.sort.snpdb155.hg38_multianno.vcf &&\
#	 bcftools index -t ${PREFIX}.VQSR.sort.snpdb155.hg38_multianno.vcf.gz
#echo "Mission completes at `date`"
#####################################
