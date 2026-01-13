#!/usr/bin/perl

use warnings;
use strict;

my $namepre = shift or die "perl $0 [Sample prefix]: $!";
my $thread = shift or die $!;
my $memory = "-Xmx" . $thread * 4 . "g";

print <<EOF
#!/bin/bash

#SBATCH -J $namepre\_germline
#SBATCH -p cpu
#SBATCH -n $thread
#SBATCH -o $namepre.germline.slurm.o.%j
#SBATCH -e $namepre.germline.slurm.e.%j

uname -n &&\\

GENOME=/lustre/home/acct-medkkw/medlyb/database/annotation/gatk_ann/hg38/bwaindex2/Homo_sapiens_assembly38.fasta
TMP=/lustre/home/acct-medkkw/medlyb/tmp
ANNO=/lustre/home/acct-medkkw/medlyb/database/annotation/gatk_ann/hg38/hg38bundle
INT=/lustre/home/acct-medty/medty-c/project/02.WES_infect/03.GATK_call_split/bak.GATK_call/SHPHC_infect.intervals

PREFIX=$namepre

########################################
echo "Mission $namepre SNP_recal starts at  `date`" &&\\
gatk --java-options "$memory" SelectVariants -R \$GENOME -V \${PREFIX}.vcf.gz --tmp-dir \$TMP -O \${PREFIX}.raw.snp.vcf.gz --exclude-non-variants --select-type-to-include SNP &&\\
gatk --java-options "$memory" VariantRecalibrator -R \$GENOME -V \${PREFIX}.raw.snp.vcf.gz \\
	 --resource:hapmap,known=false,training=true,truth=true,prior=15.0 \$ANNO/hapmap_3.3.hg38.vcf.gz  \\
	 --resource:omni,known=false,training=true,truth=false,prior=12.0 \$ANNO/1000G_omni2.5.hg38.vcf.gz \\
	 --resource:1000G,known=false,training=true,truth=false,prior=10.0 \$ANNO/1000G_phase1.snps.high_confidence.hg38.vcf \\
	 --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 \$ANNO/dbsnp_146.hg38.vcf \\
	 -an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum \\
	 -mode SNP \\
	 --output \${PREFIX}.recalibrate_SNP.recal \\
	 --tranches-file \${PREFIX}.recalibrate_SNP.tranches \\
	 --rscript-file \${PREFIX}.recalibrate_SNP_plots.R &&\\
gatk --java-options "$memory" ApplyVQSR -R \$GENOME -V \${PREFIX}.raw.snp.vcf.gz --tranches-file \${PREFIX}.recalibrate_SNP.tranches --recal-file \${PREFIX}.recalibrate_SNP.recal -O \${PREFIX}.VQSR.snp.vcf.gz -mode SNP &&\\
rm \${PREFIX}.raw.snp.vcf.gz \${PREFIX}.raw.snp.vcf.gz.tbi &&\\
echo "Mission completes at `date`"
########################################
echo "Mission $namepre INDEL_recal starts at  `date`" &&\\
gatk --java-options "$memory" SelectVariants -R \$GENOME -V \${PREFIX}.vcf.gz --tmp-dir \$TMP -O \${PREFIX}.raw.indel.vcf.gz --exclude-non-variants --select-type-to-include INDEL &&\\
gatk --java-options "$memory" VariantRecalibrator -R \$GENOME -V \${PREFIX}.raw.indel.vcf.gz \\
		 --resource:mills,known=true,training=true,truth=true,prior=12.0 \$ANNO/Mills_and_1000G_gold_standard.indels.hg38.vcf  \\
		 --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 \$ANNO/dbsnp_146.hg38.vcf \\
		 -an SOR -an MQ -an MQRankSum -an ReadPosRankSum --max-gaussians 4 \\
		 -mode INDEL \\
		 -O \${PREFIX}.recalibrate_INDEL.recal \\
		 --tranches-file \${PREFIX}.recalibrate_INDEL.tranches \\
		 --rscript-file \${PREFIX}.recalibrate_INDEL_plots.R &&\\
gatk --java-options "$memory" ApplyVQSR -R \$GENOME -V \${PREFIX}.raw.indel.vcf.gz --tranches-file \${PREFIX}.recalibrate_INDEL.tranches --recal-file \${PREFIX}.recalibrate_INDEL.recal -O \${PREFIX}.VQSR.indel.vcf.gz -mode INDEL &&\\
rm \${PREFIX}.raw.indel.vcf.gz  \${PREFIX}.raw.indel.vcf.gz.tbi &&\\
echo "Mission completes at `date`"
#######################################
echo "Mission $namepre Merge starts at  `date`" &&\\
gatk --java-options "$memory" SortVcf -I \${PREFIX}.VQSR.indel.vcf.gz -O \${PREFIX}.VQSR.indel.sort.vcf.gz &&\\
gatk --java-options "$memory" SortVcf -I \${PREFIX}.VQSR.snp.vcf.gz -O \${PREFIX}.VQSR.snp.sort.vcf.gz &&\\
gatk --java-options "$memory" MergeVcfs -I \${PREFIX}.VQSR.snp.sort.vcf.gz -I \${PREFIX}.VQSR.indel.sort.vcf.gz -O \${PREFIX}.VQSR.sort.vcf.gz &&\\
rm \${PREFIX}.VQSR.indel.vcf.gz \${PREFIX}.VQSR.indel.vcf.gz.tbi &&\\
rm \${PREFIX}.VQSR.indel.sort.vcf.gz   \${PREFIX}.VQSR.indel.sort.vcf.gz.tbi &&\\
rm \${PREFIX}.VQSR.snp.vcf.gz \${PREFIX}.VQSR.snp.vcf.gz.tbi &&\\
rm \${PREFIX}.VQSR.snp.sort.vcf.gz \${PREFIX}.VQSR.snp.sort.vcf.gz.tbi &&\\
rm -r \${PREFIX}_DB &&\\
rm \${PREFIX}.recalibrate_INDEL.recal \${PREFIX}.recalibrate_INDEL.recal.idx &&\\
rm \${PREFIX}.recalibrate_SNP.recal \${PREFIX}.recalibrate_SNP.recal.idx &&\\
rm \${PREFIX}.vcf.gz \${PREFIX}.vcf.gz.tbi &&\\
echo "Mission completes at `date`"
#####################################
EOF
