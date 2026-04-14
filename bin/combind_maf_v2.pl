#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;

my ($maf, $retain_mode, $help);
GetOptions(
	"maf:s" => \$maf,
	"retain_mode:s" => \$retain_mode,
	"h"=>\$help,
);

my $usage=<<USAGE;
Discription: Combind the maf resut from all mutation callers (Mutect2, VaiDict, LoFreq, Pisces, HapCaller);
usage:perl $0 -maf <multicaller.merge.maf> [-retain_mode all|exonic]
USAGE

die $usage if (!$maf || $help);
$retain_mode = defined $retain_mode ? lc($retain_mode) : "exonic";
die "Unsupported -retain_mode: $retain_mode\n$usage" unless $retain_mode =~ /^(all|exonic)$/;

my %hash;
open IN, "$maf" or die $!;
my $header = <IN>;
chomp($header);
print "$header\tMean_VAF\tN_Callers\n";
my @header = split(/\t/, $header);

maf_input($maf, \%hash);

for my $i (sort keys %hash){
	#next if ($hash{$i}{"COMMON"} == 1);
	next if ($retain_mode eq "exonic" && $hash{$i}{"Variant_Classification"} !~ /(Splice|Mutation|Frame)/i);
	my $caller_num = scalar(@{$hash{$i}{'Mutation_Status'}});
	#	next if ($caller_num < 2);
	for my $j(@header){
		if($j eq "Mutation_Status" || $j eq "DEPTH" || $j  eq "GT" || $j eq "FILTER"  || $j eq "REF_R" || $j eq "ALT_R" || $j eq "VAF" ){
			print join("|",@{$hash{$i}{$j}})  . "\t";
		}else{
			print "$hash{$i}{$j}\t";
		}
	}
	my $merge_vaf = 0;
	for my $v (@{$hash{$i}{'VAF'}}){
		$merge_vaf +=  $v;
	}
	my $mean_vaf = sprintf("%.6f",$merge_vaf / scalar($caller_num));
	print "$mean_vaf\t$caller_num\n";
}


sub maf_input{
	my $maf = shift;
	my $h = shift;
	open IN, "$maf" or die $!;
	my $header = <IN>;
	chomp($header);
	my @header = split(/\t/, $header);
	while(<IN>){
		chomp;	
		next if (/^Hugo_Symbol/);
		my @tmp = split /\t/;
		#		next unless ($tmp[42] eq "PASS");
		my $pos = join("\t", @tmp[4..6,10..12,14]);
		for (my $n =0; $n<@tmp;$n++){
			if($header[$n] eq "Mutation_Status" || $header[$n] eq "DEPTH" || $header[$n] eq "GT" || $header[$n] eq "FILTER" 
				|| $header[$n] eq "REF_R" || $header[$n] eq "ALT_R" || $header[$n] eq "VAF" ){                   
				#Mutation_Status,DEPTH,GT,FILTER,REF_R,ALT_R,VAF
				push @{$h->{$pos}->{$header[$n]}}, $tmp[$n];
			}else{
				$h->{$pos}->{$header[$n]} = $tmp[$n];
			}
		}
	}
	close IN;

}


# 1	Hugo_Symbol
# 2	Entrez_Gene_Id
# 3	Center
# 4	NCBI_Build
# 5	Chromosome
# 6	Start_Position
# 7	End_Position
# 8	Strand
# 9	Variant_Classification
# 10	Variant_Type
# 11	Reference_Allele
# 12	Tumor_Seq_Allele1
# 13	Tumor_Seq_Allele2
# 14	dbSNP_RS
# 15	Tumor_Sample_Barcode
# 16	Mutation_Status
# 17	AAChange
# 18	Transcript_Id
# 19	TxChange
# 20	Exon_Number
# 21	DetailGene
# 22	RMSK
# 23	GeneHancer
# 24	COSMIC
# 25	ClinID
# 26	ClinDN
# 27	ClinSIG
# 28	Interpro_domain
# 29	GTEx_V6p_tissue
# 30	MCAP
# 31	REVEL
# 32	CADD
# 33	Polyphen2_HDIV
# 34	Polyphen2_HVAR
# 35	FATHMM
# 36	ExAC_ALL
# 37	ExAC_EAS
# 38	ExAC_SAS
# 39	gnomad312_AF_all
# 40	gnomad312_AF_eas
# 41	AF
# 42	DEPTH
# 43	FILTER
# 44	COMMON
# 45	REF
# 46	ALT
# 47	GT
# 48	REF_R
# 49	ALT_R
# 50	VAF
