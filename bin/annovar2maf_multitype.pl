#!/usr/bin/perl

use warnings;
use strict;

my $file = shift or die $!;
my $prefix = shift or die $!;
my $type = shift or die $!;  #Mutect2 Pindel LoFreq VarDict Pisces HapCaller
my $source = shift or die $!;

print "Hugo_Symbol\tEntrez_Gene_Id\tCenter\tNCBI_Build\tChromosome\tStart_Position\tEnd_Position\tStrand\tVariant_Classification\tVariant_Type\tReference_Allele\tTumor_Seq_Allele1\tTumor_Seq_Allele2\t";
print "dbSNP_RS\tTumor_Sample_Barcode\tMutation_Status\tAAChange\tTranscript_Id\tTxChange\tExon_Number\tDetailGene\t";
print "RMSK\tGeneHancer\tCOSMIC\tClinID\tClinDN\tClinSIG\tInterpro_domain\tGTEx_V6p_tissue\t";
print "MCAP\tREVEL\tCADD\tPolyphen2_HDIV\tPolyphen2_HVAR\tFATHMM\t";
print "ExAC_ALL\tExAC_EAS\tExAC_SAS\tgnomad312_AF_all\tgnomad312_AF_eas\t";
print "AF\tDEPTH\tFILTER\tCOMMON\tREF\tALT\tGT\tREF_R\tALT_R\tVAF\n";

my %sam;
$sam{$prefix} = $source; #Center

my $entrez_gene_id = "NA";
my $ncbi_build = "hg38";
my $mutation_status = $type;
my $strand = "+";

my %variant = (
		"frameshift deletion" => "Frame_Shift_Del",
		"frameshift insertion" => "Frame_Shift_Ins",
		"frameshift block substitution" => "Frameshift_INDEL",
		"frameshift substitution" => "Frameshift_INDEL",
		"nonframeshift deletion" => "In_Frame_Del",
		"nonframeshift insertion" => "In_Frame_Ins",
		"nonframeshift block substitution" => "Inframe_INDEL",
		"nonframeshift substitution" => "Inframe_INDEL",
		"nonsynonymous SNV" => "Missense_Mutation",
		"stopgain" => "Nonsense_Mutation",
		"stoploss" => "Nonstop_Mutation",
		"startloss" => "Translation_Start_Site",
		"startgain" => "Unknown", 
		"synonymous SNV" => "Silent",
		"unknown" => "UNKNOWN",
		"UNKNOWN" => "UNKNOWN",
		"NA" => "UNKNOWN",
		"downstream" => "3'Flank",
		"upstream" => "5'Flank",
		"splicing" => "Splice_Site",
		"UTR3" => "3'UTR",
		"UTR5" => "5'UTR",
		"intergenic" => "IGR",
		"intronic" => "Intron",
		"ncRNA_exonic" => "RNA",
		"ncRNA_intronic" => "RNA",
		"ncRNA_splicing" => "RNA",
		"ncRNA_UTR5" => "RNA",
		"ncRNA_UTR3" => "RNA",
		"ncRNA" => "RNA",
		"exonic" => "RNA"
);

if ($file =~ /gz$/){
	open IN, "gzip -dc $file |" or die $!;
}else{
	open IN, "$file" or die $!;
}

while(<IN>){
	chomp;
	next if (/^Chr/);
	my @tmp = split /\t/;
	next if ($tmp[4] eq "0");
	for (my $i=0;$i<@tmp;$i++){
		if($tmp[$i] eq "."){
			$tmp[$i] = "NA";
		}
	}
	my $hugo_symbol = $tmp[6];
	my $chr = $tmp[0];
	my $start = $tmp[1];
	my $end = $tmp[2];
	my $ref_allele = $tmp[3];
	my $var_allele1 = $tmp[3];
	my $var_allele2 = $tmp[4];
	my $genedetail = $tmp[7];
	my $rmsk = $tmp[10];
	my $genehancer = $tmp[11];
	my $cosmic = $tmp[12];
	my $mcap = $tmp[60];
	my $revel = $tmp[63];
	my $cadd = $tmp[67];
	my $Polyphen2_HDIV = $tmp[31];
	my $Polyphen2_HVAR = $tmp[34];
	my $FATHMM = $tmp[46];
	my $clinid = $tmp[15];
	my $clindn = $tmp[16];
	my $clinsig = $tmp[19];
	my $Interpro_domain = $tmp[95];
	my $GTEx_V6p_tissue = $tmp[97];
	my $ExAC_ALL = $tmp[20]; 
	my $ExAC_EAS = $tmp[23];
	my $ExAC_SAS = $tmp[27];
	my $gnomad312_AF_all = $tmp[98];
	my $gnomad312_AF_eas = $tmp[108];
	my $filter = $tmp[123];
	my $snp_rs = $tmp[119];
	my $AF_all = $tmp[114];
	my $common = 0;
	if($tmp[124] =~ /COMMON/){
		$common = 1;
	}
	if($tmp[124] =~ /RS=(\S+?);/){
		$snp_rs = "rs$1";
	}
	my $ref = $tmp[120];
	my $alt = $tmp[121];
	my $variant_class;
	my ($aa_change, $tr_id, $tx_change, $exon_num) = ("NA", "NA", "NA", "NA");
	my @func  = split(/;/, $tmp[5]);
	if(@func > 1){
		my @gene = split(/;/, $tmp[6]);
		$hugo_symbol = $gene[0];
	}
	if ($func[0] eq "exonic"){
		my @exonic = split(/;/, $tmp[8]);
		$variant_class = $variant{$exonic[0]};
		$genedetail = $tmp[9];
		$tmp[9] = ":NA:NA:NA:NA" if ($tmp[9] eq "UNKNOWN");
		my @change = split(/:/, $tmp[9]);
		if ($tmp[9] =~ /whole/){
			push @change, "NA";
			push @change, "NA";
		}
		$aa_change = $change[-1];
		$tr_id = $change[-4];
		$tx_change = $change[-2];
		$exon_num = $change[-3];
	}else{
		$variant_class = $variant{$func[0]};
	}
	my $variant_type;
	if ($tmp[3] eq "-"){
		$variant_type = "INS";
	}elsif($tmp[4] eq "-"){
		$variant_type = "DEL";
	}else{
		$variant_type = "SNP";
	}
	my @sam_id = sort keys %sam;
	for (my $i=126; $i<@tmp; $i++){
		my $j = $i-126;
		my $gt = $tmp[$i];
		my $center = $sam{$sam_id[$j]};
		my $samid = $sam_id[$j];
		my @tmp2 = split (/:/, $tmp[$i]);
		my @format = split (/:/, $tmp[125]);
		my %hash;
		for (my $g = 0; $g< @format; $g++){
			$hash{$format[$g]} = $tmp2[$g]
		}
		#		next if ($tmp2[0] eq "0/0" ||  $tmp2[0] eq "\./\." || $tmp2[0] eq "\.|\." || $tmp2[0] eq "0|0" || $tmp2[0] eq "\.");
		my ($ref_r, $alt_r, $vaf, $dep);
		if($type eq "Mutect2"){
			$vaf = $hash{'AF'};
			($ref_r,$alt_r) = split(/,/, $hash{'AD'});
			$dep = $hash{'DP'};
		}elsif($type eq "HapCaller"){
			($ref_r,$alt_r) = split(/,/, $hash{'AD'});
			$dep = $hash{'DP'};
			next if ($dep == 0);
			$vaf = sprintf("%.6f", $alt_r / $dep);
		}elsif($type eq "LoFreq"){
			$vaf = $hash{'AF'};
			$dep = $hash{'DP'};
			my @DP4 = split(/,/,$hash{'DP4'});
			$ref_r = $DP4[0] + $DP4[1];
			$alt_r = $DP4[2] + $DP4[3];
		}elsif($type eq "Pindel"){
			($ref_r,$alt_r) = split(/,/, $hash{'AD'});
			$dep = $ref_r + $alt_r;
			next if ($dep == 0);
			$vaf = sprintf("%.6f", $alt_r / $dep);
		}elsif($type eq "VarDict"){
			$dep = $hash{'DP'};
			$vaf = $hash{'AF'};
			($ref_r,$alt_r) = split(/,/, $hash{'AD'});
		}elsif($type eq "Pisces"){
			$dep = $hash{'DP'};
			($ref_r,$alt_r) = split(/,/, $hash{'AD'});
			$vaf = $hash{'VF'};
		}else{
			die $!;
		}
		$var_allele1 = $tmp[4] if($vaf > 0.99);
		print "$hugo_symbol\t$entrez_gene_id\t$center\t$ncbi_build\t$chr\t$start\t$end\t$strand\t$variant_class\t$variant_type\t$ref_allele\t$var_allele1\t$var_allele2\t";
		print "$snp_rs\t$samid\t$mutation_status\t$aa_change\t$tr_id\t$tx_change\t$exon_num\t$genedetail\t";
		print "$rmsk\t$genehancer\t$cosmic\t$clinid\t$clindn\t$clinsig\t$Interpro_domain\t$GTEx_V6p_tissue\t";
		print "$mcap\t$revel\t$cadd\t$Polyphen2_HDIV\t$Polyphen2_HVAR\t$FATHMM\t";
		print "$ExAC_ALL\t$ExAC_EAS\t$ExAC_SAS\t$gnomad312_AF_all\t$gnomad312_AF_eas\t";
		print "$AF_all\t$dep\t$filter\t$common\t$ref\t$alt\t$gt\t$ref_r\t$alt_r\t$vaf\n";
	}
}
close IN;
