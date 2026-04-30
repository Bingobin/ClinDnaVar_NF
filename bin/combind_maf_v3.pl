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

open my $IN, "<", $maf or die $!;
my $header = <$IN>;
chomp($header);
print "$header\tMean_VAF\tN_Callers\n";
my @header = split(/\t/, $header);
my %col_index = map { $header[$_] => $_ } 0..$#header;
my @key_columns = qw(Chromosome Start_Position End_Position Reference_Allele Tumor_Seq_Allele1 Tumor_Seq_Allele2 Tumor_Sample_Barcode);
for my $column (@key_columns, "Variant_Classification") {
	die "Missing required MAF column: $column\n" unless exists $col_index{$column};
}

my %multi_value = map { $_ => 1 } qw(Mutation_Status DEPTH GT FILTER REF_R ALT_R VAF);
my ($current_key, $current_record);

while (my $line = <$IN>) {
	chomp($line);
	next if ($line =~ /^Hugo_Symbol/);
	my @tmp = split(/\t/, $line, -1);
	my $key = variant_key(\@tmp, \%col_index, \@key_columns);
	if (defined $current_key && $key ne $current_key) {
		flush_record($current_record, \@header, \%multi_value, $retain_mode);
		$current_record = undef;
	}
	$current_key = $key;
	merge_record(\@tmp, \@header, \%multi_value, \$current_record);
}
close $IN;

flush_record($current_record, \@header, \%multi_value, $retain_mode) if defined $current_record;

sub variant_key {
	my ($row, $col_index, $key_columns) = @_;
	return join("\t", map { defined $row->[$col_index->{$_}] ? $row->[$col_index->{$_}] : "" } @$key_columns);
}

sub merge_record {
	my ($row, $header, $multi_value, $record_ref) = @_;
	if (!defined $$record_ref) {
		$$record_ref = {};
	}
	my $record = $$record_ref;
	for (my $n = 0; $n < @$header; $n++) {
		my $column = $header->[$n];
		my $value = defined $row->[$n] ? $row->[$n] : "";
		if ($multi_value->{$column}) {
			push @{$record->{$column}}, $value;
		} else {
			$record->{$column} = $value;
		}
	}
}

sub flush_record {
	my ($record, $header, $multi_value, $retain_mode) = @_;
	return unless defined $record;
	#next if ($record->{"COMMON"} == 1);
	return if ($retain_mode eq "exonic" && $record->{"Variant_Classification"} !~ /(Splice|Mutation|Frame)/i);
	my $caller_num = scalar(@{$record->{'Mutation_Status'} || []});
	return if $caller_num == 0;
	#	next if ($caller_num < 2);
	for my $column (@$header) {
		if ($multi_value->{$column}) {
			print join("|", @{$record->{$column} || []}) . "\t";
		} else {
			print "$record->{$column}\t";
		}
	}
	my $merge_vaf = 0;
	for my $v (@{$record->{'VAF'} || []}) {
		$merge_vaf += $v;
	}
	my $mean_vaf = sprintf("%.6f", $merge_vaf / $caller_num);
	print "$mean_vaf\t$caller_num\n";
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
