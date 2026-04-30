#!/usr/bin/perl

use strict;
use warnings;
use File::Temp qw(tempfile);
use Getopt::Long;

my ($tmpdir, $help);
GetOptions(
	"tmpdir:s" => \$tmpdir,
	"h"        => \$help,
);

my $file = shift;
my $usage = <<"USAGE";
Usage: perl $0 [--tmpdir DIR] <input.maf>

Sort a MAF/TSV file by chromosome, position, allele, sample, and caller.
This version writes sortable keys to a temporary file and delegates sorting to
the system sort command, so memory use stays low for large MAF files.
USAGE

die $usage if ($help || !$file);
$tmpdir = "." unless defined $tmpdir;

open my $IN, "<", $file or die "Cannot open $file: $!\n";
my $header = <$IN>;
die "Empty input file: $file\n" unless defined $header;
chomp($header);
print "$header\n";

my @header = split(/\t/, $header, -1);
my %col_index = map { $header[$_] => $_ } 0..$#header;
my @required_columns = qw(
	Chromosome
	Start_Position
	End_Position
	Reference_Allele
	Tumor_Seq_Allele1
	Tumor_Seq_Allele2
	Tumor_Sample_Barcode
);

for my $column (@required_columns) {
	die "Missing required MAF column: $column\n" unless exists $col_index{$column};
}

my $mutation_status_index = $col_index{"Mutation_Status"};
my @temp_args = (UNLINK => 1, DIR => $tmpdir);
my ($TMP, $temp_file) = tempfile("maf_sort_by_pos_v2.XXXXXX", @temp_args);

while (my $line = <$IN>) {
	chomp($line);
	next if $line =~ /^Hugo_Symbol/;
	my @row = split(/\t/, $line, -1);
	my $chromosome = value_at(\@row, \%col_index, "Chromosome");
	my $chromosome_rank = chromosome_rank($chromosome);
	my $start = numeric_value(value_at(\@row, \%col_index, "Start_Position"));
	my $end = numeric_value(value_at(\@row, \%col_index, "End_Position"));
	my $caller = defined $mutation_status_index ? ($row[$mutation_status_index] // "") : "";
	print $TMP join(
		"\t",
		$chromosome_rank,
		$chromosome,
		$start,
		$end,
		value_at(\@row, \%col_index, "Reference_Allele"),
		value_at(\@row, \%col_index, "Tumor_Seq_Allele1"),
		value_at(\@row, \%col_index, "Tumor_Seq_Allele2"),
		value_at(\@row, \%col_index, "Tumor_Sample_Barcode"),
		$caller,
		$line
	), "\n";
}
close $IN;
close $TMP or die "Cannot close temporary file $temp_file: $!\n";

my @sort_cmd = (
	"sort",
	"-t", "\t",
	"-k1,1n",
	"-k2,2",
	"-k3,3n",
	"-k4,4n",
	"-k5,5",
	"-k6,6",
	"-k7,7",
	"-k8,8",
	"-k9,9",
	$temp_file,
);

open my $SORTED, "-|", @sort_cmd or die "Cannot run sort: $!\n";
while (my $line = <$SORTED>) {
	chomp($line);
	my @fields = split(/\t/, $line, 10);
	print "$fields[9]\n" if defined $fields[9];
}
close $SORTED or die "sort failed for $file\n";

sub value_at {
	my ($row, $col_index, $column) = @_;
	my $index = $col_index->{$column};
	return defined $row->[$index] ? $row->[$index] : "";
}

sub numeric_value {
	my ($value) = @_;
	return $value =~ /^-?\d+(?:\.\d+)?$/ ? $value : 0;
}

sub chromosome_rank {
	my ($chromosome) = @_;
	$chromosome = "" unless defined $chromosome;
	my $normalized = uc($chromosome);
	$normalized =~ s/^CHR//;
	return sprintf("%03d", $normalized) if $normalized =~ /^\d+$/;
	return "023" if $normalized eq "X";
	return "024" if $normalized eq "Y";
	return "025" if $normalized eq "M" || $normalized eq "MT";
	return "999";
}
