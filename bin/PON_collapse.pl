#!/usr/bin/perl

use warnings;
use strict;

my $maf = shift or die $!;


my %hash;
open IN, "$maf" or die $!;
my @header = split(/\t/,<IN>);
print join("\t",@header[0..39]) . "\tREF_R\tALT_R\tNS\tNS2\tDepth\n";
while(<IN>){
	chomp;
	my @tmp = split /\t/;
	$tmp[14] = "PON_combind"; 
	$tmp[15] = "PON16";
	my $mut = join("\t",@tmp[4,5,6,10,12]);
	$hash{$mut}{"INFO"} = join("\t",@tmp[0..39]);
	$hash{$mut}{"NS"} += 1;
	if($tmp[50] >= 0.02){
		$hash{$mut}{"NS2"} += 1;
	}
	my @ref_c = split(/[|,]/, $tmp[47]);
	my @alt_c = split(/[|,]/, $tmp[48]);
	@ref_c = sort {$b <=> $a} @ref_c;
	@alt_c = sort {$b <=> $a} @alt_c;
	$hash{$mut}{"REF_R"} += $ref_c[0];
	$hash{$mut}{"ALT_R"} += $alt_c[0];
}
close IN;

for my $i (sort {$a cmp $b} keys %hash){
	$hash{$i}{"NS2"} = 0 unless(defined $hash{$i}{"NS2"});
	my @mut = split("\t", $i);
	my $i16 = `bcftools mpileup -f /lustre/home/acct-medkkw/medlyb/database/annotation/gatk_ann/hg38/bwaindex2/Homo_sapiens_assembly38.fasta -b PON_bam.list -S PON_sample.list -r $mut[0]:$mut[1] -O v  -d 100000 -Q 5 -q 5 | bcftools query -f '%I16'`;
	my @dep = split(/,/, $i16);
	my $depth = join(",", @dep[0..3]);
	print "$hash{$i}{'INFO'}\t$hash{$i}{'REF_R'}\t$hash{$i}{'ALT_R'}\t$hash{$i}{'NS'}\t$hash{$i}{'NS2'}\t$depth\n";
}

