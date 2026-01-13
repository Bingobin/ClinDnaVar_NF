#!/usr/bin/perl

use strict;
use warnings;

my $vcf = shift or die $!;
my $sample = shift or die $!;

open IN, "$vcf" or die $!;
while(<IN>){
	chomp;
	if(/^##/){
		print "$_\n";
	}
}
close IN;

print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
print "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Raw Depth\">\n";
print "##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">\n";
print "##FORMAT=<ID=SB,Number=1,Type=Integer,Description=\"Phred-scaled strand bias at this position\">\n";
print "##FORMAT=<ID=DP4,Number=4,Type=Integer,Description=\"Counts for ref-forward bases, ref-reverse, alt-forward and alt-reverse bases\">\n";
print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sample\n";

open IN, "$vcf" or die $!;
while(<IN>){
	chomp;
	next if (/^#/);
	my @tmp = split /\t/;
	my %hash;
	$hash{'GT'} = "0/1";
	my @info = split (/;/,$tmp[7]);
	for my $i (@info){
		my @e = split(/=/, $i);
		$hash{$e[0]} = $e[1];
	}
	print "$_\tGT:DP:AF:SB:DP4\t$hash{'GT'}:$hash{'DP'}:$hash{'AF'}:$hash{'SB'}:$hash{'DP4'}\n";
}
close IN;
