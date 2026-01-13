#!/usr/bin/perl

use strict;
use warnings;

my $file = shift or die $!;

my %hash;
open IN, "$file" or die $!;
my $head = <IN>;
print $head;
while(<IN>){
	chomp;
	my @tmp = split /\t/;
	$hash{$tmp[4]}{$tmp[5]}{$tmp[6]}{$tmp[10]}{$tmp[12]}{$tmp[14]}{$tmp[15]} = $_;
}
close IN;

my @chr;
open IN, "/lustre/home/acct-medkkw/medlyb/chr.list" or die $!;
while(<IN>){
	chomp;
	push @chr, $_;
}
close IN;

for my $c (@chr){
	for my $s (sort {$a <=> $b} keys %{$hash{$c}}){
		for my $e (sort {$a <=> $b} keys %{$hash{$c}{$s}}){
			for my $r (sort {$a cmp $b} keys %{$hash{$c}{$s}{$e}}){
				for my $v (sort {$a cmp $b} keys %{$hash{$c}{$s}{$e}{$r}}){
					for my $i (sort {$a cmp $b} keys %{$hash{$c}{$s}{$e}{$r}{$v}}){
						for my $m (sort {$a cmp $b} keys %{$hash{$c}{$s}{$e}{$r}{$v}{$i}}){
							 print "$hash{$c}{$s}{$e}{$r}{$v}{$i}{$m}\n";
						}
					}
				}
			}
		}
	}
}

close IN;
