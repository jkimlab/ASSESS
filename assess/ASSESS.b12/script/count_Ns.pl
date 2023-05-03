#!/usr/bin/perl

use strict;
use warnings;

my $seq_f = shift;
my $size_f = shift;

my %size = ();
open(F, $size_f);
while(<F>) {
	chomp;
	my ($scf, $size) = split(/\s+/);
	$size{$scf} = $size;
}
close(F);

my $scf_n = "";
my $pre_n = 0;
my $pre_scf = "";
my $base_pos = 0;
my %N_pos = ();
my ($n_start, $n_end) = (0,0);
open(F, $seq_f);
while(<F>) {
	chomp;
	if ($_ =~ /^>(\S+)/) {
		$scf_n = $1;
		$base_pos = 0;
		($n_start, $n_end) = (0,0);
	} else {
		my @bases = split(//, $_);
		foreach my $base (@bases) {
			$base_pos++;
			if ($base eq "N" || $base eq "n") {
				if ($pre_n == 0) { $n_start = $base_pos; }
				else { $n_end = $base_pos; }
				if ($n_end == $size{$scf_n}) {
					$N_pos{$scf_n}{$n_start} = $n_end;	
					print "$scf_n\t$n_start\t$n_end\t", $n_end-$n_start+1, "\n";
				}
				$pre_n = 1;
			} else {
				if ($pre_n == 1) { 
					$N_pos{$scf_n}{$n_start} = $n_end; 
					print "$scf_n\t$n_start\t$n_end\t", $n_end-$n_start+1, "\n";
				}
				$pre_n = 0;
			}
		}
	}
}
close(F);

=p
foreach my $scf (sort keys %N_pos) {
	foreach my $s (sort {$a<=>$b} %{$N_pos{$scf}}) {
		my $e = $N_pos{$scf}{$s};

	}
}
=cut
