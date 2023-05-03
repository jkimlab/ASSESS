#!/usr/bin/perl

use strict;
use warnings;
use FindBin qw($Bin);

my $sstart = shift;
my $send  = shift;
my $block_f = shift;
my $flanksize = shift;
my $windowsize = shift;
my $cov_dir = shift;
my $out_f = shift; 

my $stepsize = int($windowsize/2);

my %hs_blocks = ();
open(F,"$block_f");
while(<F>) {
	chomp;
	my ($scf, $start, $end) = split(/\s+/);
	$hs_blocks{$scf}{$start} = $end;
}
close(F);

my @pct_values = ();
my @raw_values = ();
my %hs_avgcov = ();

my @scfs = sort keys %hs_blocks;
for (my $i = $sstart; $i <= $send; $i++) {
	my $scf = $scfs[$i];
	my $rhs = $hs_blocks{$scf};
	my @starts = sort {$a<=>$b} keys %$rhs;

	# read coverage data
	my $f = "$cov_dir/$scf.cov";
	if (!(-f $f)) {
		print STDERR "File $f doesn't exist...\n"; next; 
	}

	my %hs_cov = ();
	my $scfsize = 0;
	my $gsum = 0;
	my $gcnt = 0;
	open(F,"$f");
	while(<F>) {
		chomp;
		my ($pos, $cov) = split(/\s+/);
		$hs_cov{$pos} = $cov;
		$scfsize = $pos;
		if ($cov != 0) { 
			# ignore coverage in potential gap regions
			$gsum += $cov;
			$gcnt++;
		}
	}
	close(F);

	# skip small scaffolds
	if ($scfsize <= 2*$flanksize) { next; }

    my $gavg = 0;
    if ($gcnt > 0) { $gavg = $gsum/$gcnt; }
	$hs_avgcov{$scf} = $gavg;

    if ($gcnt == 0) { next; }

	for (my $p = $flanksize+1; $p < $scfsize - $flanksize; $p += $stepsize) {
		my ($wstart, $wend) = ($p, $p+$windowsize-1);
		if ($wend >= $scfsize-$flanksize) { $wend = $scfsize-$flanksize-1; $p = $scfsize-$flanksize; }
		
		my $sum = 0;
		my $cnt = 0;
		for (my $pp = $wstart; $pp <= $wend; $pp++) {
			my $cov = $hs_cov{$pp};
			$sum += $cov;
			$cnt++;
		}
		
		if ($sum == 0) {
			# ignore coverage in potential gap regions
			next;
		}
		
		my $avg = $sum/$cnt;
		my $pct = $avg/$gavg;
		push(@pct_values, $pct);
		push(@raw_values, $avg);
	}
}

# Output to a file
open(O,">$out_f");
for (my $i = 0; $i <= $#pct_values; $i++) {
	my $pct = $pct_values[$i];
	my $raw = $raw_values[$i];
	print O "$pct\t$raw\n";
}

foreach my $scf (sort keys %hs_avgcov) {
	my $avgcov = $hs_avgcov{$scf};
	print O "#\t$scf\t$avgcov\n";
}
close(O);
