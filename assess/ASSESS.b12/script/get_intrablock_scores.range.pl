#!/usr/bin/perl

use strict;
use warnings;

my $sstart = shift;
my $send = shift;
my $block_f = shift;
my $flanksize = shift;
my $windowsize = shift;
my $cov_dir = shift;
my $avgcov_f = shift;
my $out_f = shift;

my $avgcov_cutoff = 2;	# minimum coverage cutoff

my $stepsize = int($windowsize/2);

my %hs_blocks = ();
my %hs_blockschr = ();
open(F,"$block_f");
while(<F>) {
	chomp;
	my @ar = split(/\s+/);
	my ($scf, $start, $end, $rchr) = ($ar[0], $ar[1], $ar[2], $ar[5]);
	$hs_blocks{$scf}{$start} = $end;
	$hs_blockschr{$scf}{$start}{$end} = $rchr;
}
close(F);

my %hs_avgcov = ();
open(F, $avgcov_f);
while(<F>) {
	chomp;
	my ($scf, $avgcov) = split(/\s+/);
	$hs_avgcov{$scf} = $avgcov;
}
close(F);

my @finals = ();

my @scfs = sort keys %hs_blocks;
for (my $i = $sstart; $i <= $send; $i++) {
	my $scf = $scfs[$i];
	my $rhs = $hs_blocks{$scf};
	my @starts = sort {$a<=>$b} keys %$rhs;
	if (scalar(@starts) == 1) { next; }

	# read coverage data
	my $f = "$cov_dir/$scf.cov";
	if (!(-f $f)) {
		print STDERR "\nWarning: File $f not exist.\n"; next; 
	}

	my $gavg = $hs_avgcov{$scf};
    if ($gavg == 0) {
		print STDERR "\nWarning: No coverage data in file $f.\n"; next; 
    }

	my %hs_cov = ();
	open(F,"$f");
	while(<F>) {
		chomp;
		my ($pos, $cov) = split(/\s+/);
		$hs_cov{$pos} = $cov;
	}
	close(F);

	my %hs_lowcov_regions = ();

	for (my $i = 0; $i < $#starts; $i++) {
		my $start1 = $starts[$i];
		my $end1 = $$rhs{$start1};
		my $start2 = $starts[$i+1];
		my $end2 = $$rhs{$start2};

		my $rchr1 = $hs_blockschr{$scf}{$start1}{$end1};
		my $rchr2 = $hs_blockschr{$scf}{$start2}{$end2};

		my ($brstart, $brend) = ($end1, $start2);
		if ($end1 > $start2) { ($brstart, $brend) = ($start2, $end1); }

		$brstart -= $windowsize * 5;
		my $mid1 = $end1 - int(($end1 - $start1)/3);
		if ($brstart < $mid1) { $brstart = $mid1; }

		$brend += $windowsize * 5;
		my $mid2 = $start2 + int(($end2 - $start2)/3);
		if ($brend > $mid2) { $brend = $mid2; }

		my $minavg = -1;
		my ($min_wstart, $min_wend) = (0,0);
		for (my $p = $brstart; $p < $brend; $p += $stepsize) {
			my ($wstart, $wend) = ($p, $p+$windowsize);
			if ($wend > $brend) { $wend = $brend; $p = $brend; }

			my $sum = 0;
			my $cnt = 0;
			for (my $pp = $wstart+1; $pp<=$wend; $pp++) {
				my $cov = $hs_cov{$pp};
				$sum += $cov;
				$cnt++;
			}
			
            my $avgcov = $sum/$cnt;
			if ($minavg == -1 || $avgcov < $minavg) {
				$minavg = $avgcov;
				($min_wstart, $min_wend) = ($wstart, $wend);
			}

			if ($avgcov < $avgcov_cutoff) {
				$hs_lowcov_regions{$wstart}{$wend} = $avgcov;
			}
		}
		my $type = "INTER";
		if ($rchr1 eq $rchr2) { $type = "INTRA"; }
		my $minpct = $minavg/$gavg;
		my $out = "$scf:$start1-$end1 +\t$scf:$start2-$end2 +\t$minpct\t$brstart-$brend\t$min_wstart-$min_wend\tMIN=$minavg\t$gavg\t$type\t$rchr1 $rchr2\t";
		
		my $pavgcov = -1;
		my ($prstart, $prend) = (-1, -1);
		foreach my $rstart (sort {$a<=>$b} keys %hs_lowcov_regions) {
			my $rhs = $hs_lowcov_regions{$rstart};
			foreach my $rend (sort {$a<=>$b} keys %$rhs) {
				my $avgcov = $$rhs{$rend};
				if ($pavgcov == -1) {
					($prstart, $prend) = ($rstart, $rend);
					$pavgcov = $avgcov;
				} elsif ($avgcov == $pavgcov) {
					if ($rstart < $prend) {
						$prend = $rend;
					} else {
						$out .= "($prstart $prend $pavgcov) ";
						($prstart, $prend) = ($rstart, $rend);
						$pavgcov = $avgcov;
					}
				} else {
					$out .= "($prstart $prend $pavgcov) ";
					($prstart, $prend) = ($rstart, $rend);
					$pavgcov = $avgcov;
				}
			}
		}
		if ($prstart != -1) {
			$out .= "($prstart $prend $pavgcov) ";
		}
		push(@finals, $out);
	} 	
}

open(O,">$out_f");
foreach my $out (@finals) {
	print O "$out\n";
}
close(O);
