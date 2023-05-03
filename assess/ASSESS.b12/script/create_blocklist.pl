#!/usr/bin/perl

use strict;
use warnings;

my $rspc = shift;
my $spc = shift;
my $fasize_f = shift;
my $data_dir = shift;

my %data = ();
my %data_blockid = ();
my %data_dir = ();
my %data_rchr = ();
my $f = "$data_dir/Conserved.Segments";
my $o = "$data_dir/block_list.txt";
my $oext = "$data_dir/block_list.ext.txt";

my %hs_scfsize = ();
open(F, "$fasize_f");
while(<F>) {
	chomp;
	my ($scf, $size) = split(/\s+/);
	$hs_scfsize{$scf} = $size;
}
close(F);

my $blockid = -1;
my $rchr = "";
open(F,"$f");
while(<F>) {
	chomp;

	if ($_ =~ /^>/) {
		$blockid = substr($_,1);
	} elsif ($_ =~ /^$rspc\.(\S+):(\S+)\-(\S+) (\S+)/) {
		$rchr = $1;
	} elsif ($_ =~ /^$spc\.(\S+):(\S+)\-(\S+) (\S+)/) {
		my ($scf, $start, $end, $dir) = ($1,$2,$3,$4);
		$data{$scf}{$start} = $end;	
		$data_blockid{$scf}{$start} = $blockid;	
		$data_dir{$scf}{$start} = $dir;	
		$data_rchr{$scf}{$start} = $rchr;	
	}
}
close(F);

open(OE, ">$oext");
open(O,">$o");
foreach my $scf (sort keys %data) {
	my $rhs = $data{$scf};
	$scf =~ /\D+(\d+)/;
	my $scfsize = $hs_scfsize{$scf};

	my $slast = scalar(keys %$rhs);
	my $flag_fullscf = ($slast == 1) ? 1 : 0;
	
	my $scnt = 1;
	my @starts = sort {$a<=>$b} keys %$rhs;
	#foreach my $start (sort {$a<=>$b} keys %$rhs) {
	for (my $i = 0; $i <= $#starts; $i++) {
		my $start = $starts[$i];
		my $end = $$rhs{$start};

		my $blockid = $data_blockid{$scf}{$start};
		my $dir = $data_dir{$scf}{$start};
		my $rchr = $data_rchr{$scf}{$start};
		print O "$scf\t$start\t$end\t$dir\t$blockid\t$rchr\n";

		if ($flag_fullscf == 1) {
			print OE "$scf\t0\t$scfsize\t$dir\t$blockid\n";
		} else {
			# check overlap
			my $newend = $end;
			if ($scnt < $slast) {
				my $nextstart = $starts[$i+1];
				if ($nextstart < $end) {
					$newend = int(($end + $nextstart)/2);		
				}
			}
			
			my $newstart = $start;
			if ($scnt > 1) {
				my $prevend = $$rhs{$starts[$i-1]};
				if ($prevend > $start) {
					$newstart = int(($prevend + $start)/2);	
				}
			}

			if ($scnt == 1) {
				print OE "$scf\t0\t$newend\t$dir\t$blockid\n";
			} elsif ($scnt == $slast) {
				print OE "$scf\t$newstart\t$scfsize\t$dir\t$blockid\n";
			} else {
				print OE "$scf\t$newstart\t$newend\t$dir\t$blockid\n";
			}
		}

		$scnt++;
	} 
}
close(O);
close(OE);
