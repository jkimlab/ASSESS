#!/usr/bin/perl

use strict;
use warnings;

my $blist_f = shift;
my $fjoinbreak_f = shift;

my %hs_crds = ();
my %hs_bids = ();
my %hs_dirs = ();
open(F, $blist_f);
while(<F>) {
	chomp;
	my ($scf, $start, $end, $dir, $absbid) = split(/\s+/);
	$hs_crds{$scf}{$start} = $end;
	my $bid = $absbid;
	if ($dir eq "-") { $bid = -1*$absbid; }
	$hs_bids{$scf}{$start} = $bid;
	$hs_dirs{$absbid} = $dir;
}
close(F);

open(F, $fjoinbreak_f);
while(<F>) {
	chomp;
	if (length($_) == 0 || $_ =~ /^#/) { next; }
	if ($_ =~ /^SF\s+(\S+)\s+(\S+)\s+(\S+)/) {
		my ($bid1, $bid2, $label) = ($1, $2, uc($3));
		if ($label eq "KEEP") { $label = "JOIN"; }

		print "$bid1\t$bid2\t$label\n";
	} else {
		my ($crd, $label) = split(/\s+/);
		$label = uc($label);
		if ($label eq "KEEP") { $label = "JOIN"; }

		if ($crd =~ /(\S+):(\S+)\-(\S+)/) {
			# format: scf13:3451-5000
			my ($scf, $tstart, $tend) = ($1, $2, $3);
			my $rhs = $hs_crds{$scf};
			if (!defined($rhs)) { next; }
	
			my @starts = sort {$a<=>$b} keys %$rhs; 
			if (scalar(@starts) == 1) { next; }

			my %hs_break_crds = ();
			my %hs_break_bids = ();
			for (my $i = 0; $i < $#starts; $i++) {
				my $start1 = $starts[$i];
				my $end1 = $$rhs{$start1};
				my $start2 = $starts[$i+1];
				my $end2 = $$rhs{$start2};

				my $bid1 = $hs_bids{$scf}{$start1};
				my $bid2 = $hs_bids{$scf}{$start2};

				my ($bstart, $bend) = ($end1, $start2);
				if ($bend < $bstart) { ($bstart, $bend) = ($start2, $end1); }
			
				$hs_break_crds{$bstart} = $bend;
				$hs_break_bids{$bstart} = "$bid1\t$bid2";
			}

			foreach my $start (sort {$a<=>$b} keys %hs_break_crds) {
				my $end = $hs_break_crds{$start};
				if ($tstart <= $end && $start <= $tend) {
					my $bidpair = $hs_break_bids{$start};
					my ($bid1, $bid2) = split(/\t/, $bidpair);
					print "$bid1\t$bid2\t$label\n";
				}
			}
		} else {
			# format: scf13	
			my $rhs = $hs_crds{$crd};
			if (!defined($rhs)) { next; }
			
			my @bids = ();
			foreach my $start (sort {$a<=>$b} keys %$rhs) {
				my $bid = $hs_bids{$crd}{$start};
				push(@bids, $bid);
			}

			if (scalar(@bids) == 1) { next; }
		
			for (my $i = 0; $i < $#bids; $i++) {
				my $bid1 = $bids[$i];
				my $bid2 = $bids[$i+1];

				print "$bid1\t$bid2\t$label\n"; 
			}
		}
	}
}
close(F);
