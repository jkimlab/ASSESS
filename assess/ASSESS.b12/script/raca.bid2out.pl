#!/usr/bin/perl

use strict;
use warnings;

my $oraca_f = shift;	# raca.bid.out 
my $blist_f = shift;	# block_list
my $score_f = shift;    
my $scf_score_f = shift;

my %hs_bids = ();
my %hs_dir = ();
open(F, "$blist_f");
while(<F>) {
	chomp;
	my ($scf, $start, $end, $dir, $bid) = split(/\s+/);
	my $crd = "$scf:$start-$end";
	$hs_bids{$bid} = $crd;
	$hs_dir{$bid} = $dir;
}
close(F);

my %hs_scores = ();
open(F, $score_f);
while(<F>) {
    chomp;
    my ($bid1, $bid2, $score) = split(/\s+/);
    my ($rbid1, $rbid2) = (-1*$bid1, -1*$bid2);
    $hs_scores{"$bid1:$bid2"} = $score;
    $hs_scores{"$rbid2:$rbid1"} = $score;
}
close(F);

if($scf_score_f ne "NA"){
	open(F, $scf_score_f);
	while(<F>) {
	    chomp;
	    my ($bid1, $bid2, $score) = split(/\s+/);
	    my ($rbid1, $rbid2) = (-1*$bid1, -1*$bid2);
	    $hs_scores{"$bid1:$bid2"} = $score;
	    $hs_scores{"$rbid2:$rbid1"} = $score;
	}
	close(F);
}

open(F, $oraca_f);
while(<F>) {
    chomp;
    if (length($_) == 0) { next; }

    if ($_ =~ /^#(\S+)/) {
        my $chrid = $1; 
        print ">$chrid\n";
    } else {
        my @ar = split(/\s+/);
        pop(@ar);
        for (my $i = 0; $i < $#ar; $i++) {
            my $bid1 = $ar[$i];
            my $bid2 = $ar[$i+1];

            my $score = $hs_scores{"$bid1:$bid2"};

            my $crd1 = $hs_bids{abs($bid1)};
            my $crd2 = $hs_bids{abs($bid2)};
            my $dir1 = $hs_dir{abs($bid1)};
            my $dir2 = $hs_dir{abs($bid2)};

            my $fdir1 = $dir1;
            if ($bid1 < 0) {
                if ($fdir1 eq "+") { $fdir1 = "-"; }
                else { $fdir1 = "+"; }
            }
            my $fdir2 = $dir2;
            if ($bid2 < 0) {
                if ($fdir2 eq "+") { $fdir2 = "-"; }
                else { $fdir2 = "+"; }
            }
			if (defined($score)) {
            	print "$crd1 $fdir1\t$crd2 $fdir2\t$score\n";
			} else {
            	print "$crd1 $fdir1\t$crd2 $fdir2\tNA\n";
			}
        }
        print "\n";
    }
}
close(F);

