#!/usr/bin/perl

use strict;
use warnings;

my $inter_f = shift;
my $intra_f = shift;
my $min_con_read = shift;
my $ibthr_f = shift;        
my $outbid_f = shift;

my $gMaxScore = 1.0;
my $gMinScore = 0.0;
my $gDiff = $gMaxScore - $gMinScore;

my %hs_pct_thrs = ();
my %hs_raw_thrs = ();
my ($pct_sd, $raw_sd) = (0, 0);
my ($pct_50, $raw_50) = (0, 0);
my ($pct_2_5, $raw_2_5) = (0, 0);
my ($pct_16, $raw_16) = (0, 0);
my ($pct_84, $raw_84) = (0, 0);
my ($pct_97_5, $raw_97_5) = (0, 0);
open(F,"$ibthr_f");
while(<F>) {
	chomp;
    if ($_ =~ /^#/) {
        my @ar = split(/\s+/);
        ($pct_sd, $raw_sd) = ($ar[2], $ar[4]);
    } else {
        my ($perc, $pct_thr, $raw_thr) = split(/\s+/);
        $hs_pct_thrs{$pct_thr} = $perc;
        $hs_raw_thrs{$raw_thr} = $perc;
        
        if ($perc == 50) {
            $pct_50 = $pct_thr;
            $raw_50 = $raw_thr;
        } elsif ($perc == 2.5) {
            $pct_2_5 = $pct_thr;
            $raw_2_5 = $raw_thr;
        } elsif ($perc == 16) {
            $pct_16 = $pct_thr;
            $raw_16 = $raw_thr;
        } elsif ($perc == 84) {
            $pct_84 = $pct_thr;
            $raw_84 = $raw_thr;
        } elsif ($perc == 97.5) {
            $pct_97_5 = $pct_thr;
            $raw_97_5 = $raw_thr;
        } 
    }
}
close(F);

my @pct_thrs = sort {$a<=>$b} keys %hs_pct_thrs;
my @raw_thrs = sort {$a<=>$b} keys %hs_raw_thrs;

open(OB, ">$outbid_f");
print OB "#raw_2.5: $raw_2_5\traw_16: $raw_16\traw_50: $raw_50\traw_84: $raw_84\traw_97.5: $raw_97_5\traw_sd: $raw_sd\n";
my $diff = 0;

if($intra_f ne "NA"){
	open(F,"$intra_f");
	while(<F>) {
	    chomp;
		my @ar = split(/\s+/);
	    my ($block1, $dir1, $block2, $dir2, $avgscore) = split(/\s+/);
		my ($bid1, $bid2) = ($ar[$#ar-1], $ar[$#ar]);
	    $ar[7] =~ /MIN=(\S+)/;
	    my $score = $1;
	
	    if ($score < $min_con_read) { next; }

    	my $newscore = convert_score($score, $raw_50, \@raw_thrs, \%hs_raw_thrs);
		printf OB "$bid1 $bid2\t%.10f\t$score\n", $newscore; 
	}
	close(F);
}

open(F, "$inter_f");
while(<F>) {
	chomp;
	my ($block1, $dir1, $block2, $dir2, $count, $bid1, $bid2) = split(/\s+/);
	if ($count < $min_con_read) { next; }
	
    my $newscore = convert_score($count, $raw_50, \@raw_thrs, \%hs_raw_thrs);
    printf OB "$bid1 $bid2\t%.10f\t$count\n", $newscore; 
}
close(F);
close(OB);

##########
sub convert_score {
	my $org_score = shift;
    my $avg_thrs = shift;
	my $rar_thrs = shift;
	my $rhs_thrs = shift;

	my $i;
	my $perc;
	for ($i = 0; $i < scalar(@$rar_thrs); $i++) {
		my $thr = $$rar_thrs[$i];
		$perc = $$rhs_thrs{$thr};
		if ($thr > $org_score) { last; } ##
	}
		
	my $new_score = -1;
    if ($perc > 97.5) { # Over two sds from mean, for handling repetitive regions
        $new_score = 0;
    } elsif ($org_score >= $avg_thrs) {
		$new_score = 1;
	} else {
        $new_score = $org_score/$avg_thrs;
	} 

	return $new_score;
}
