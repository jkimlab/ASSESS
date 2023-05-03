#!/usr/bin/perl

use strict;
use warnings;

my $raca_dir = shift;
my $tarspc = shift;
my $gap_f = shift;

my $block_list_f = "$raca_dir/SFs/block_list.txt";
my %hs_block = ();
open(F, $block_list_f);
while(<F>) {
	chomp;
	my @cols = split(/\s+/);
	$hs_block{$cols[4]} = "$cols[0]:$cols[1]-$cols[2]";
}	
close(F);

my %gap_info = ();
open(F, $gap_f);
while(<F>) {
	my ($scf, $start, $end, $size) = split(/\s+/);
	$gap_info{$scf}{$start} = $end;
}
close(F);


my $adjinfo_all_f = "$raca_dir/adjinfo.all.txt";
my %label = ();
my %intra_score = ();
open(F, $adjinfo_all_f);
<F>;
while(<F>) {
	chomp;
	my @cols = split(/\s+/);
	$cols[0] =~ /(\S+)\(\S+\):(\S+)\(\S+\)/;
	my ($bid1, $bid2) = ($1, $2);
	my $label = "$cols[-3]\t$cols[-2]";
	$label{$bid1}{$bid2} = $label;
	$label{$bid2*(-1)}{$bid1*(-1)} = $label;
	if ($_ =~ /intra:(\S+)/) {
		$intra_score{$bid1}{$bid2} = $1;
		$intra_score{$bid2*(-1)}{$bid1*(-1)} = $1;
	}
}	
close(F);

my $raca_bid_out_f = "$raca_dir/raca2.flt.bid.out";
#my $raca_bid_out_f = "$raca_dir/raca2.flt.ext.bid.out";
my $raca_frag_n = "";
my %frag_bid_order = ();
open(F, $raca_bid_out_f);
while(<F>) {
	chomp;
	if ($_ =~ /^#(\S+)/) { 
		$raca_frag_n = $1;
		next;
	}
	my @bids = split(/\s+/, substr($_,0,-1));
	for (my $i=0; $i<$#bids; $i++) {
		my ($bid1, $bid2) = ($bids[$i], $bids[$i+1]);
		$frag_bid_order{$bid1}{$bid2} = "";
		$frag_bid_order{$bid2*(-1)}{$bid1*(-1)} = "";
	}
}
close(F);

print"bid1:bid2\tintra_score\tn_count\tspc_label\tr_label\n";
my $g_order_f = "$raca_dir/SFs/Genomes.Order";
my $check = 0;
my %g_order = ();
my $scf_n = "";
open(F, $g_order_f);
while(<F>) {
	chomp;
	if ($_ =~ />(\S+)/) {
		$check = ($tarspc eq $1)? 1 : 0;
	} else {
		if ($check == 0) { next ;} 
		if ($_ =~ /^# (\S+)/) {
			$scf_n = $1;
		} else {
			my @bids = split(/\s+/, substr($_,0,-1));
			if (scalar(@bids) == 1) { next; }
			for (my $i=0; $i<$#bids; $i++) {
				my ($bid1, $bid2) = ($bids[$i], $bids[$i+1]);
				if (! exists $frag_bid_order{$bid1}{$bid2}) {
					my @bid1_coords = split(/[:|-]/, $hs_block{abs($bid1)});
					my @bid2_coords = split(/[:|-]/, $hs_block{abs($bid2)});
					my ($bid1_e, $bid2_s) = ($bid1_coords[2], $bid2_coords[1]+1);
					my $N_num = 0;
					foreach my $gap_s (sort {$a<=>$b} keys %{$gap_info{$scf_n}}) {
						my $gap_e = $gap_info{$scf_n}{$gap_s};
						my $gap_size = 0;
						if ($gap_s >= $bid1_e && $gap_e <= $bid2_s) {
							$gap_size = $gap_e - $gap_s + 1;
							$N_num += $gap_size;
						} elsif ($gap_s >= $bid1_e && $gap_s <= $bid2_s) {
							$gap_size = $bid2_s - $gap_s + 1;
							$N_num += $gap_size;
						} elsif ($gap_e >= $bid1_e && $gap_e <= $bid2_s) {
							$gap_size = $gap_e - $bid1_e + 1;
							$N_num += $gap_size;
						} 
					}
					print "$bid1:$bid2\t$intra_score{$bid1}{$bid2}\t$N_num\t$label{$bid1}{$bid2}\n";
				}
			}
		}
	}
}
close(F);

