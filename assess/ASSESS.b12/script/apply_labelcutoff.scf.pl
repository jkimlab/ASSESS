#!/usr/bin/perl

use strict;
use warnings;

my $type = shift;
my $labelcutoff_f = shift;
my $raca_dir = shift;

my $sr_linkoutbid_f = "$raca_dir/scf_block_linkscores.SR.flt.bid.txt";
my $linkr_linkoutbid_f = "$raca_dir/scf_block_linkscores.10xR.flt.bid.txt";
my $lr_linkoutbid_f = "$raca_dir/scf_block_linkscores.LR.flt.bid.txt";

my %label_cutoff = ();
open(F, $labelcutoff_f);
while(<F>) {
    chomp;
    my ($label, $value) = split(/\s+/);
    my ($clabel, $rlabel) = split(/:/, $label);

	$label_cutoff{$clabel}{$rlabel} = $value;
}
close(F);

my %bid_pairs = ();
my %sr_link_scores = ();
my %sr_link_counts = ();
my %sr_stat = ();

if ($type =~ /SR/) {
    read_linkscores("$raca_dir/scf_block_linkscores.SR.bid.txt", \%sr_link_scores, \%sr_link_counts, \%sr_stat, \%bid_pairs);
}

my %linkr_link_scores = ();
my %linkr_link_counts = ();
my %linkr_stat = ();

if ($type =~ /LnR/) {
    read_linkscores("$raca_dir/scf_block_linkscores.10xR.bid.txt", \%linkr_link_scores, \%linkr_link_counts, \%linkr_stat, \%bid_pairs);
}

my %lr_link_scores = ();
my %lr_link_counts = ();
my %lr_stat = ();
if ($type =~ /LR/) {
    read_linkscores("$raca_dir/scf_block_linkscores.LR.bid.txt", \%lr_link_scores, \%lr_link_counts, \%lr_stat, \%bid_pairs);
}

if ($type =~ /SR/) { open(OSR, ">$sr_linkoutbid_f"); }
if ($type =~ /LnR/) { open(OLnR, ">$linkr_linkoutbid_f"); }
if ($type =~ /LR/) { open(OLR, ">$lr_linkoutbid_f"); }

my %used_pairs = ();
foreach my $pair (sort keys %bid_pairs) {
    if (defined($used_pairs{$pair})) { next; }
    $used_pairs{$pair} = 1;

	my ($bid1, $bid2) = split(/:/, $pair);
    my $rkey = -1*$bid2 . ":" . -1*$bid1;	
    $used_pairs{$rkey} = 1;
    
    my $clabel = "CN";

    my $sr_count = $sr_link_counts{$pair};
    my $sr_score = $sr_link_scores{$pair};
    my $sr_label = get_rlabel($sr_count, \%sr_stat);
    
    my $linkr_count = $linkr_link_counts{$pair};
    my $linkr_score = $linkr_link_scores{$pair};
    my $linkr_label = get_rlabel($linkr_count, \%linkr_stat);

    my $lr_count = $lr_link_counts{$pair};
    my $lr_score = $lr_link_scores{$pair};
    my $lr_label = get_rlabel($lr_count, \%lr_stat);

    if ($type eq "SR") {
        if ($label_cutoff{$clabel}{$sr_label} eq "yes") {
            if (defined($sr_score)) {
                print OSR "$bid1\t$bid2\t$sr_score\t$sr_label\n"; 
            }
        } 
	} elsif ($type eq "LnR") {
        if ($label_cutoff{$clabel}{$linkr_label} eq "yes") {
            if (defined($linkr_score)) {
                print OLnR "$bid1\t$bid2\t$linkr_score\t$linkr_label\n"; 
            }
        } 
    } elsif ($type eq "LR") {
        if ($label_cutoff{$clabel}{$lr_label} eq "yes") {
            if (defined($lr_score)) {
                print OLR "$bid1\t$bid2\t$lr_score\t$lr_label\n"; 
            }
        } 
    } elsif ($type eq "SR:LR") {
        if ($label_cutoff{$clabel}{$sr_label} eq "yes" || $label_cutoff{$clabel}{$lr_label} eq "yes") {
            if (defined($sr_score)) {
                print OSR "$bid1\t$bid2\t$sr_score\t$sr_label\n"; 
            }
            
            if (defined($lr_score)) {
                print OLR "$bid1\t$bid2\t$lr_score\t$lr_label\n"; 
            }
        }
    } elsif ($type eq "SR:LnR") {
        if ($label_cutoff{$clabel}{$sr_label} eq "yes" || $label_cutoff{$clabel}{$linkr_label} eq "yes") {
            if (defined($sr_score)) {
                print OSR "$bid1\t$bid2\t$sr_score\t$sr_label\n"; 
            }
            
            if (defined($linkr_score)) {
                print OLnR "$bid1\t$bid2\t$linkr_score\t$linkr_label\n"; 
            }
        }
    } elsif ($type eq "LnR:LR") {
        if ($label_cutoff{$clabel}{$linkr_label} eq "yes" || $label_cutoff{$clabel}{$lr_label} eq "yes") {
            if (defined($linkr_score)) {
                print OLnR "$bid1\t$bid2\t$linkr_score\t$linkr_label\n"; 
            }
            
            if (defined($lr_score)) {
                print OLR "$bid1\t$bid2\t$lr_score\t$lr_label\n"; 
            }
        }
	}
}
close(F);
if ($type =~ /SR/) { close(OSR); }
if ($type =~ /LnR/) { close(OLnR); }
if ($type =~ /LR/) { close(OLR); }

#################################################################
sub get_clabel {
    my $i_cnt = shift;
    my $o_cnt = shift;
    
    my $c_label = "";
    if ($i_cnt == 0) {
        if ($o_cnt >= 1) { 
            $c_label = "CW"; 
        } else {
            $c_label = "CN";
        }
    } elsif ($i_cnt == 1) {
        if ($o_cnt == 0) { 
            $c_label = "CW";  
        } else {
            $c_label = "CM";
        }
    } elsif ($i_cnt == 2) {
        if ($o_cnt == 0) {
            $c_label = "CM";
        } else {
            $c_label = "CS"; 
        }
    } else {
        $c_label = "CS";
    }

    return $c_label;
}

sub get_rlabel {
    my $value = shift;
    my $rhs_stat = shift;

    my $rlabel = "RN";
    if (!defined($value) || $value == 0) { return $rlabel; }
    
    my $ml = $$rhs_stat{"#raw_2.5:"};
    my $sl = $$rhs_stat{"raw_16:"};
    my $sr = $$rhs_stat{"raw_84:"};
    my $mr = $$rhs_stat{"raw_97.5:"};

    if ($value > $sl && $value < $sr) { 
        $rlabel = "RS";
    } elsif ($value > $ml && $value < $mr) { 
        $rlabel = "RM"; 
    } else { 
        $rlabel = "RW"; 
    }

    return $rlabel;
}

sub read_linkscores {
    my $f = shift;
    my $rhs_link_scores = shift;
    my $rhs_link_counts = shift;
    my $rhs_stat = shift;
    my $rhs_bid_pairs = shift;

    my $out = `head -1 $f`;
    chomp($out);
    my @ar = split(/\s+/, $out);
    for (my $i = 0; $i < $#ar; $i+=2) {
        my $key = $ar[$i];
        my $val = $ar[$i+1];
        $$rhs_stat{$key} = $val;
    }
    
    open(F, $f);
    while(<F>) {
	    chomp;
	    if ($_ =~ /^#/ || $_ eq "") { next; }
	    my ($bid1, $bid2, $score, $count) = split(/\s+/);

	    my ($rbid1, $rbid2) = (-1*$bid1, -1*$bid2);
	    $$rhs_link_scores{"$bid1:$bid2"} = $score;
	    $$rhs_link_scores{"$rbid2:$rbid1"} = $score;
	    
        $$rhs_link_counts{"$bid1:$bid2"} = $count;
	    $$rhs_link_counts{"$rbid2:$rbid1"} = $count;

        $$rhs_bid_pairs{"$bid1:$bid2"} = 1;
        $$rhs_bid_pairs{"$rbid2:$rbid1"} = 1;
    }
    close(F);
}
