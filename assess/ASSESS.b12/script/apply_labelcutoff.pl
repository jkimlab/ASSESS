#!/usr/bin/perl

use strict;
use warnings;

my $type = shift;
my $tspc = shift;
my $labelcutoff_f = shift;
my $raca_dir = shift;

my $consout_f = "$raca_dir/SFs/block_consscores.flt.bid.txt";
my $sr_linkoutbid_f = "$raca_dir/SFs/block_linkscores.SR.flt.bid.txt";
my $linkr_linkoutbid_f = "$raca_dir/SFs/block_linkscores.10xR.flt.bid.txt";
my $lr_linkoutbid_f = "$raca_dir/SFs/block_linkscores.LR.flt.bid.txt";

my %label_cutoff = ();
open(F, $labelcutoff_f);
while(<F>) {
    chomp;
    my ($label, $value) = split(/\s+/);
    my ($clabel, $rlabel) = split(/:/, $label);

    if ($type ne "CG") {
        $label_cutoff{$clabel}{$rlabel} = $value;
    } else {
        $label_cutoff{$clabel} = $value;
    }
}
close(F);

my %bid_pairs = ();

my %cons_scores = ();
open(F, "$raca_dir/SFs/block_consscores.bid.txt");
while(<F>) {
	chomp;
	if ($_ =~ /^#/ || $_ eq "") { next; }
	my ($bid1, $bid2, $score) = split(/\s+/);
	my ($rbid1, $rbid2) = (-1*$bid1, -1*$bid2);
	$cons_scores{"$bid1:$bid2"} = $score;
	$cons_scores{"$rbid2:$rbid1"} = $score;
	
    $bid_pairs{"$bid1:$bid2"} = 1;
	$bid_pairs{"$rbid2:$rbid1"} = 1;
}
close(F);

my %ns_blocks = ();
my $tflag = 0;
open(F, "$raca_dir/SFs/Genomes.Order");
while(<F>) {
	chomp;
	if ($_ =~ /^>(\S+)/) {
		if ($1 eq $tspc) { $tflag = 1; }	
		else { $tflag = 0; }
	} elsif (length($_) == 0) {
		if ($tflag == 1) { $tflag = 0; }
	} elsif ($_ !~ /^#/ && $tflag == 1) {
		my @bids = split(/\s+/);
		pop(@bids);
		if (scalar(@bids) > 1) {
			foreach my $bid (@bids) {
				$ns_blocks{abs($bid)} = 1;		
			}
		}
	}
}
close(F);

my %hs_spcs = ();
my @spcs = ();
open(F, "$raca_dir/SFs/config.file");
while(<F>) {
	chomp;
	if ($_ =~ /^(\S+)\s+(\d)$/) {
		my ($spc, $ord) = ($1, $2);	
		$hs_spcs{$spc} = $ord;
		push(@spcs, $spc);
	}
}
close(F);

my %joins = ();
foreach my $spc (@spcs) {
	open(F, "$raca_dir/SFs/$spc.joins");
	while(<F>) {
		chomp;
		if ($_ =~ /^#/) { next; }
		my $line = $_;
		$line =~ s/^\s+//;
		$line =~ s/\s+$//;
		my ($bid1, $bid2) = split(/\s+/, $line);
		my ($rbid1, $rbid2) = (-1*$bid1, -1*$bid2);
		$joins{$spc}{"$bid1:$bid2"} = 1;
		$joins{$spc}{"$rbid2:$rbid1"} = 1;
	}
	close(F);
}

my %block_rchr = ();
open(F, "$raca_dir/SFs/block_list.txt");
while(<F>) {
	chomp;
	my ($scf, $start, $end, $dir, $bid, $rchr) = split(/\s+/);
	$block_rchr{$bid} = $rchr;
}
close(F);

my %sr_link_scores = ();
my %sr_link_counts = ();
my %sr_stat = ();
if ($type =~ /SR/) {
    read_linkscores("$raca_dir/SFs/block_linkscores.SR.bid.txt", \%sr_link_scores, \%sr_link_counts, \%sr_stat, \%bid_pairs);
}

my %linkr_link_scores = ();
my %linkr_link_counts = ();
my %linkr_stat = ();
if ($type =~ /LnR/) {
    read_linkscores("$raca_dir/SFs/block_linkscores.10xR.bid.txt", \%linkr_link_scores, \%linkr_link_counts, \%linkr_stat, \%bid_pairs);
}

my %lr_link_scores = ();
my %lr_link_counts = ();
my %lr_stat = ();
if ($type =~ /LR/) {
    read_linkscores("$raca_dir/SFs/block_linkscores.LR.bid.txt", \%lr_link_scores, \%lr_link_counts, \%lr_stat, \%bid_pairs);
}

open(OC, ">$consout_f");
if ($type =~ /SR/ ) { open(OSR, ">$sr_linkoutbid_f"); }
if ($type =~ /LnR/ ) { open(OLnR, ">$linkr_linkoutbid_f"); }
if ($type =~ /LR/ ) { open(OLR, ">$lr_linkoutbid_f"); }

my %used_pairs = ();
foreach my $pair (sort keys %bid_pairs) {
    if (defined($used_pairs{$pair})) { next; }
    $used_pairs{$pair} = 1;

	my ($bid1, $bid2) = split(/:/, $pair);
    my $rkey = -1*$bid2 . ":" . -1*$bid1;	
    $used_pairs{$rkey} = 1;
    
    my $cscore = $cons_scores{$pair};
    my $clabel = "CN";
    if (defined($cscore)) {
        my $icnt = 0;
        my $ocnt = 0;
	    foreach my $spc (@spcs) {
		    my $join = $joins{$spc}{$pair};	
		    if (defined($join)) {
                if ($spc eq $tspc) {  }
                elsif ($hs_spcs{$spc} == 2) { $ocnt++; }
                else { $icnt++; }
    	    } 
	    }

        $clabel = get_clabel($icnt, $ocnt);
    }

    my $sr_count = $sr_link_counts{$pair};
    my $sr_score = $sr_link_scores{$pair};
    my $sr_label = get_rlabel($sr_count, \%sr_stat);

    my $linkr_count = $linkr_link_counts{$pair};
    my $linkr_score = $linkr_link_scores{$pair};
    my $linkr_label = get_rlabel($linkr_count, \%linkr_stat);
    
    my $lr_count = $lr_link_counts{$pair};
    my $lr_score = $lr_link_scores{$pair};
    my $lr_label = get_rlabel($lr_count, \%lr_stat);

    if ($type eq "CG") {
        if ($label_cutoff{$clabel} eq "yes") {
            print OC "$bid1\t$bid2\t$cscore\t$clabel\n";
        }
    } elsif ($type eq "CG:SR") {
        if ($label_cutoff{$clabel}{$sr_label} eq "yes") {
            if (defined($cscore)) {
                print OC "$bid1\t$bid2\t$cscore\t$clabel\n";
            }

            if (defined($sr_score)) {
                print OSR "$bid1\t$bid2\t$sr_score\t$sr_label\n"; 
            }
        } 
    } elsif ($type eq "CG:LnR") {
        if ($label_cutoff{$clabel}{$linkr_label} eq "yes") {
            if (defined($cscore)) {
                print OC "$bid1\t$bid2\t$cscore\t$clabel\n";
            }

            if (defined($linkr_score)) {
                print OLnR "$bid1\t$bid2\t$linkr_score\t$linkr_label\n"; 
            }
        } 
    } elsif ($type eq "CG:LR") {
        if ($label_cutoff{$clabel}{$lr_label} eq "yes") {
            if (defined($cscore)) {
                print OC "$bid1\t$bid2\t$cscore\t$clabel\n";
            }

            if (defined($lr_score)) {
                print OLR "$bid1\t$bid2\t$lr_score\t$lr_label\n"; 
            }
        } 
	} elsif ($type eq "CG:SR:LR") {
        if ($label_cutoff{$clabel}{$sr_label} eq "yes" || $label_cutoff{$clabel}{$lr_label} eq "yes") {
            if (defined($cscore)) {
                print OC "$bid1\t$bid2\t$cscore\t$clabel\n";
            }

            if (defined($sr_score)) {
                print OSR "$bid1\t$bid2\t$sr_score\t$sr_label\n"; 
            }
            
            if (defined($lr_score)) {
                print OLR "$bid1\t$bid2\t$lr_score\t$lr_label\n"; 
            }
        }
	} elsif ($type eq "CG:SR:LnR") {
        if ($label_cutoff{$clabel}{$sr_label} eq "yes" || $label_cutoff{$clabel}{$linkr_label} eq "yes") {
            if (defined($cscore)) {
                print OC "$bid1\t$bid2\t$cscore\t$clabel\n";
            }

            if (defined($sr_score)) {
                print OSR "$bid1\t$bid2\t$sr_score\t$sr_label\n"; 
            }
            
            if (defined($linkr_score)) {
                print OLnR "$bid1\t$bid2\t$linkr_score\t$linkr_label\n"; 
            }
        }
	} elsif ($type eq "CG:LnR:LR") {
        if ($label_cutoff{$clabel}{$linkr_label} eq "yes" || $label_cutoff{$clabel}{$lr_label} eq "yes") {
            if (defined($cscore)) {
                print OC "$bid1\t$bid2\t$cscore\t$clabel\n";
            }

            if (defined($linkr_score)) {
                print OLnR "$bid1\t$bid2\t$linkr_score\t$linkr_label\n"; 
            }
            
            if (defined($lr_score)) {
                print OLR "$bid1\t$bid2\t$lr_score\t$lr_label\n"; 
            }
        }
    } else { 	# CG:SR:LnR:LR type
        if ($label_cutoff{$clabel}{$sr_label} eq "yes" || $label_cutoff{$clabel}{$linkr_label} eq "yes" || $label_cutoff{$clabel}{$lr_label} eq "yes") {
            if (defined($cscore)) {
                print OC "$bid1\t$bid2\t$cscore\t$clabel\n";
            }

            if (defined($sr_score)) {
                print OSR "$bid1\t$bid2\t$sr_score\t$sr_label\n"; 
            }

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
close(OC);
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
