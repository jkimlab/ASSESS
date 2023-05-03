#!/usr/bin/perl

use strict;
use warnings;

my $type = shift;
my $tspc = shift;
my $sr_inter_f = shift;
my $sr_intra_f = shift;
my $linkr_inter_f = shift;
my $linkr_intra_f = shift;
my $lr_inter_f = shift;
my $lr_intra_f = shift;
my $raca_dir = shift;

my %cons_scores = ();
open(F, "$raca_dir/SFs/adjacencies.prob");
while(<F>) {
	chomp;
	if ($_ =~ /^#/ || $_ eq "") { next; }
	my ($bid1, $bid2, $score) = split(/\s+/);
	my ($rbid1, $rbid2) = (-1*$bid1, -1*$bid2);
	$cons_scores{"$bid1:$bid2"} = $score;
	$cons_scores{"$rbid2:$rbid1"} = $score;
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

my $aout = `head -1 $raca_dir/SFs/adjacencies.prob`;
my $numsfs = substr($aout, 1);

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

my %sr_stat = ();
my %sr_link_scores = ();
my %sr_inter_reads = ();
my %sr_intra_reads = ();
if ($type =~ /SR/) {
    read_link_data("$raca_dir/SFs/block_linkscores.SR.bid.txt", $sr_inter_f, $sr_intra_f, \%sr_stat, \%sr_link_scores, 
                    \%sr_inter_reads, \%sr_intra_reads);
}

my %linkr_stat = ();
my %linkr_link_scores = ();
my %linkr_inter_reads = ();
my %linkr_intra_reads = ();
if ($type =~ /LnR/) {
    read_link_data("$raca_dir/SFs/block_linkscores.10xR.bid.txt", $linkr_inter_f, $linkr_intra_f, \%linkr_stat, \%linkr_link_scores, 
                    \%linkr_inter_reads, \%linkr_intra_reads);
}

my %lr_stat = ();
my %lr_link_scores = ();
my %lr_inter_reads = ();
my %lr_intra_reads = ();
if ($type =~ /LR/) {
    read_link_data("$raca_dir/SFs/block_linkscores.LR.bid.txt", $lr_inter_f, $lr_intra_f, \%lr_stat, \%lr_link_scores, 
                    \%lr_inter_reads, \%lr_intra_reads);
}

my @plines = ();
for (my $bid1 = -1 * $numsfs; $bid1 <= $numsfs; $bid1++) {
    if ($bid1 == 0) { next; }

    for (my $bid2 = -1 * $numsfs; $bid2 <= $numsfs; $bid2++) {
        if ($bid2 == 0) { next; }
        if (abs($bid1) == abs($bid2)) { next; }

        if ($type eq "CG") {
            if (!defined($cons_scores{"$bid1:$bid2"})) { next; }
        } elsif ($type eq "CG:SR" || $type eq "SR") {
            if (!defined($cons_scores{"$bid1:$bid2"}) && !defined($sr_link_scores{"$bid1:$bid2"})) { next; }
        } elsif ($type eq "CG:LnR" || $type eq "LnR") {
            if (!defined($cons_scores{"$bid1:$bid2"}) && !defined($linkr_link_scores{"$bid1:$bid2"})) { next; }
        } elsif ($type eq "CG:LR" || $type eq "LR") {
            if (!defined($cons_scores{"$bid1:$bid2"}) && !defined($lr_link_scores{"$bid1:$bid2"})) { next; }
        } elsif ($type eq "CG:SR:LR" || $type eq "SR:LR") {
            if (!defined($cons_scores{"$bid1:$bid2"}) && !defined($sr_link_scores{"$bid1:$bid2"}) &&
                    !defined($lr_link_scores{"$bid1:$bid2"})) { next; }
        } elsif ($type eq "CG:SR:LnR" || $type eq "SR:LnR") {
            if (!defined($cons_scores{"$bid1:$bid2"}) && !defined($sr_link_scores{"$bid1:$bid2"}) &&
                    !defined($linkr_link_scores{"$bid1:$bid2"})) { next; }
        } elsif ($type eq "CG:LnR:LR" || $type eq "LnR:LR") {
            if (!defined($cons_scores{"$bid1:$bid2"}) && !defined($linkr_link_scores{"$bid1:$bid2"}) &&
                    !defined($lr_link_scores{"$bid1:$bid2"})) { next; }
        } elsif ($type eq "CG:SR:LnR:LR" || $type eq "SR:LnR:LR") {
            if (!defined($cons_scores{"$bid1:$bid2"}) && !defined($sr_link_scores{"$bid1:$bid2"}) &&
                    !defined($linkr_link_scores{"$bid1:$bid2"}) && !defined($lr_link_scores{"$bid1:$bid2"})) { next; }
		}

		my $line = "";
		my ($ns1, $ns2) = ("S", "S");
		if (defined($ns_blocks{abs($bid1)})) { $ns1 = "N"; }
		if (defined($ns_blocks{abs($bid2)})) { $ns2 = "N"; }
		$line .= "$bid1($ns1):$bid2($ns2)\t";

		my $cscore = $cons_scores{"$bid1:$bid2"};
		if (defined($cscore)) {
			$line .= sprintf "%.20f\t", $cscore;
		} else {
			$line .= "N/A\t";
			$cscore = 0;
		}

		my $icnt = 0;
        my $ocnt = 0;
		foreach my $spc (@spcs) {
			my $join = $joins{$spc}{"$bid1:$bid2"};	
			if (defined($join)) {
				$line .= "O\t";

                if ($spc eq $tspc) {  }
                elsif ($hs_spcs{$spc} == 2) { $ocnt++; }
                else { $icnt++; }
			} else {
				$line .= "X\t";
			}
		}
		$line .= "$icnt\t$ocnt\t";

		my ($rchr1, $rchr2) = ($block_rchr{abs($bid1)}, $block_rchr{abs($bid2)});
		if ($rchr1 eq $rchr2) {	
			$line .= "INTRA:$rchr1:$rchr2\t";
		} else {
			$line .= "INTER:$rchr1:$rchr2\t";
		}

        my $sr_label = "RN";
        my $sr_score = 0;
        if ($type =~ /SR/) {
		    $sr_score = $sr_link_scores{"$bid1:$bid2"};
		    if (defined($sr_score)) {
			    $line .= sprintf "%.20f\t", $sr_score;
		    } else {
			    $line .= "N/A\t";
		    }
        
		    my $count = $sr_inter_reads{"$bid1:$bid2"};
		    if (defined($count)) {
			    $line .= "inter:$count\t";
		    } else {
			    $count = $sr_intra_reads{"$bid1:$bid2"};
			    if (defined($count)) {
                    $line .= "intra:$count\t";
			    } else {
				    $line .= "N/A\t";
			    }
		    }

            $sr_label = get_rlabel($count, \%sr_stat);
        }
       
        my $linkr_label = "RN";
        my $linkr_score = 0;
        if ($type =~ /LnR/) {
		    $linkr_score = $linkr_link_scores{"$bid1:$bid2"};
		    if (defined($linkr_score)) {
			    $line .= sprintf "%.20f\t", $linkr_score;
		    } else {
			    $line .= "N/A\t";
		    }
        
		    my $count = $linkr_inter_reads{"$bid1:$bid2"};
		    if (defined($count)) {
			    $line .= "inter:$count\t";
		    } else {
			    $count = $linkr_intra_reads{"$bid1:$bid2"};
			    if (defined($count)) {
                    $line .= "intra:$count\t";
			    } else {
				    $line .= "N/A\t";
			    }
		    }

            $linkr_label = get_rlabel($count, \%linkr_stat);
        }

        my $lr_label = "RN";
        my $lr_score = 0;
        if ($type =~ /LR/) {
		    $lr_score = $lr_link_scores{"$bid1:$bid2"};
		    if (defined($lr_score)) {
			    $line .= sprintf "%.20f\t", $lr_score;
		    } else {
			    $line .= "N/A\t";
		    }
		    
            my $count = $lr_inter_reads{"$bid1:$bid2"};
		    if (defined($count)) {
			    $line .= "inter:$count\t";
		    } else {
			    $count = $lr_intra_reads{"$bid1:$bid2"};
			    if (defined($count)) {
                    $line .= "intra:$count\t";
			    } else {
				    $line .= "N/A\t";
			    }
		    }

            $lr_label = get_rlabel($count, \%lr_stat);
        }

        my $clabel = get_clabel($icnt, $ocnt);
        $line .= "$clabel\t";

        if ($type =~ /SR/) { $line .= "$sr_label\t"; }
        if ($type =~ /LnR/) { $line .= "$linkr_label\t"; }
        if ($type =~ /LR/) { $line .= "$lr_label\t"; }
		
        my $fscore = -1;
        if ($type eq "CG") {
            if (defined($cscore)) { $fscore = $cscore; }
        } elsif ($type eq "CG:SR" || $type eq "SR") {
            if (defined($cscore)) { $fscore = $cscore; }
            if (defined($sr_score)) { $fscore += $sr_score; }
            if ($fscore >= 0) { $fscore = $fscore/2; }
        } elsif ($type eq "CG:LnR" || $type eq "LnR") {
            if (defined($cscore)) { $fscore = $cscore; }
            if (defined($linkr_score)) { $fscore += $linkr_score; }
            if ($fscore >= 0) { $fscore = $fscore/2; }
        } elsif ($type eq "CG:LR" || $type eq "LR" ) {
            if (defined($cscore)) { $fscore = $cscore; }
            if (defined($lr_score)) { $fscore += $lr_score; }
            if ($fscore >= 0) { $fscore = $fscore/2; }
        } elsif ($type eq "CG:SR:LR" || $type eq "SR:LR") {
            if (defined($cscore)) { $fscore = $cscore; }
            if (defined($sr_score)) { $fscore += $sr_score; }
            if (defined($lr_score)) { $fscore += $lr_score; }
            if ($fscore >= 0) { $fscore = $fscore/3; }
        } elsif ($type eq "CG:SR:LnR" || $type eq "SR:LnR") {
            if (defined($cscore)) { $fscore = $cscore; }
            if (defined($sr_score)) { $fscore += $sr_score; }
            if (defined($linkr_score)) { $fscore += $linkr_score; }
            if ($fscore >= 0) { $fscore = $fscore/3; }
        } elsif ($type eq "CG:LnR:LR" || $type eq "LnR:LR") {
            if (defined($cscore)) { $fscore = $cscore; }
            if (defined($linkr_score)) { $fscore += $linkr_score; }
            if (defined($lr_score)) { $fscore += $lr_score; }
            if ($fscore >= 0) { $fscore = $fscore/3; }
        } elsif ($type eq "CG:SR:LnR:LR" || $type eq "SR:LnR:LR") {
            if (defined($cscore)) { $fscore = $cscore; }
            if (defined($sr_score)) { $fscore += $sr_score; }
            if (defined($linkr_score)) { $fscore += $linkr_score; }
            if (defined($lr_score)) { $fscore += $lr_score; }
            if ($fscore >= 0) { $fscore = $fscore/4; }
        }

		$line .= sprintf "%.20f\t", $fscore;

		push(@plines, $line);
	}
}
close(F);

open(O, ">$raca_dir/adjinfo.all.txt");
print O "#bid1:bid2\tcons_score\t";
foreach my $spc (@spcs) {
	print O "$spc\t";
}

print O "num_in_agr\tnum_out_agr\tref_chrs\t";

if ($type =~ /SR/) { print O "SR_link_score\tSR_num_reads\t"; }
if ($type =~ /LnR/) { print O "10xR_link_score\t10xR_num_reads\t"; }
if ($type =~ /LR/) { print O "LR_link_score\tLR_num_reads\t"; }

print O "cg_label\t";

if ($type =~ /SR/) { print O "SR_label\t"; }
if ($type =~ /SR/) { print O "10xR_label\t"; }
if ($type =~ /LR/) { print O "LR_label\t"; }

print O "final_score\n";

foreach my $line (@plines) {
	print O "$line\n";
}
close(O);

open(O, ">$raca_dir/adjinfo_sorted.all.txt");
print O "#bid1:bid2\tcons_score\t";
foreach my $spc (@spcs) {
	print O "$spc\t";
}

print O "num_in_agr\tnum_out_agr\tref_chrs\t";

if ($type =~ /SR/) { print O "SR_link_score\tSR_num_reads\t"; }
if ($type =~ /LnR/) { print O "10xR_link_score\t10xR_num_reads\t"; }
if ($type =~ /LR/) { print O "LR_link_score\tLR_num_reads\t"; }

print O "cg_label\t";

if ($type =~ /SR/) { print O "SR_label\t"; }
if ($type =~ /SR/) { print O "10xR_label\t"; }
if ($type =~ /LR/) { print O "LR_label\t"; }

print O "final_score\n";

my @exlines = ();
foreach my $line (sort by_adjscore @plines) {
	print O "$line\n";
	my @ar = split(/\s+/, $line);
	
	if ($ar[3] eq "X") { 
		push(@exlines, $line);  
	}
}
close(O);

open(O, ">$raca_dir/adjinfo_nonref.all.txt");
print O "#bid1:bid2\tcons_score\t";
foreach my $spc (@spcs) {
	print O "$spc\t";
}

print O "num_in_agr\tnum_out_agr\tref_chrs\t";

if ($type =~ /SR/) { print O "SR_link_score\tSR_num_reads\t"; }
if ($type =~ /LnR/) { print O "10xR_link_score\t10xR_num_reads\t"; }
if ($type =~ /LR/) { print O "LR_link_score\tLR_num_reads\t"; }

print O "cg_label\t";

if ($type =~ /SR/) { print O "SR_label\t"; }
if ($type =~ /SR/) { print O "10xR_label\t"; }
if ($type =~ /LR/) { print O "LR_label\t"; }

print O "final_score\n";

foreach my $line (@exlines) {
	print O "$line\n";
}

close(O);

sub by_adjscore {
	my @ar1 = split(/\s+/, $a); 
	my @ar2 = split(/\s+/, $b); 

	my ($score1, $score2) = ($ar1[-1], $ar2[-1]);

	return ($score1 <=> $score2); 
}

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

sub read_link_data {
    my $f = shift;
    my $inter_f = shift;
    my $intra_f = shift;
    my $rhs_stat = shift;
    my $rhs_scores = shift;
    my $rhs_inter = shift;
    my $rhs_intra = shift;

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
        $$rhs_scores{"$bid1:$bid2"} = $score;
        $$rhs_scores{"$rbid2:$rbid1"} = $score;
    }
    close(F);

    open(F, "$inter_f");
    while(<F>) {
        chomp;
        my ($scf1, $dir1, $scf2, $dir2, $score, $bid1, $bid2) = split(/\s+/);

        my ($rbid1, $rbid2) = (-1*$bid1, -1*$bid2);
        $$rhs_inter{"$bid1:$bid2"} = $score;
        $$rhs_inter{"$rbid2:$rbid1"} = $score;
    }
    close(F);

    open(F, "$intra_f");
    while(<F>) {
        chomp;
        my @ar = split(/\s+/);
        my ($scf1, $dir1, $scf2, $dir2) = ($ar[0], $ar[1], $ar[2], $ar[3]);
        $ar[7] =~ /MIN=(\S+)/;
        my $score = $1;

        my ($bid1, $bid2) = ($ar[$#ar-1], $ar[$#ar]);
        my ($rbid1, $rbid2) = (-1*$bid1, -1*$bid2);
        $$rhs_intra{"$bid1:$bid2"} = $score;
        $$rhs_intra{"$rbid2:$rbid1"} = $score;
    }
    close(F);
}
