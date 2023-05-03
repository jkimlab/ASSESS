#!/usr/bin/perl

use strict;
use warnings;
use FindBin qw($Bin);

my $samtools = "$Bin/../util/samtools-1.5/samtools";

my $map_f = shift;		# sam or bam file
my $size_f = shift;
my $new_f = shift;
my $min_mapq = shift;

my ($fr, $rf, $ee) = (0, 0 ,0);
my %fr_tlens = ();
my %rf_tlens = ();
my ($fr_sum, $fr_2_sum) = (0, 0);
my ($rf_sum, $rf_2_sum) = (0, 0);
my $no_skipped_rname = 0;
my $no_skipped_crd_samestart = 0;
my $no_skipped_nonscf = 0;
my $total = 0;

my %scfnames = ();
open(F, $size_f);
while(<F>) {
    chomp;
    my ($name, $size) = split(/\s+/);
    $scfnames{$name} = 1;
}
close(F);

my %unpaired = ();
open(F, "$samtools view -F 0xF0C -q $min_mapq $map_f |");
while(<F>) {
	chomp;
	if ($_ =~ /^@/) { next; }
	my ($qname, $flag) = split(/\s+/);
	if ($flag & 0x4 || $flag & 0x8) { next; }

    if (defined($unpaired{$qname})) {
        delete $unpaired{$qname};
    } else {
        $unpaired{$qname} = 1;
    }
}
close(F);

my $upcnt = 0;
my $pcnt = 0;
open(F, "$samtools view -F 0xF0C -q $min_mapq $map_f |");
while(<F>) {
	chomp;
	if ($_ =~ /^@/) { next; }
	my ($qname, $flag, $rname, $pos, $mapq, $cigar, $rnext, $pnext, $tlen) = split(/\s+/);
	if ($flag & 0x80 || $flag & 0x4 || $flag & 0x8) { next; }
    $pcnt++; 
    if (defined($unpaired{$qname})) { $upcnt++; next; }

	my ($dir1, $dir2) = ("+", "+");
	if ($flag & 0x10) { $dir1 = "-"; }
	if ($flag & 0x20) { $dir2 = "-"; }
	
	my ($rname1, $rname2) = ($rname, $rnext);
	if ($rname2 eq "=") { $rname2 = $rname1; }

	$total++;

	if ($rname1 ne $rname2) {
		$no_skipped_rname++;
		next;	
	}

	if (!defined($scfnames{$rname1})) {
		$no_skipped_nonscf++;
		next;
	}

    if ($pos == $pnext) {
		$no_skipped_crd_samestart++;
		next;	
	}

	if ($dir1 eq "+" && $dir2 eq "-") {
		if ($pos < $pnext) {
			$fr++;
            my $abs_tlen = abs($tlen);
			$fr_tlens{$abs_tlen} = (exists $fr_tlens{$abs_tlen}) ? $fr_tlens{$abs_tlen} + 1 : 1;
            $fr_sum += $abs_tlen;
            $fr_2_sum += ($abs_tlen * $abs_tlen); 
		 } else {
		    $rf++;
            my $abs_tlen = abs($tlen);
			$rf_tlens{$abs_tlen} = (exists $rf_tlens{$abs_tlen}) ? $rf_tlens{$abs_tlen} + 1 : 1;
            $rf_sum += $abs_tlen;
            $rf_2_sum += ($abs_tlen * $abs_tlen);
	     }
	} elsif ($dir1 eq "-" && $dir2 eq "+") {
		if ($pos > $pnext) {
			$fr++;
            my $abs_tlen = abs($tlen);
			$fr_tlens{$abs_tlen} = (exists $fr_tlens{$abs_tlen}) ? $fr_tlens{$abs_tlen} + 1 : 1;
            $fr_sum += $abs_tlen;
            $fr_2_sum += ($abs_tlen * $abs_tlen); 
        } else {
            $rf++;
            my $abs_tlen = abs($tlen);
			$rf_tlens{$abs_tlen} = (exists $rf_tlens{$abs_tlen}) ? $rf_tlens{$abs_tlen} + 1 : 1;
            $rf_sum += $abs_tlen;
            $rf_2_sum += ($abs_tlen * $abs_tlen);
        }
	} else {
		$ee++;
	}
}
close(F);

print STDERR "No. total: $pcnt\n";
print STDERR "No. unpaired: $upcnt\n";

my ($fr_mean, $fr_stdev) = get_stats($fr, $fr_sum, $fr_2_sum);
my ($rf_mean, $rf_stdev) = get_stats($rf, $rf_sum, $rf_2_sum);
my $fr_median = cal_median(\%fr_tlens, $fr);
my $rf_median = cal_median(\%rf_tlens, $rf);

my $f_rname = $no_skipped_rname / $total;
my $f_crds = $no_skipped_crd_samestart / $total;
my $f_nonscf = $no_skipped_nonscf / $total;

print STDERR "[$map_f]: Skipped reads out of total ($total)\n";
print STDERR "\tType\tCount (%)\n";
printf STDERR "\tRNAME\t$no_skipped_rname (%.2f %%)\n", $f_rname*100;
printf STDERR "\tCRD_SAMESTART\t$no_skipped_crd_samestart (%.2f %%)\n", $f_crds*100;
printf STDERR "\tNOTUSED_SCF\t$no_skipped_nonscf (%.2f %%)\n", $f_nonscf*100;

my $r_total = $fr + $rf + $ee;
my ($f_fr, $f_rf, $f_ee) = (0, 0, 0);
if ($r_total > 0) {
	$f_fr = $fr/$r_total;
	$f_rf = $rf/$r_total;
	$f_ee = $ee/$r_total;
}
print STDERR "\n[$map_f]: Read orientation stats.\n";
print STDERR "\tOT\tCount (%)\tMean_length (Stdev)\tMedian_length\n";
if ($fr > 0) { printf STDERR "\tFR\t%d (%.2f %%)\t%.2f (%.2f)\t%.2f\n", $fr, $f_fr*100, $fr_mean, $fr_stdev, $fr_median; }
else { printf STDERR "\tFR\t%d (%.2f %%)\t%.2f (%.2f)\t%.2f\n", 0, 0, 0, 0, 0; }
if ($rf > 0) { printf STDERR "\tRF\t%d (%.2f %%)\t%.2f (%.2f)\t%.2f\n", $rf, $f_rf*100, $rf_mean, $rf_stdev, $rf_median; }
else { printf STDERR "\tRF\t%d (%.2f %%)\t%.2f (%.2f)\t%.2f\n", 0, 0, 0, 0, 0; }
if ($ee > 0) { printf STDERR "\tOther\t%d (%.2f %%)\n", $ee, $f_ee*100; }

my $final_ot = ($fr > $rf) ? "fr" : "rf";
my $final_ot_reads = ($final_ot eq "fr") ? $fr : $rf;

my ($max_len, $min_len) = (0, 0);
my $pass_cnt = 0;
my $pass_ratio = 0;

for (my $len_ext = 0.25; $len_ext <= 1; $len_ext += 0.05) {
	my $cnt = 0;
	if ($final_ot eq "fr") {
		$max_len = $fr_median * (1 + $len_ext);
		$min_len = $fr_median * (1 - $len_ext);
		foreach my $tlen (sort {$a <=> $b} keys %fr_tlens) {
			$cnt += ($tlen >= $min_len && $tlen <= $max_len) ? $fr_tlens{$tlen} : 0;
		}
	} else {
		$max_len = $rf_median * (1 + $len_ext);
		$min_len = $rf_median * (1 - $len_ext);
		foreach my $tlen (sort {$a <=> $b} keys %rf_tlens) {
			$cnt += ($tlen >= $min_len && $tlen <= $max_len) ? $rf_tlens{$tlen} : 0;
		}
	}

	$pass_cnt = $cnt;
	$pass_ratio = $pass_cnt / $final_ot_reads;
	if ($pass_ratio >= 0.9) { 
		last; 
	}
}

if ($fr > $rf) { printf "%.2f\t%.2f\t%.2f\t%.2f\t%.2f\tfr", $fr_median, $fr_mean, $fr_stdev, $min_len, $max_len; }
else { printf "%.2f\t%.2f\t\t%.2f\t%.2f\t%.2f\trf", $rf_median, $rf_mean, $rf_stdev, $min_len, $max_len; }

print STDERR "\n[$map_f]: Read filtering\n";
printf STDERR "\tMin. and max. fragment size: %.2f, %.2f\n", $min_len, $max_len;
printf STDERR "\tPassed reads: %d (%.2f)\n", $pass_cnt, $pass_ratio;

open(O, ">$new_f");
open(F, "$samtools view -F 0xF0C -q $min_mapq $map_f |");
while(<F>) {
	chomp;
	if ($_ =~ /^@/) { next; }
	my ($qname, $flag, $rname, $pos, $mapq, $cigar, $rnext, $pnext, $tlen) = split(/\s+/);
	if ($flag & 0x80 || $flag & 0x4 || $flag & 0x8) { next; }

	my ($dir1, $dir2) = ("+", "+");
	if ($flag & 0x10) { $dir1 = "-"; }
	if ($flag & 0x20) { $dir2 = "-"; }
	
	my ($rname1, $rname2) = ($rname, $rnext);
	if ($rname2 eq "=") { $rname2 = $rname1; }
	
	if ($rname1 ne $rname2) { next;	}
	if (!defined($scfnames{$rname1})) { next; }
    if ($pos == $pnext) { next;	}

	my $final_scf = $rname1;
	my ($final_start, $final_end) = (-1, -1);
	my $local_ot = "";

	if ($dir1 eq "+" && $dir2 eq "-") {
		if ($pos < $pnext) {
			$local_ot = "fr";
			($final_start, $final_end) = ($pos, $pos + abs($tlen) - 1);
		 } else {
			$local_ot = "rf";
			($final_start, $final_end) = ($pnext, $pnext + abs($tlen) - 1);
	     }
	} elsif ($dir1 eq "-" && $dir2 eq "+") {
		if ($pos > $pnext) {
			$local_ot = "fr";
			($final_start, $final_end) = ($pnext, $pnext + abs($tlen) - 1);
        } else {
			$local_ot = "rf";
			($final_start, $final_end) = ($pos, $pos + abs($tlen) - 1);
        }
	} 
	
	my $len = $final_end - $final_start;
    if ($final_ot ne $local_ot || $len > $max_len || $len < $min_len) { next; }
	
	print O "$final_scf\t$final_start\t$final_end\n"; 
}
close(F);
close(O);


##########
sub get_stats {
    my $cnt = shift;
    my $sum = shift;
    my $sum2 = shift;

    if ($cnt == 0) { return (0, 0); }

    my $mean = $sum/$cnt;
    my $var = $sum2/$cnt - ($mean * $mean);
    my $stdev = sqrt($var);

	return ($mean, $stdev);
}

sub cal_median {
	my %hs_tlens = %{$_[0]};
	my $cnt = $_[1];
	my $median = 0;

	my $idx_cnt = 0;
	my @tlens = sort {$a<=>$b} keys %hs_tlens;
	if ($cnt == 0) { return 0; }
	if ($cnt % 2 == 0) {
		my ($mid1_idx, $mid2_idx) = ($cnt/2, ($cnt+2)/2);
		for (my $i=0; $i<=$#tlens; $i++) {
			$idx_cnt += $hs_tlens{$tlens[$i]};
			if ($idx_cnt == $mid1_idx) {
				$median = ($tlens[$i] + $tlens[$i+1])/2;
				return $median;
			} elsif ($idx_cnt > $mid1_idx) {
				$median = $tlens[$i];
				return $median;
			}
		}
	} else {
		my $mid_idx = ($cnt + 1) / 2;
		for (my $i=0; $i<=$#tlens; $i++) {
			$idx_cnt += $hs_tlens{$tlens[$i]};
			if ($idx_cnt >= $mid_idx) {
				$median = $tlens[$i];
				return $median;
			}
		}
	}
}
