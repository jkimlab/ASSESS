#!/usr/bin/perl

use strict;
use warnings;
use POSIX;
use FindBin qw($Bin);
use lib "$Bin/../lib/perl";
use Parallel::ForkManager;

my $block_f = shift;
my $flanksize = shift;
my $windowsize = shift;
my $cov_dir = shift;
my $num_p = shift;
my $out_f = shift;
my $out_avgcov_f = shift;
my $out_dir = shift; 

my $tmp_dir = "$out_dir/r_tmp";
`mkdir -p $tmp_dir`;

my %hs_blocks = ();
open(F,"$block_f");
while(<F>) {
	chomp;
	my ($scf, $start, $end) = split(/\s+/);
	$hs_blocks{$scf}{$start} = $end;
}
close(F);

my @scfs = sort keys %hs_blocks;
my $scnt = scalar(@scfs);
my $cbin = ceil($scnt/$num_p);
if ($cbin == 0) { $cbin = 1; }

my $pm = new Parallel::ForkManager($num_p);
for (my $pi = 0; $pi < $num_p; $pi++) {
	my $sstart = $cbin * $pi;
	if ($sstart > $scnt - 1) { last; }

	my $send = $sstart + $cbin - 1;
	if ($send > $scnt-1) { $send = $scnt-1; }
    if ($pi == $num_p - 1) { $send = $scnt-1; }

	$pm->start and next;
	`$Bin/estimate_intrablock_thr.range.pl $sstart $send $block_f $flanksize $windowsize $cov_dir $tmp_dir/thr$pi.txt`; 
	$pm->finish;
}
$pm->wait_all_children;

my @pct_values = ();
my @raw_values = ();
my @weights = ();
my %hs_avgcov = ();

my ($pct_sum, $raw_sum) = (0, 0);
my ($pct_sum2, $raw_sum2) = (0, 0);

for (my $pi = 0; $pi < $num_p; $pi++) {
	my $sstart = $cbin * $pi;
	if ($sstart > $scnt - 1) { last; }

    my $f = "$tmp_dir/thr$pi.txt";
    if (!(-f $f)) { next; }

	open(F, "$f");
	while(<F>) {
		chomp;
		if ($_ =~ /^#/) {
			my @ar = split(/\s+/);
			$hs_avgcov{$ar[1]} = $ar[2];
		} else {
			my ($pct, $raw) = split(/\s+/);
			push(@pct_values, $pct);
			push(@raw_values, $raw);
			push(@weights, 1);

            $pct_sum += $pct;
            $raw_sum += $raw;
            $pct_sum2 += ($pct * $pct);
            $raw_sum2 += ($raw * $raw);
		}
	}
	close(F);
}

`rm -f $tmp_dir/thr*.txt`;

@pct_values = sort {$a<=>$b} @pct_values;
@raw_values = sort {$a<=>$b} @raw_values;
my $num_values = scalar(@pct_values);

my $pct_sd = calc_sd($num_values, $pct_sum, $pct_sum2); 
my $raw_sd = calc_sd($num_values, $raw_sum, $raw_sum2); 

my $val01 = percentile(0.01, $num_values, \@pct_values);
my $val05 = percentile(0.05, $num_values, \@pct_values);
print "$val05\t$val01\n";

# Output to a file
open(O,">$out_f");
print O "# pct_sd: $pct_sd\traw_sd: $raw_sd\n";
for (my $pthr = 1; $pthr <= 999999; $pthr++) {
	my $pct = percentile($pthr/1000000, $num_values, \@pct_values);
	my $raw = percentile($pthr/1000000, $num_values, \@raw_values);
	print O $pthr/10000, "\t$pct\t$raw\n";
}
close(O);

open(O, ">$out_avgcov_f");
foreach my $scf (sort keys %hs_avgcov) {
	my $avgcov = $hs_avgcov{$scf};
	print O "$scf\t$avgcov\n";
}
close(O);

sub calc_sd {
    my $num = shift;
    my $sum = shift;
    my $sum2 = shift;

    my $mean = $sum/$num;
    my $var = $num/($num-1) * ($sum2/$num - $mean * $mean);
    return sqrt($var);
}

sub percentile {
	my $p = shift;
	my $total = shift;
	my $aref = shift;
	
	my $value = "";
	my $n = $p * $total;

	if ($n =~ m/^[d]*$/) { # integer
		$value = ($$aref[$n-1] + $$aref[$n])/2;
	} else {	# non-integer
		my $i = floor($n) + 1;
		$i--;
		$value = $$aref[$i];
	}

	return $value;
}
