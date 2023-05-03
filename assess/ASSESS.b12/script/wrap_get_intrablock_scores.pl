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
my $avgcov_f = shift;
my $out_dir = shift;

my $tmp_dir = "$out_dir/r_tmp";
`mkdir -p $tmp_dir`;

my %hs_blocks = ();
my %hs_blockschr = ();
my %hs_blocks_map = ();
my %hs_blocks_dir = ();
open(F,"$block_f");
while(<F>) {
	chomp;
	my ($scf, $start, $end, $dir, $bid, $rchr) = split(/\s+/);
	$hs_blocks{$scf}{$start} = $end;
	$hs_blockschr{$scf}{$start}{$end} = $rchr;
	$hs_blocks_map{"$scf:$start-$end"} = $bid;
	$hs_blocks_dir{$bid} = $dir;
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

    #print STDERR "$sstart\t$send\t$scnt\n";
    $pm->start and next;
    `$Bin/get_intrablock_scores.range.pl $sstart $send $block_f $flanksize $windowsize $cov_dir $avgcov_f $tmp_dir/intra$pi.txt`;
    $pm->finish;
}
$pm->wait_all_children;

my @finals = ();

for (my $pi = 0; $pi < $num_p; $pi++) {
    my $sstart = $cbin * $pi;
    if ($sstart > $scnt - 1) { last; }

    my $f = "$tmp_dir/intra$pi.txt";
    if (!(-f $f)) { next; }

    open(F, "$f");
    while(<F>) {
        chomp;
        push(@finals, $_);
    }
    close(F);
}

#`rm -rf $tmp_dir`;

open(O,">$out_f");
foreach my $out (@finals) {
	my ($crd1, $dir1, $crd2, $dir2) = split(/\s+/, $out);
	my $bid1 = $hs_blocks_map{$crd1};
    my $orgdir1 = $hs_blocks_dir{$bid1};
    if ($orgdir1 ne $dir1) { $bid1 *= -1; }

	my $bid2 = $hs_blocks_map{$crd2};
	my $orgdir2 = $hs_blocks_dir{$bid2};
	if ($orgdir2 ne $dir2) { $bid2 *= -1; }

	print O "$out\t$bid1\t$bid2\n";
}
close(O);
