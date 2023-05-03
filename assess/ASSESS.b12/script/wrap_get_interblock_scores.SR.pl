#!/usr/bin/perl

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../lib/perl";
use File::Basename;
use List::Util qw (min max);
use Parallel::ForkManager;

my $insert_f = shift;
my $blist_f = shift;
my $blist_ext_f = shift;
my $size_f = shift;
my $data_dir = shift;
my $out_f = shift;
my $out_dir = shift;
my $numt = shift;
my $min_mapq = shift;

my $tmp_dir = "$out_dir/sr_tmp";
`mkdir -p $tmp_dir`;

# process paired-ends mappings
my @files = <$data_dir/*>;
my $numfiles = scalar(@files);

my $pm = new Parallel::ForkManager($numt);
for (my $pi = 0; $pi < $numfiles; $pi++) {
	my $file = $files[$pi];
	my $bname = basename($file); 
	if (!(-f $file)) { next; }
	if ($file !~ /\.sam$/ && $file !~ /\.bam$/) { next; }

	$pm->start and next;
	`$Bin/get_interblock_scores.single.pl $insert_f $blist_f $blist_ext_f $size_f $file $tmp_dir/inter$pi.txt $min_mapq & $tmp_dir/$bname.log 2>&1`;
    $pm->finish;

}
$pm->wait_all_children;

my %mappings = ();

for (my $pi = 0; $pi < $numfiles; $pi++) {
	my $file = $files[$pi];
	if (!(-f $file)) { next; }
	if ($file !~ /\.sam$/ && $file !~ /\.bam$/) { next; }
	
	open(F, "$tmp_dir/inter$pi.txt");
	while(<F>) {
		chomp;
		my ($crd1, $dir1, $crd2, $dir2, $cnt, $bid1, $bid2) = split(/\s+/);

		my $key1 = "$crd1 $dir1";
		my $key2 = "$crd2 $dir2";
	
		if (defined($mappings{$key1}{$key2})) {
			$mappings{$key1}{$key2} += $cnt;
		} else {
			$mappings{$key1}{$key2} = $cnt;
		}
	}
	close(F);
}

#`rm -rf $tmp_dir`;

my %hs_blocks_map = ();
my %hs_blocks_dir = ();
open(F,"$blist_f");
while(<F>) {
    chomp;
    my ($scfnum, $start, $end, $dir, $bid) = split(/\s+/);
	$hs_blocks_map{"$scfnum:$start-$end"} = $bid;
	$hs_blocks_dir{$bid} = $dir;
}
close(F);

# print results
open(O,">$out_f");
foreach my $key1 (sort keys %mappings) {
	my $rhs = $mappings{$key1};
	my ($crd1, $dir1) = split(/\s+/, $key1);
	my $bid1 = $hs_blocks_map{$crd1};
	my $orgdir1 = $hs_blocks_dir{$bid1};
	if ($orgdir1 ne $dir1) { $bid1 *= -1; }

	foreach my $key2 (sort keys %$rhs) {
		my $cnt = $$rhs{$key2};
		my ($crd2, $dir2) = split(/\s+/, $key2);
		my $bid2 = $hs_blocks_map{$crd2};
		my $orgdir2 = $hs_blocks_dir{$bid2};
		if ($orgdir2 ne $dir2) { $bid2 *= -1; }
		print O "$key1\t$key2\t$cnt\t$bid1\t$bid2\n";
	}
}


