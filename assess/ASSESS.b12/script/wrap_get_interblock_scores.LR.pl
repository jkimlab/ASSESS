#!/usr/bin/perl

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../lib/perl";
use File::Basename;
use List::Util qw (min max);
use Parallel::ForkManager;

my $blist_f = shift;
my $blist_ext_f = shift;
my $data_dir = shift;
my $out_f = shift;
my $out_dir = shift;
my $numt = shift;
my $min_mapq = shift;

my $tmp_dir = "$out_dir/lr_tmp";
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
	`$Bin/get_interblock_scores.LR.pl $blist_f $blist_ext_f $file $tmp_dir/inter$pi.txt $min_mapq > $tmp_dir/$bname.log 2>&1`;
    $pm->finish;
}
$pm->wait_all_children;

my %sf_info = ();
my %mappings = ();

for (my $pi = 0; $pi < $numfiles; $pi++) {
	my $file = $files[$pi];
	if (!(-f $file)) { next; }
	if ($file !~ /\.sam$/ && $file !~ /\.bam$/) { next; }
	
	open(F, "$tmp_dir/inter$pi.txt");
	while(<F>) {
		chomp;
		my ($crd1, $dir1, $crd2, $dir2, $cnt, $bid1, $bid2) = split(/\s+/);

		$sf_info{$bid1} = "$crd1 $dir1";
		$sf_info{$bid2} = "$crd2 $dir2";

		if (defined($mappings{$bid1}{$bid2})) {
			$mappings{$bid1}{$bid2} += $cnt;
		} else {
			$mappings{$bid1}{$bid2} = $cnt;
		}
	}
	close(F);
}

# print results
open(O,">$out_f");
foreach my $bid1 (sort keys %mappings) {
	my $crd1 = $sf_info{$bid1};

	foreach my $bid2 (sort keys %{$mappings{$bid1}}) {
		my $cnt = $mappings{$bid1}{$bid2};
		my $crd2 = $sf_info{$bid2};
		
		print O "$crd1\t$crd2\t$cnt\t$bid1\t$bid2\n";
	}
}


