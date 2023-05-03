#!/usr/bin/perl

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../lib/perl";
use Parallel::ForkManager;
use File::Basename;

my $samtools = "$Bin/../util/samtools-1.5/samtools";
my $bedtoolsbin = "$Bin/../util/bedtools2/bin";

my $src_dir = shift;		# input directory containing SAM or BAM files
my $num_p = shift;			# number of threads
my $out_dir = shift;		# output directory
my $cov_dir = shift;		# coverage file directory
my $size_f = shift;			# assembly size file
my $stat_f = shift;
my $min_mapq = shift;


`mkdir -p $out_dir`;
`mkdir -p $cov_dir`;

my @mapfiles = <$src_dir/*>;
my $numfiles = scalar(@mapfiles);

print STDERR "Processing read mapping files...\n";
open(O, ">$stat_f");
my $pm = new Parallel::ForkManager($num_p);
for (my $pi = 0; $pi < $numfiles; $pi++) {
	my $mapfile = $mapfiles[$pi];
	if (!(-f $mapfile)) { next; }
	if ($mapfile !~ /\.(sam|SAM)$/ && $mapfile !~ /\.(bam|BAM)$/) { next; }
	my $bname = basename($mapfile);
	$pm->start and next;
	
	print STDERR "[P$pi Checking insert size and read orientation] $bname ...\n";
	my $outstr = `$Bin/check_read_lenot.pl $mapfile $size_f $out_dir/$bname.nsorted.checked.bed $min_mapq`;
	chomp($outstr);
	my ($median, $mean, $stdev, $min_len, $max_len, $read_ot) = split(/\s+/, $outstr);
	print O "$bname\t$median\t$mean($stdev)\t$min_len\t$max_len\t$read_ot\n";

	$pm->finish;
}
$pm->wait_all_children;
close(O);

print STDERR "Merging and sorting bed files...\n";
`cat $out_dir/*.nsorted.checked.bed | sort -k1,1 -k2n,2 -k3n,3 -u > $out_dir/all.bed`;
`rm -f $out_dir/*.nsorted.checked.bed`;

print STDERR "Computing physical coverage...\n";
`$bedtoolsbin/genomeCoverageBed -i $out_dir/all.bed -d -g $size_f > $out_dir/all_cov.txt`;
`rm -f $out_dir/all.bed`;

open(F, "$out_dir/all_cov.txt");
my $pname = "";
my $fh;
while(<F>) {
    chomp;
    my ($name, $pos, $cov) = split(/\s+/);
    if ($name ne $pname) {
        if (defined($fh)) { close($fh); }
        open($fh, ">$cov_dir/$name.cov");
    }
    print $fh "$pos\t$cov\n";
    $pname = $name;
}
close($fh);
close(F);
`rm -f $out_dir/all_cov.txt`;
print STDERR "Done.\n";

