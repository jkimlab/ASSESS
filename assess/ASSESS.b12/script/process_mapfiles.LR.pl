#!/usr/bin/perl

use strict;
use warnings;
use FindBin qw($Bin);
use File::Basename;
use Parallel::ForkManager;

my $samtools = "$Bin/../util/samtools-1.5/samtools";

my $src_dir = shift;		# input directory containing SAM or BAM files
my $threads = shift;
my $out_dir = shift;		# output directory
my $cov_dir = shift;
my $min_mapq = shift;
my $seqsize_f = shift;

`mkdir -p $out_dir`;
`mkdir -p $cov_dir`;

my $tmp_bam_dir = "$out_dir/bam/";
`mkdir -p $tmp_bam_dir`;

print STDERR "Processing read mapping files...\n";
my @bams = <$src_dir/*>;
my @sort_bams = ();

foreach my $bam (@bams){
	if ($bam !~ /\.(sam|SAM)$/ && $bam !~ /\.(bam|BAM)$/) { next; }
	my $cmd = "$samtools view -H $bam";
	my $out = `$cmd`;

	my $is_sort = 0;
	if($out =~ /SO:coordinate/){
		$is_sort = 1;
	}

	if($is_sort){
		my $index_cmd = "$samtools index -@ $threads $bam";
		`$index_cmd`;
		push(@sort_bams, $bam);
	}else{
		my $prefix = basename($bam);
		my $sort_cmd = "$samtools sort -@ $threads $bam > $tmp_bam_dir/sort.$prefix";
		`$sort_cmd`;
		my $index_cmd = "$samtools index -@ $threads $tmp_bam_dir/sort.$prefix";
		`$index_cmd`;
		push(@sort_bams, "$tmp_bam_dir/sort.$prefix");
	}
}

my @seq_ns = ();
open(F, $seqsize_f);
while(<F>) {
	chomp;
	my ($seq_n, $size) = split(/\s+/);
	push(@seq_ns, $seq_n);
}
close(F);

my @split_bams = ();
foreach my $sorted_bam (@sort_bams) {
	foreach my $chr (@seq_ns) {
		my $prefix = basename($sorted_bam);
		my $split_cmd = "$samtools view -Sbh -@ $threads $sorted_bam $chr > $tmp_bam_dir/$chr.$prefix";
		print STDERR "$split_cmd\n";
		`$split_cmd`;
		push(@split_bams, "$tmp_bam_dir/$chr.$prefix");
	}
}

my $fm = new Parallel::ForkManager($threads);
foreach my $bam (@split_bams)
{
	if ($bam !~ /\.(sam|SAM)$/ && $bam !~ /\.(bam|BAM)$/) { next; }

	my $prefix = basename ($bam);
	my $cmd = "$samtools depth -Q $min_mapq -a $bam > $out_dir/$prefix.depth.txt";
	print STDERR "$cmd\n";
	my $fm_qid = $fm -> start($cmd) and next;
	system($cmd);
	$fm -> finish(0, \$cmd);
}
$fm -> wait_all_children;

`$Bin/base_cov.LR.pl $out_dir $cov_dir`; 
#`$Bin/estimate_depth_dist.pl $out_dir/LR.w$wsize.cov $out_dir/LR.w$wsize.dist $out_dir/LR.w$wsize.stat`;

print STDERR "Done.\n";
