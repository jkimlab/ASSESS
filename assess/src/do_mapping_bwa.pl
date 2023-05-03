#!/usr/bin/perl
use strict;
use warnings;
use Cwd 'abs_path';

my $prefix = shift;
my $ref_fa = shift;
my $out_dir = shift;
my $f_fa = shift;
my $r_fa = shift;
my $threads = shift;
my $index_dir = shift;
=pod
if($threads >= 2){
		$threads = int($threads/2);
}
=cut

`mkdir -p $out_dir`;
print STDERR "# ---# Read mapping to initial assembly with BWA\n";
my $time = `date`;
chomp($time);
print STDERR "# \t---# [START] $time\n";

print STDERR "# \t---# check if index exists\n";
if(!-f "$index_dir/done.indexing"){
		print STDERR "# \t---# index file is not prepared, start indexing with initial assembly\n";
		`bwa index -p $index_dir/index $ref_fa`;
		`touch $index_dir/done.indexing`;
}
else{
		print STDERR "# \t---# index file is already prepared (skipped)\n";
}
print STDERR "# \t---# mapping start (BWA-MEM)\n";
`bwa mem -t $threads $index_dir/index $f_fa $r_fa | samtools view -Sb - > $out_dir/$prefix.bam`;
if(-e "$out_dir/$prefix.bam"){
		$time = `date`;
		chomp($time);
		print STDERR "# \t---# $time\n";
		print STDERR "# ---# BWA is finished ($out_dir/$prefix.bam)\n";
		`rm -f $out_dir/ERRORCHECK`;
}
else{
		$time = `date`;
		chomp($time);
		print STDERR "# \t---# [ END ] $time\n";
		print STDERR "# ---# [ERROR] BWA fail\n";
}

