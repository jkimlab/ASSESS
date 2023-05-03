#!/usr/bin/perl
use strict;
use warnings;
use Cwd 'abs_path';
use FindBin '$Bin';
use File::Basename;


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

`mkdir -p $out_dir/data $out_dir/fastqs $out_dir/bam`;
`cp $ref_fa $out_dir/data/reference.fa`;
$ref_fa = "$out_dir/data/reference.fa";
if($f_fa =~ /.gz$/){ 
		print "$f_fa is a compressed file --> copying\n";
		`cp $f_fa $out_dir/fastqs/$prefix\_S1_L001_R1_001.fastq.gz`;
		$f_fa = "$out_dir/fastqs/$prefix\_S1_L001_R1_001.fastq.gz";
}
else{
		print "$f_fa is not a compressed file --> compressing\n";
		`gzip -c $f_fa > $out_dir/fastqs/$prefix\_S1_L001_R1_001.fastq.gz`;
		$f_fa = "$out_dir/fastqs/$prefix\_S1_L001_R1_001.fastq.gz";
}
if($r_fa =~ /.gz$/){
		print "$r_fa is a compressed file --> copying\n";
		`cp $r_fa $out_dir/fastqs/$prefix\_S1_L001_R2_001.fastq.gz`;
		$r_fa = "$out_dir/fastqs/$prefix\_S1_L001_R2_001.fastq.gz";
}
else{ 
		print "$r_fa is not a compressed file --> compressing\n";
		`gzip -c $r_fa > $out_dir/fastqs/$prefix\_S1_L001_R2_001.fastq.gz`;
		$r_fa = "$out_dir/fastqs/$prefix\_S1_L001_R2_001.fastq.gz";
}
print STDERR "# \t---# [CMD] samtools faidx $ref_fa -o $ref_fa.fai\n";
`samtools faidx $ref_fa -o $ref_fa.fai`;
my $ref_name = basename("$ref_fa");
print STDERR "# ---# Read mapping to initial assembly with LongRanger mapping pipeline\n";
my $time = `date`;
chomp($time);
print STDERR "# \t---# [START] $time\n";

print STDERR "# \t---# mapping start (LongRanger pipeline)\n";
print STDERR "# \t---# make reference\n";
chdir("$out_dir");
print STDERR "# \t---# [CMD] longranger mkref ./data/reference.fa\n";
`longranger mkref ./data/reference.fa 2> log.mkref.txt 2>&1`;
print STDERR "# \t---# longranger align\n";
print STDERR "# \t---# [CMD] longranger align --fastqs=./fastqs --id=$prefix --reference=./refdata-reference --sample=$prefix\n";
`longranger align --fastqs=./fastqs --id=$prefix --reference=./refdata-reference --sample=$prefix > log.align.txt 2>&1`;
`cp $prefix/outs/*.bam $out_dir/bam/$prefix.bam`;
