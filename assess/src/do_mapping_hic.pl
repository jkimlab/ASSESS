#!/usr/bin/perl
use strict;
use warnings;
use Cwd 'abs_path';
use FindBin '$Bin';


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

`mkdir -p $out_dir/data`;
`cp $ref_fa $out_dir/data/reference.fa`;
$ref_fa = "$out_dir/data/reference.fa";
if($f_fa =~ /.gz$/){ 
		`cp $f_fa $out_dir/data/$prefix\_1.fastq.gz`;
		$f_fa = "$out_dir/data/$prefix\_1.fastq.gz";
}
else{
		`gzip -c $f_fa > $out_dir/data/$prefix\_1.fastq.gz`;
		$f_fa = "$out_dir/data/$prefix\_1.fastq.gz";
}
if($r_fa =~ /.gz$/){
		`cp $r_fa $out_dir/data/$prefix\_2.fastq.gz`;
		$r_fa = "$out_dir/data/$prefix\_2.fastq.gz";
}
else{ 
		`gzip -c $r_fa > $out_dir/data/$prefix\_2.fastq.gz`;
		$r_fa = "$out_dir/data/$prefix\_2.fastq.gz";
}
`samtools faidx $ref_fa -o $ref_fa.fai`;
open(FSH, ">$out_dir/01_mapping_arima.sh");

print STDERR "# ---# Read mapping to initial assembly with Arima mapping pipeline\n";
my $time = `date`;
chomp($time);
print STDERR "# \t---# [START] $time\n";

print FSH "#! /bin/bash\n";
print FSH "SRA=\'"."$prefix"."\'\n";
print FSH "LABEL=\'"."$prefix"."\'\n";
print FSH "REP_LABEL=\$LABEL\_rep1\n";
print FSH "IN_DIR=\'"."$out_dir/data"."\'\n";
print FSH "REF=\'"."$ref_fa"."\'\n";
print FSH "FAIDX=\'"."$ref_fa.fai"."\'\n";
print FSH "PREFIX=\'"."$ref_fa"."\'\n";
print FSH "RAW_DIR=\'"."$out_dir/bam/raw"."\'\n";
print FSH "FILT_DIR=\'"."$out_dir/bam/filt"."\'\n";
print FSH "TMP_DIR=\'"."$out_dir/tmp"."\'\n";
print FSH "PAIR_DIR=\'"."$out_dir/bam/pair"."\'\n";
print FSH "REP_DIR=\'"."$out_dir/bam/dup"."\'\n";
print FSH "MERGE_DIR=\'"."$out_dir/bam/merge"."\'\n";
close(FSH);
`cat $Bin/../sources/arima_pipeline >> $out_dir/01_mapping_arima.sh`;

print STDERR "# \t---# mapping start (Arima pipeline)\n";
chdir("$out_dir");
`chmod +x ./01_mapping_arima.sh`;
`./01_mapping_arima.sh`;

`rm -f $out_dir/bam/filt/*.bam $out_dir/bam/raw/*.bam $out_dir/bam/merge/*.bam`;
=pod
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
=cut
