#!/usr/bin/env perl
use strict;
use warnings;

my $genome_id = shift;
my $out_dir = shift;

print STDERR "# ---# Reference genome preparation \n";
my $time = `date`;
chomp($time);
print STDERR "# \t---# [START] $time\n";

print STDERR "# \t---# check if file exists\n";
if(!-f "$out_dir/$genome_id.fa.masked.gz"){
		print STDERR "# \t---# $genome_id sequence is not prepared\n";
		print STDERR "# \t---# downloading start\n";
		`wget https://hgdownload.soe.ucsc.edu/goldenPath/$genome_id/bigZips/$genome_id.2bit -O $out_dir/tmp.$genome_id.2bit > $out_dir/log.download.txt 2>&1`;
		`twoBitToFa $out_dir/tmp.$genome_id.2bit $out_dir/tmp.$genome_id.fa.masked`;
		`gzip $out_dir/tmp.$genome_id.fa.masked`;

		print STDERR "# \t---# converting unrecognized characters\n";
		open(FFA, "gunzip -c $out_dir/tmp.$genome_id.fa.masked.gz|");
		open(FOUT, ">$out_dir/$genome_id.fa.masked");
		my $id = "";
		while(<FFA>){
				chomp;
				if($_ =~ /^>(\S+)/){
						$id = $1;
						if($id !~ /_/){ print FOUT ">$id\n"; }
						else{ $id = ""; }
				}
				else{
						if($id ne ""){
								my $seq = $_;
								print FOUT "$seq\n";
						}
				}
		}
		close(FFA);
		close(FOUT);
		print STDERR "# \t---# compressing\n";
		`gzip $out_dir/$genome_id.fa.masked`;
		`rm -f $out_dir/tmp.$genome_id.fa.masked.gz`;
		$time = `date`;
		chomp($time);
		print STDERR "# \t---# $time\n";
		print STDERR "# ---# genome sequence downloading is finished ($out_dir/$genome_id.fa.masked.gz)\n";
}
else{
		$time = `date`;
		chomp($time);
		print STDERR "# \t---# [ END ] $time\n";
		print STDERR "# ---# $genome_id sequence is already prepared (skipped)\n";
}
