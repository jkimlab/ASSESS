#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';

my $out_dir = shift;
`mkdir -p $out_dir`;
`touch $out_dir/ERRORCHECK`;

my $genome_size;
GetOptions (
		"g=s" => \$genome_size
);
### Default options

print STDERR "# ---# Long read assembly with Canu\n";
my $time = `date`;
chomp($time);
print STDERR "# \t---# [START] $time\n";
print STDERR "# \t---# Final output will be stored at $out_dir/scaf.fa\n";
my @input_lib = @ARGV;
my %hs_lib = ();
my $nanopore_used = 0;
my $pacbio_used = 0;
my $hifi_used = 0;
my $nanopore_files = "";
my $pacbio_files = "";
my $hifi_files = "";
for(my $i = 0; $i <= $#input_lib; $i++){
		chomp($input_lib[$i]);
		my ($type, $num_lib, $file, $lib_type, $tmp) = split(/:/,$input_lib[$i]);
		if(!-f "$file"){
				print STDERR "# \t---# [ERROR] there is no sequence file ($file).\n";
#				exit();
		}
		if($type ne "LR"){ next; }
		else{
				if($lib_type eq "nanopore"){
						$nanopore_files .= " $file";
						$nanopore_used = 1;
				}
				elsif($lib_type eq "hifi"){
						$hifi_files .= " $file";
						$hifi_used = 1;
				}
				else{ 
						$pacbio_files .= " $file";
						$pacbio_used = 1;
				}
		}
}
chdir("$out_dir");
my $cmd_option = "-p canuout -d $out_dir/canu.out genomeSize=$genome_size ";
if($hifi_used == 1){ $cmd_option .= "-pacbio-hifi$hifi_files"; }
elsif($pacbio_used == 1){ $cmd_option .= "-pacbio$pacbio_files"; }
elsif($nanopore_used == 1){ $cmd_option .= "-nanopore$nanopore_files"; }

print "/RACA-w/bin/canu/build/bin/canu $cmd_option useGrid=false\n";
`/RACA-w/bin/canu/build/bin/canu $cmd_option useGrid=false > $out_dir/log.canu.txt 2>&1`;
`cp $out_dir/canu.out/canuout.contigs.fasta $out_dir/scaf.fa`;

`rm -f $out_dir/ERRORCHECK`;
$time = `date`;
chomp($time);
print STDERR "# \t---# [ END ] $time\n";
print STDERR "# ---# Canu is finished (assembly out: $out_dir/scaf.fa)\n";
