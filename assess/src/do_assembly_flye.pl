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
		"g=i" => \$genome_size
);
### Default options

print STDERR "# ---# Long read assembly with Flye\n";
my $time = `date`;
chomp($time);
print STDERR "# \t---# [START] $time\n";
print STDERR "# \t---# Final output will be stored at $out_dir/scaf.fa\n";
my @input_lib = @ARGV;
my %hs_lib = ();
my $nanopore_used = 0;
my $pacbio_used = 0;
my $hifi_used = 0;
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
						if($file =~ /.gz$/){ `gunzip -c $file >> $out_dir/nanopore.fa`; }
						else{ `cat $file >> $out_dir/nanopore.fa`; }
						$nanopore_used = 1;
				}
				elsif($lib_type eq "hifi"){
						if($file =~ /.gz$/){ `gunzip -c $file >> $out_dir/hifi.fa`; }
						else{ `cat $file >> $out_dir/hifi.fa`; }
						$hifi_used = 1;
				}
				else{ 
						if($file =~ /.gz$/){ `gunzip -c $file >> $out_dir/pacbio.fa`; }
						else{ `cat $file >> $out_dir/pacbio.fa`; }
						$pacbio_used = 1;
				}
		}
}
chdir("$out_dir");
my $cmd_option = "";
if($hifi_used == 1){ $cmd_option .= "--pacbio-hifi $out_dir/hifi.fa "; }
elsif($pacbio_used == 1){ $cmd_option .= "--pacbio-raw $out_dir/pacbio.fa "; }
elsif($nanopore_used == 1){ $cmd_option .= "--nano-raw $out_dir/nanopore.fa "; }
if(defined($genome_size)){ $cmd_option .= "-g $genome_size "; }
print "flye $cmd_option -t 1 --out-dir $out_dir/Flye.out\n";
`flye $cmd_option -t 1 --out-dir $out_dir/Flye.out > $out_dir/log.flye.txt 2>&1`;
`cp $out_dir/Flye.out/assembly.fasta $out_dir/scaf.fa`;

`rm -f $out_dir/ERRORCHECK`;
$time = `date`;
chomp($time);
print STDERR "# \t---# [ END ] $time\n";
print STDERR "# ---# Flye is finished (assembly out: $out_dir/scaf.fa)\n";
