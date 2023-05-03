#!/usr/bin/env perl
use strict;
use warnings;

my $out_dir = shift;
`mkdir -p $out_dir`;
`touch $out_dir/ERRORCHECK`;

### Default options
my $max_rd_len = 150;
my $avg_ins = 300;
my $reverse_seq = 0;
my $asm_flags = 3;
my $rd_len_cutoff = 100;
my $rank = 1;
my $pair_num_cutoff = 3;
my $map_len = 32;
my $kmer = 31;

print STDERR "# ---# Short read assembly with SOAPdenovo2\n";
my $time = `date`;
chomp($time);
print STDERR "# \t---# [START] $time\n";
print STDERR "# \t---# Final output will be stored at $out_dir/scaf.fa\n";
my @input_lib = @ARGV;
my %hs_lib = ();
my $total_mean_ins = 0;
for(my $i = 0; $i <= $#input_lib; $i++){
		chomp($input_lib[$i]);
		my ($type, $num_lib, $file, $mean_ins, $stdev_ins) = split(/:/,$input_lib[$i]);
		$total_mean_ins += $mean_ins;
#		my ($type, $num_lib, $file) = split(/:/,$input_lib[$i]);
		if(!-f "$file"){
				print STDERR "# \t---# [ERROR] there is no sequence file ($file).\n";
#				exit();
		}
		my ($n_lib, $n_fq) = split(/-/,$num_lib);
		$hs_lib{$n_lib}{$n_fq} = $file;
}

$avg_ins = $total_mean_ins/(scalar(@input_lib));
print STDERR "# \t---# making configFile\n";
open(FOUT, ">$out_dir/configFile");
print FOUT "max_rd_len=$max_rd_len\n";
print FOUT "[LIB]\n";
print FOUT "avg_ins=$avg_ins\n";
print FOUT "reverse_seq=$reverse_seq\n";
print FOUT "asm_flags=$asm_flags\n";
print FOUT "rd_len_cutoff=$rd_len_cutoff\n";
print FOUT "rank=$rank\n";
print FOUT "pair_num_cutoff=$pair_num_cutoff\n";
print FOUT "map_len=$map_len\n";

print FOUT "#paired-end fastq files, read 1 file should always be followed by read 2 file\n";
foreach my $n_lib (sort {$a<=>$b} keys(%hs_lib)){
		print FOUT "q$n_lib=$hs_lib{$n_lib}{1}\n";
		print FOUT "q$n_lib=$hs_lib{$n_lib}{2}\n";
}
close(FOUT);
print STDERR "# \t---# configFile is prepared now\n";
print STDERR "# \t---# Running SPAPdenovo2\n";
`SOAPdenovo-63mer all -s $out_dir/configFile -o $out_dir/outputGraph -K $kmer > $out_dir/log.running_soapdenovo2.txt 2>&1`;
`cp $out_dir/outputGraph.scafSeq $out_dir/scaf.fa`;

`rm -f $out_dir/ERRORCHECK`;
$time = `date`;
chomp($time);
print STDERR "# \t---# [ END ] $time\n";
print STDERR "# ---# SOAPdenovo2 is finished (assembly out: $out_dir/scaf.fa)\n";
