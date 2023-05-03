#!/usr/bin/env perl
use strict;
use warnings;

my $out_dir = shift;
`mkdir -p $out_dir`;
`touch $out_dir/ERRORCHECK`;

### Default options

print STDERR "# ---# Short read assembly with masurca\n";
my $time = `date`;
chomp($time);
print STDERR "# \t---# [START] $time\n";
print STDERR "# \t---# Final output will be stored at $out_dir/scaf.fa\n";
print STDERR "# \t---# making configFile\n";
open(FOUT, ">$out_dir/configFile");
print FOUT "DATA\n";
my @input_lib = @ARGV;
my %hs_lib = ();
my $total_mean_ins = 0;
my $nanopore_used = 0;
my $pacbio_used = 0;
for(my $i = 0; $i <= $#input_lib; $i++){
		chomp($input_lib[$i]);
		my ($type, $num_lib, $file, $mean_ins, $stdev_ins) = split(/:/,$input_lib[$i]);
		if(!-f "$file"){
				print STDERR "# \t---# [ERROR] there is no sequence file ($file).\n";
#				exit();
		}
		if(defined($stdev_ins)){
				my ($n_lib, $n_fq) = split(/-/,$num_lib);
				$hs_lib{$n_lib}{$n_fq} = $file;
				$hs_lib{$n_lib}{'mean'} = $mean_ins;
				$hs_lib{$n_lib}{'stdev'} = $stdev_ins;
		}
		else{
				if($mean_ins eq "nanopore"){
						if($file =~ /.gz$/){ `gunzip -c $file >> $out_dir/nanopore.fa`; }
						else{ `cat $file >> $out_dir/nanopore.fa`; }
						$nanopore_used = 1;
				}
				else{ 
						if($file =~ /.gz$/){ `gunzip -c $file >> $out_dir/pacbio.fa`; }
						else{ `cat $file >> $out_dir/pacbio.fa`; }
						$pacbio_used = 1;
				}
		}
}
my @arr_letter = ("a" .. "z");
my ($m, $n) = (0, 0);
foreach my $n_lib (sort {$a<=>$b} keys(%hs_lib)){
		if($n == $#arr_letter){
				$m += 1;
				$n = 0;
		}
		my $prefix = "$arr_letter[$m]$arr_letter[$n]";
		print FOUT "PE= $prefix $hs_lib{$n_lib}{mean} $hs_lib{$n_lib}{stdev} $hs_lib{$n_lib}{1} $hs_lib{$n_lib}{2}\n";
		$n++;
}
if($nanopore_used == 1){ print FOUT "NANOPORE=$out_dir/nanopore.fa\n"; }
if($pacbio_used == 1){ print FOUT "PACBIO=$out_dir/pacbio.fa\n"; }
print FOUT "END\n";

print FOUT "PARAMETERS\n";
print FOUT "EXTEND_JUMP_READS=0\n";
print FOUT "GRAPH_KMER_SIZE = auto\n";
print FOUT "USE_LINKING_MATES = 0\n";
print FOUT "USE_GRID=0\n";
print FOUT "GRID_ENGINE=SGE\n";
print FOUT "GRID_QUEUE=all.q\n";
print FOUT "GRID_BATCH_SIZE=500000000\n";
print FOUT "LHE_COVERAGE=25\n";
print FOUT "MEGA_READS_ONE_PASS=0\n";
print FOUT "LIMIT_JUMP_COVERAGE = 300\n";
print FOUT "CA_PARAMETERS =  cgwErrorRate=0.15\n";
print FOUT "CLOSE_GAPS=1\n";
print FOUT "NUM_THREADS = 16\n";
print FOUT "JF_SIZE = 200000000\n";
print FOUT "SOAP_ASSEMBLY=0\n";
print FOUT "END\n";
close(FOUT);
print STDERR "# \t---# configFile is prepared now\n";
print STDERR "# \t---# Running Masurca\n";
chdir("$out_dir");
`masurca $out_dir/configFile`;
`bash $out_dir/assemble.sh`;
`cp $out_dir/CA/final.genome.scf.fasta $out_dir/scaf.fa`;

`rm -f $out_dir/ERRORCHECK`;
$time = `date`;
chomp($time);
print STDERR "# \t---# [ END ] $time\n";
print STDERR "# ---# Masurca is finished (assembly out: $out_dir/scaf.fa)\n";
