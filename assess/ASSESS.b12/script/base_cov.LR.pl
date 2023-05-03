#!/usr/bin/perl

use warnings;
use strict;

my $in_d = shift;
my $out_dir = shift;

my %depth = ();
my @depth_fs = <$in_d/*.depth.txt>;

#open(WF1, "> $out_dir/avgcov_scfs.txt");
if($#depth_fs >= 1){
	foreach my $fle (@depth_fs){
		open(RF, $fle);
		while(<RF>){
			chomp;

			my @col = split(/\s+/);
			if( exists $depth{$col[0]}{$col[1]} ){
				$depth{$col[0]}{$col[1]} += $col[2];
			}else{
				$depth{$col[0]}{$col[1]} = $col[2];
			}
		}
		close RF;
	}

	foreach my $chr (sort keys %depth){
		my $chr_size = 0;
		my $depth_sum = 0;

		open(WF2, ">$out_dir/$chr.cov");
		foreach my $base (sort {$a <=> $b} keys %{$depth{$chr}}){
			if($depth{$chr}{$base} != 0){
				++$chr_size;
				$depth_sum += $depth{$chr}{$base};
			}

			print "$chr\t$base\t$depth{$chr}{$base}\n";
			print WF2 "$base\t$depth{$chr}{$base}\n";
		}
		close WF2;

#		print WF1 "$chr\t".($depth_sum/$chr_size)."\n";
	}
}else{
	my $chr = "";
	my $line_stdout = "";
	my $line_file = "";
	my $chr_size = 0;
	my $depth_sum = 0;

	open ( RF, $depth_fs[0] );
	while( <RF> ){
		chomp;

		my @col = split(/\s+/);
		if($chr ne $col[0]){
			if($chr ne ""){	
#				print WF1 "$chr\t".($depth_sum / $chr_size)."\n";
				$chr_size = 0;
				$depth_sum = 0;

				open(WF2, ">$out_dir/$chr.cov");
				print WF2 $line_file;
				print $line_stdout;
				close WF2;
			}

			$chr = $col[0];
			$line_file = "";
			$line_stdout = "";
		}
		
		if($col[2] != 0){
			++$chr_size;
			$depth_sum += $col[2];
		}

		$line_file .= "$col[1]\t$col[2]\n";
		$line_stdout .= "$_\n";
	}
#	print WF1 "$chr\t".($depth_sum / $chr_size)."\n";
	
	open(WF2, ">$out_dir/$chr.cov");
	print WF2 $line_file;
	print $line_stdout;
	close WF2;
}
#close WF1;
