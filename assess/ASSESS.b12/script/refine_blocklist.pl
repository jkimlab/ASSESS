#!/usr/bin/perl

use strict;
use warnings;

my $block_f = shift;
my $block_ext_f = shift;
my $fasize_f = shift;
my $out_dir = shift;

`mkdir -p $out_dir`;
my $scf_block_f = "$out_dir/block_list.update_scf.txt";
my $scf_block_ext_f = "$out_dir/block_list.update_scf.ext.txt";

my %built_sf = ();
my $last_bid = 0;

open (WF1, ">$scf_block_f");
open (WF2, ">$scf_block_ext_f");

open (RF, $block_f);
while(<RF>){
	chomp;
	print WF1 "$_\n";

	my ( $scf, $scf_beg, $scf_end, $dir, $bid ) = split(/\s+/);
	$built_sf{$scf} = 1;

	if($bid > $last_bid){ $last_bid = $bid; }
}
close RF;

open (RF, $block_ext_f);
while(<RF>){
	chomp;
	print WF2 "$_\n";
}
close RF;


open (RF, $fasize_f);
while(<RF>){
	chomp;
	
	my ( $scf, $len ) = split(/\s+/);
	if(! exists $built_sf{$scf}){
		print WF1 "$scf\t0\t$len\t+\t".++$last_bid."\tNA\n";	
		print WF2 "$scf\t0\t$len\t+\t".$last_bid."\n";	
	}
}
close RF;

close WF1;
close WF2;

