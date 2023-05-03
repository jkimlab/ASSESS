#!/usr/bin/perl

use strict;
use warnings;
use FindBin qw($Bin);

my $resolution = shift;
my $rec_spc = shift;
my $tar_spc = shift;
my $ref_spc = shift;
my $src_dir = shift;
my $src_f = shift;
my $out_dir = shift;
my $fa_f = shift;
my $ef = shift;

my $out_f = "$src_dir/outgroup.txt";
#my $block_f = "$src_dir/block_list.txt";
my $block_f = ($ef)? "$src_dir/../fill/block_list.update_scf.txt" : "$src_dir/block_list.txt";
my $seg_f = "$src_dir/Conserved.Segments";
my $ingroup_f = "$src_dir/ingroup.txt";

open(F,"$out_f");
my @outspcs = <F>;
close(F);
chomp(@outspcs);
system("mkdir -p $out_dir");
system("mkdir -p $out_dir/tmp");

my %ingroup = ();
open(F, "$ingroup_f");
while(<F>) {
	chomp;
	$ingroup{$_} = 1;
}
close(F);

my $rec_chrs_f = "$out_dir/rec_chrs.txt";
`$Bin/create_chrs.pl $block_f $fa_f $src_f > $rec_chrs_f`;

my $tar_seg_f = "$out_dir/rec_chrs.$tar_spc.segments.txt";
`$Bin/create_mapping.pl $rec_spc $tar_spc $rec_chrs_f > $tar_seg_f`;
my $tar_seg_refined_f = "$out_dir/rec_chrs.$tar_spc.segments.refined.txt";
my $tar_seg_refined_wgaps_f = "$out_dir/rec_chrs.$tar_spc.segments.refined.wgaps.txt";
`$Bin/merge_pos2ex.pl $rec_spc $tar_spc $tar_seg_f > $tar_seg_refined_f`;
`$Bin/merge_pos2ex.wgaps.pl $resolution $rec_spc $tar_spc $tar_seg_f > $tar_seg_refined_wgaps_f`;

foreach my $ing_spc (sort keys %ingroup) {
	if ($ing_spc eq $tar_spc) { next; }
	
	my $ing_seg_f = "$out_dir/rec_chrs.$ing_spc.segments.txt";
	`$Bin/infer_mapping_ing.pl $rec_spc $tar_spc $ing_spc $tar_seg_f $seg_f > $ing_seg_f`;
	my $ing_seg_refined_f = "$out_dir/rec_chrs.$ing_spc.segments.refined.txt";
	my $ing_seg_refined_wgaps_f = "$out_dir/rec_chrs.$ing_spc.segments.refined.wgaps.txt";
	`$Bin/merge_pos2ex.pl $rec_spc $ing_spc $ing_seg_f > $ing_seg_refined_f`;
	`$Bin/merge_pos2ex.wgaps.pl $resolution $rec_spc $ing_spc $ing_seg_f > $ing_seg_refined_wgaps_f`;
}

my $rec_chrs_refined_f = "$out_dir/rec_chrs.refined.txt";
`$Bin/merge_rec.pl $rec_chrs_f > $rec_chrs_refined_f`;

`$Bin/get_recon_size.pl $rec_chrs_refined_f > $out_dir/rec_chrs.size.txt`;

#`$Bin/print_adjscores.pl $out_dir/rec_chrs.txt $out_dir/raca2.out > $out_dir/rec_chrs.adjscores.txt`;
`$Bin/print_adjscores.pl $out_dir/rec_chrs.txt $src_f > $out_dir/rec_chrs.adjscores.txt`;

if($ef){
	`$Bin/create_fasta.block.pl $rec_chrs_refined_f $fa_f $out_dir/rec_chrs.fa > $out_dir/rec_chrs_coverage.txt`; 
}else{
	`$Bin/create_fasta.block.pl $rec_chrs_refined_f $fa_f $out_dir/rec_chrs.fa > $out_dir/rec_chrs_coverage.txt`; 
}

`mv $rec_chrs_f $out_dir/tmp/`;
foreach my $ing_spc (sort keys %ingroup) {
	my $ing_seg_f = "$out_dir/rec_chrs.$ing_spc.segments.txt";
	`mv $ing_seg_f $out_dir/tmp/`;	
}
