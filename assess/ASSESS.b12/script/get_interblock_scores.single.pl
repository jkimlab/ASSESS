#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;
use List::Util qw (min max);
use FindBin qw($Bin);

my $insert_f = shift;
my $blist_f = shift;
my $blistext_f = shift;
my $size_f = shift;
my $data_f = shift;	# SAM or BAM file
my $out_f = shift;
my $min_mapq = shift;

my $len_ext = 0.25;
my $samtools = "$Bin/../util/samtools-1.5/samtools";

my %hs_readdir = ();
my %hs_maxsize = ();
my %insertsizes_cal = ();
open(F,"$insert_f");
while(<F>) {
    chomp;
	if (length($_) == 0 || $_ =~ /^#/) { next; }
	my ($insert, $calsize, $mean_sd, $min_len, $max_len, $readdir) = split(/\s+/);
	#my ($insert, $calsize, $mean_sd, $readdir) = split(/\s+/);
    $insertsizes_cal{$insert} = $calsize;
	$hs_readdir{$insert} = lc($readdir);
	$hs_maxsize{$insert} = $max_len;
}
close(F);

# read synteny blocks
my %hs_blocks = ();
my %hs_blocks_id = ();
open(F,"$blistext_f");
while(<F>) {
	chomp;
	my ($scfnum, $start, $end, $dir, $bid) = split(/\s+/);
	$hs_blocks{$scfnum}{$start} = $end;
	$hs_blocks_id{"$scfnum:$start-$end"} = $bid;
}
close(F);

my %hs_blocks_map = ();
open(F,"$blist_f");
while(<F>) {
	chomp;
	my ($scfnum, $start, $end, $dir, $bid) = split(/\s+/);
	$hs_blocks_map{$bid} = "$scfnum:$start-$end";
}
close(F);

my %scf_size = ();
open(F,"$size_f");
while(<F>) {
    chomp;
    my ($scf, $size) = split(/\s+/);
    $scf_size{$scf} = $size;
}
close(F);

my %mappings = ();

# process paired-ends mappings
my $insertname = basename($data_f);
my $insertreaddir = $hs_readdir{$insertname};
my $insertsize_cal = $insertsizes_cal{$insertname};
my $maxsize = $hs_maxsize{$insertname};
#my $maxsize = $insertsize_cal * (1+$len_ext);

my %unpaired = ();
open(F, "$samtools view -F 0xF0C -q $min_mapq $data_f |");
while(<F>) {
    chomp;
    if ($_ =~ /^@/) { next; }
    my ($readname, $flag) = split(/\s+/);

    if (defined($unpaired{$readname})) {
        delete $unpaired{$readname};
    } else {
        $unpaired{$readname} = 1;
    }
}
close(F);

my $upcnt = 0;
open(F, "$samtools view -F 0xF0C -q $min_mapq $data_f |");

my ($petotal, $peused) = (0, 0);

my %hs_used = ();
while(<F>) {
    chomp;
	if ($_ =~ /^@/) { next; } # skip header

	my @ar = split(/\s+/);
	my ($readname, $flag, $readseq) = ($ar[0], $ar[1], $ar[9]);
	my ($scfname1, $pos1, $mapq, $scfname2, $pos2) = ($ar[2], $ar[3], $ar[4], $ar[6], $ar[7]);

	# filtering conditions
	if ($scfname1 eq "*" || $scfname2 eq "*") { next; }
	if ($scfname1 eq "=" || $scfname2 eq "=") { next; }
	if ($scfname1 eq $scfname2) { next; }
	if ($flag & 0x80) { next; }

    if (defined($unpaired{$readname})) { $upcnt++; next; }

	my $readlen = length($readseq);

	my $readdir1 = "+";
	if ($flag & 0x10) { $readdir1 = "-"; }
	my $readdir2 = "+";
	if ($flag & 0x20) { $readdir2 = "-"; }

	# check the existency of scaffolds in the list
	if (!defined($hs_blocks{$scfname1} || !defined($hs_blocks{$scfname2}))) { 
		next; 
	} 	
	
	if (defined($hs_used{"$flag:$scfname1:$pos1:$scfname2:$pos2"})) { next; }
	$hs_used{"$flag:$scfname1:$pos1:$scfname2:$pos2"} = 1;

	my ($scf1len, $scf2len) = ($scf_size{$scfname1}, $scf_size{$scfname2});
	my ($pestart1, $peend1) = ($pos1, $pos1+$readlen-1);
	my ($pestart2, $peend2) = ($pos2, $pos2+$readlen-1);
	
	# search for the first read	
	my ($block_start1, $block_end1) = read_search($pestart1, $peend1, $scfname1, \%hs_blocks);
	if ($block_start1 == -1) { next; }	
	
	# search for the second read	
	my ($block_start2, $block_end2) = read_search($pestart2, $peend2, $scfname2, \%hs_blocks);
	if ($block_start2 == -1) { next; }

	# check distance between reads
	my $dist = 0;		

		if ($insertreaddir eq "fr") {
			if ($readdir1 eq "+") { $dist += ($scf1len - $pos1 + 1); }
			else { $dist += ($pos1 + $readlen - 1); }
			
			if ($readdir2 eq "-") { $dist += ($pos2 + $readlen - 1); }
			else { $dist += ($scf2len - $pos2 + 1); }
		} elsif ($insertreaddir eq "rf") {
			if ($readdir1 eq "-") { $dist += ($scf1len - $pos1 + 1); }
			else { $dist += ($pos1 + $readlen - 1); }
				
			if ($readdir2 eq "+") { $dist += ($pos2 + $readlen - 1); }
			else { $dist += ($scf2len - $pos2 + 1); }
		} else {
			print STDERR "\nError: read orientation $insertreaddir\n";
			die;
		}

	if ($dist > $maxsize) { next; }

	# check block directions
	my $key1 = "$scfname1:$block_start1-$block_end1 ";
	my $key2 = "$scfname2:$block_start2-$block_end2 ";
	my $keydir1 = "+";
	my $keydir2 = "+";
	
		if ($insertreaddir eq "fr") {
			if ($readdir1 eq "-") { $keydir1 = "-"; }
			if ($readdir2 eq "+") { $keydir2 = "-"; }
		} elsif ($insertreaddir eq "rf") {	# org: -/+
			if ($readdir1 eq "+") { $keydir1 = "-"; }
			if ($readdir2 eq "-") { $keydir2 = "-"; }
		} else {
			print STDERR "\nError: read orientation $insertreaddir\n";
			die;
		}

	$key1 .= $keydir1;
	$key2 .= $keydir2;

	# adjust keys
	if ($scfname1 eq $scfname2) {
		if ($block_start2 < $block_start1) {
			my $keytmp = $key1;
			$key1 = $key2;
			$key2 = $keytmp;
		}
	} elsif ($scfname1 gt $scfname2) {
		if ($keydir1 eq "+") { $keydir1 = "-"; }
		else { $keydir1 = "+"; }
			
		if ($keydir2 eq "+") { $keydir2 = "-"; }
		else { $keydir2 = "+"; }
	
		$key1 = "$scfname2:$block_start2-$block_end2 $keydir2";
		$key2 = "$scfname1:$block_start1-$block_end1 $keydir1";
	}
	
	$peused++;
	
    # store
	if (defined($mappings{$key1}{$key2})) {
		$mappings{$key1}{$key2}++;
	} else {
		$mappings{$key1}{$key2} = 1;
	}
			
	#$preadname = $readname;
}
close(F);

printf STDERR "%s: %d reads used, %d unpaired\n", $insertname, $peused, $upcnt;

# print results
open(O,">$out_f");
foreach my $key1 (sort keys %mappings) {
	my $rhs = $mappings{$key1};
	my ($crd1, $dir1) = split(/\s+/, $key1);
	my $bid1 = $hs_blocks_id{$crd1};
	my $newcrd1 = $hs_blocks_map{$bid1};

	foreach my $key2 (sort keys %$rhs) {
		my $cnt = $$rhs{$key2};
		my ($crd2, $dir2) = split(/\s+/, $key2);
		my $bid2 = $hs_blocks_id{$crd2};
		my $newcrd2 = $hs_blocks_map{$bid2};
	
		print O "$newcrd1 $dir1\t$newcrd2 $dir2\t$cnt\n";
	}
}

##################################################################
sub read_search {
	my $pestart = shift;
	my $peend = shift;
	my $scfnum = shift;
	my $rhs_blocks = shift;
			
	my $rhs = $$rhs_blocks{$scfnum};
	my ($block_start, $block_end) = (-1,-1);
	foreach my $start (sort {$a<=>$b} keys %$rhs) {
		my $end = $$rhs{$start};
		if ($pestart >= $start && $peend <= $end) {
			$block_start = $start;
			$block_end = $end;
			last;
		}
	}

	return ($block_start, $block_end);
}
