#!/usr/bin/perl

use strict;
use warnings;

my $block_f = shift;
my $fa_f = shift;
my $src_f = shift;

my $seed_len = 5;
my $flank_size = 100;
my $gapsize = 100;
my $pid_thr = 1;

my %sequences_ft = ();
my %sequences_bk = ();
my %hs_scfsize = ();
my $name = "";
my $seq = "";
open(F,"$fa_f");
while(<F>) {
    chomp;
    if ($_ =~ /^>/) {
        if (length($seq) > 0) {
			my $front = substr($seq, 0, $flank_size);
			my $back = substr($seq, -$flank_size);
            $sequences_ft{$name} = $front;
            $sequences_bk{$name} = $back;
			$hs_scfsize{$name} = length($seq);
        }

        my @ar = split(/\s+/);
        $name = substr($ar[0],1);
        $seq = "";
    } else {
        $seq .= $_;
    }
}
close(F);
my $front = substr($seq, 0, $flank_size);
my $back = substr($seq, -$flank_size);
$sequences_ft{$name} = $front;
$sequences_bk{$name} = $back;
$hs_scfsize{$name} = length($seq);

my %hs_blocks = ();
my %hs_blocks_details = ();
open(F,"$block_f");
while(<F>) {
	chomp;
	my ($scf, $start, $end, $dir, $bid) = split(/\s+/);
	my $key = "$scf $start $end";
	$hs_blocks{$key} = 1;	
	$hs_blocks_details{$scf}{$start} = $end;	
}
close(F);

my $chr = "";
my $newstart = 0;
my ($pscf, $pstart, $pend, $pdir) = ("",0,0,"");
my $poverlap = 0;
open(F,"$src_f");
while(<F>) {
	chomp;
	if (length($_) == 0) { 
		# obtain scaffold info.
		my ($scfstart, $scfend, $blk_order, $blk_last) = get_scfpos($pscf, $pstart, $pend); 
		
		my $newend = $newstart + ($scfend-$scfstart) - $poverlap;
		if ($pdir eq "+") { $scfstart += $poverlap; }
		else { $scfend -= $poverlap; }	
		print "$chr\t$newstart\t$newend\t$pscf\t$scfstart\t$scfend\t$pdir\t$pstart\t$pend\n";
		$newstart = $newend;
	
		my $key = "$pscf $pstart $pend";
		delete $hs_blocks{$key};
	
		next; 
	}

	if ($_ =~ /^>(\S+)/) {
		$chr = $1;
		$newstart = 0;
		($pscf, $pstart, $pend, $pdir) = ("",0,0,"");
		$poverlap = 0;
	} else {
		my @ar = split(/\s+/);
		my ($block1, $dir1, $block2, $dir2) = ($ar[0], $ar[1], $ar[2], $ar[3]);
	
		$block1 =~ /(\S+):(\S+)\-(\S+)/;
		my ($scf, $start, $end) = ($1,$2,$3);

		# obtain scaffold info.
		my ($scfstart, $scfend, $blk_order, $blk_last) = get_scfpos($scf, $start, $end); 

		my $newend = $newstart + ($scfend-$scfstart) - $poverlap;	
		if ($dir1 eq "+") { $scfstart += $poverlap; }
		else { $scfend -= $poverlap; }
		print "$chr\t$newstart\t$newend\t$scf\t$scfstart\t$scfend\t$dir1\t$start\t$end\n";
		$newstart = $newend;
		
		my $key = "$scf $start $end";
		delete $hs_blocks{$key};
	
		$block2 =~ /(\S+):(\S+)\-(\S+)/;
		($pscf, $pstart, $pend) = ($1,$2,$3);
		$pdir = $dir2;

		# check overlap
		if ($scf eq $pscf) {
			my $order1 = get_blockorder($scf, $start, $end);
            my $order2 = get_blockorder($scf, $pstart, $pend);
            if ($dir1 ne $dir2 || abs($order1 - $order2) != 1 || ($dir1 eq "+" && !($pstart>$start && $pend>$end)) || ($dir1 eq "-" && !($pstart<$start && $pend<$end))) {
				$newend = $newstart + $gapsize;
				print "$chr\t$newstart\t$newend\tGAPS\n";
				$newstart = $newend;
				$poverlap = 0;
			} 
		} elsif (($blk_order == 0 && $dir1 eq "+") ||
				 ($blk_order == $blk_last && $dir1 eq "-") ||
				 ($blk_order > 0 && $blk_order < $blk_last)) {
			$newend = $newstart + $gapsize;
			print "$chr\t$newstart\t$newend\tGAPS\n";
			$newstart = $newend;
			$poverlap = 0;
		} else {
			$newend = $newstart + $gapsize;
			print "$chr\t$newstart\t$newend\tGAPS\n";
			$newstart = $newend;
			$poverlap = 0;
		}	
	}	
}
close(F);

foreach my $key (keys %hs_blocks) {
	$chr++;	
	$newstart = 0;
	my ($scf, $start, $end) = split(/\s+/, $key);
		
	# obtain scaffold info.
	my ($scfstart, $scfend, $blk_order, $blk_last) = get_scfpos($scf, $start, $end); 

	my $newend = $newstart + ($scfend-$scfstart);	
	print "$chr\t$newstart\t$newend\t$scf\t$scfstart\t$scfend\t+\t$start\t$end\n";
	$newstart = $newend;
}

sub get_scfpos {
	my $s_scf = shift;
	my $s_start = shift;
	my $s_end = shift;	
	
	my $s_blk_order = -1;
	my ($s_scfstart, $s_scfend) = (0,0);
	my $s_scfsize = $hs_scfsize{"$s_scf"};
	my $s_rhs_blkstarts = $hs_blocks_details{$s_scf};	
	my @s_blk_starts = sort {$a<=>$b} keys %$s_rhs_blkstarts;
	if (scalar(@s_blk_starts) == 1) {
		$s_scfstart = 0;
		$s_scfend = $s_scfsize;
	} else {
		$s_blk_order = 0;
		foreach my $s_blk_start (@s_blk_starts) {
			if ($s_blk_start == $s_start) { last; }
			$s_blk_order++; 
		}
		if ($s_blk_order == 0) {
			$s_scfstart = 0;
			my $s_next_blkstart = $s_blk_starts[$s_blk_order+1];
			$s_scfend = int(($s_end + $s_next_blkstart)/2);	
		} elsif ($s_blk_order == $#s_blk_starts) {
			$s_scfend = $s_scfsize;
			my $s_prev_blkstart = $s_blk_starts[$s_blk_order-1];
			my $s_prev_blkend = $$s_rhs_blkstarts{$s_prev_blkstart};
			$s_scfstart = int(($s_start + $s_prev_blkend)/2);	
		} else {
			my $s_prev_blkstart = $s_blk_starts[$s_blk_order-1];
			my $s_prev_blkend = $$s_rhs_blkstarts{$s_prev_blkstart};
			$s_scfstart = int(($s_start + $s_prev_blkend)/2);	
			
			my $s_next_blkstart = $s_blk_starts[$s_blk_order+1];
			$s_scfend = int(($s_end + $s_next_blkstart)/2);	
		}
	}

	return ($s_scfstart, $s_scfend, $s_blk_order, $#s_blk_starts); 
}

sub get_blockorder {
    my $s_scf = shift;
    my $s_start = shift;
    my $s_end = shift;
    my $s_rhs_blkstarts = $hs_blocks_details{$s_scf};
    my @s_blk_starts = sort {$a<=>$b} keys %$s_rhs_blkstarts;
    die if (scalar(@s_blk_starts) == 1);

    my $blk_order = 1;
    foreach my $s_blk_start (@s_blk_starts) {
        if ($s_blk_start == $s_start) { last; }
        $blk_order++;
    }

    return $blk_order;
}

