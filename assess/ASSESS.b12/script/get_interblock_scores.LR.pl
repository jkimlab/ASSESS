#!/usr/bin/perl

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../lib/perl";
use File::Basename;

my $blist_f     = shift;
my $blistext_f  = shift;
my $data_f      = shift;
my $out_f       = shift;
my $min_mapq    = shift;

my %blocks      = ();
my %sfPos       = ();
my %sfSeq       = ();
my %sfInfo      = ();
my %sizeInfo    = ();
my %mInfo_sf    = ();
my %mInfo_scf   = ();
my %counts      = ();

my $scf_ol_co  = 21;
my $ol_co      = 1;
my $samtools = "$Bin/../util/samtools-1.5/samtools";

open (RF, $blist_f);
while (<RF>){
	chomp;
	
    my @col = split(/\s+/);
    my ( $scf, $beg, $end, $dir, $id ) = ( $col[0], $col[1], $col[2], $col[3], $col[4] );
	$blocks{$id} = "$scf:$beg-$end	$dir";
}
close RF;

# 1.2 Get extended synteny position information
open( RF, $blistext_f );
while (<RF>) {
    chomp;

    my @col = split(/\s+/);
    my ( $scf, $beg, $end, $dir, $id ) = ( $col[0], $col[1], $col[2], $col[3], $col[4] );
    $sfPos{$scf}{$id}{"beg"} = $beg;
    $sfPos{$scf}{$id}{"end"} = $end;

    $sfInfo{$id}{"scf"} = $scf;
    $sfInfo{$id}{"beg"} = $beg;
    $sfInfo{$id}{"dir"} = $dir;

    $sizeInfo{"sf"}{$id} = $end - $beg;
}
close RF;

# 2 Get mapping position in scaffold/synteny
open( RF,  "$samtools view -q $min_mapq -F 0x04 $data_f |" );
while (<RF>) {
    chomp;
    my (
        $lr,       $mBeg_lr, $mEnd_lr, $scf, $mBeg_scf,
        $mEnd_scf, $flag,    $cigar,   $strand
    );

    my @col = split(/\s+/);
    ( $lr, $flag, $scf, $mBeg_scf, $cigar ) =
      ( $col[0], $col[1], $col[2], --$col[3], $col[5] );
    $strand = ( $flag & 0x0010 ) ? "-" : "+";

    my @mPos_lr  = ();
    my $scfid    = 0;
    my $mLen_scf = 0;
    my $lrLen    = 0;

    if ( my @mapPat = $cigar =~ /(\d+)([M|I|D|N|S|H|P|=|X])/g ) {
        for ( my $pi = 0 ; $pi <= $#mapPat ; $pi += 2 ) {
            my $len  = $mapPat[$pi];
            my $type = $mapPat[ $pi + 1 ];

            if ( ($type eq "M" || $type eq "=" || $type eq "X") && $#mPos_lr != 0 ) { push( @mPos_lr, $lrLen ); }
            if ( ( $type eq "S" || $type eq "H" ) && $#mPos_lr == 0 ) {
                push( @mPos_lr, $lrLen );
            }

            if (   $type eq "M"
                || $type eq "D"
                || $type eq "N"
                || $type eq "="
                || $type eq "X" )
            {
                $mLen_scf += $len;
            }
            if (   $type eq "M"
                || $type eq "I"
                || $type eq "S"
                || $type eq "H"
                || $type eq "="
                || $type eq "X" )
            {
                $lrLen += $len;
            }
        }
    }

    $mEnd_scf = $mBeg_scf + $mLen_scf;
    if ( $#mPos_lr == 0 ) { push( @mPos_lr, $lrLen ); }

    if ( $strand eq "-" ) {
        $mBeg_lr = $lrLen - $mPos_lr[1];
        $mEnd_lr = $lrLen - $mPos_lr[0];
    }
    else {
        $mBeg_lr = $mPos_lr[0];
        $mEnd_lr = $mPos_lr[1];
    }

    if ( exists $mInfo_scf{$lr} ) {
        $scfid = scalar @{ $mInfo_scf{$lr} };
        push(
            @{ $mInfo_scf{$lr} },
            {
                "scf" => $scf,
                "beg" => $mBeg_lr,
                "end" => $mEnd_lr,
                "id"  => $scfid
            }
        );
    }
    else {
        $mInfo_scf{$lr} = [
            {
                "scf" => $scf,
                "beg" => $mBeg_lr,
                "end" => $mEnd_lr,
                "id"  => $scfid
            }
        ];
    }
    $sizeInfo{"lr"}{$lr} = $lrLen;

    foreach my $sf ( keys %{ $sfPos{$scf} } ) {
        my $curBeg_lr  = $mBeg_lr;
        my $curEnd_lr  = $mEnd_lr;
        my $curBeg_scf = $mBeg_scf;
        my $curEnd_scf = $mEnd_scf;

        my $mBeg_sf_scf = $sfPos{$scf}{$sf}{"beg"};
        my $mEnd_sf_scf = $sfPos{$scf}{$sf}{"end"};

        if (   ( $mBeg_sf_scf <= $curBeg_scf && $curBeg_scf < $mEnd_sf_scf )
            || ( $mBeg_sf_scf < $curEnd_scf  && $curEnd_scf <= $mEnd_sf_scf )
            || ( $curBeg_scf <= $mBeg_sf_scf && $mBeg_sf_scf < $curEnd_scf )
            || ( $curBeg_scf < $mEnd_sf_scf  && $mEnd_sf_scf <= $curEnd_scf ) )
        {
            if ( $mBeg_scf <= $mBeg_sf_scf && $mBeg_sf_scf < $mEnd_scf ) {
                my $curPos = $mBeg_scf;
                $curBeg_scf = $mBeg_sf_scf;

                my $mapLen_sf = 0;
                my $map_beg   = 0;
                if ( my @mapPat = $cigar =~ /(\d+)([M|I|D|N|S|H|P|=|X])/g ) {
                    for ( my $pi = 0 ; $pi <= $#mapPat ; $pi += 2 ) {
                        my $len  = $mapPat[$pi];
                        my $type = $mapPat[ $pi + 1 ];

                        if ( $type eq "M" ) { $map_beg = 1; }
                        if (
                            $map_beg
                            && (   $type eq "M"
                                || $type eq "D"
                                || $type eq "N"
                                || $type eq "="
                                || $type eq "X" )
                          )
                        {
                            $curPos += $len;
                        }
                        if (   $type eq "M"
                            || $type eq "I"
                            || $type eq "S"
                            || $type eq "H"
                            || $type eq "="
                            || $type eq "X" )
                        {
                            $mapLen_sf += $len;
                        }

                        if ( $curPos >= $curBeg_scf ) {
                            if ( $type ne "N" ) {
                                $mapLen_sf -= ( $curPos - $curBeg_scf );
                            }
                            else {
                                $curBeg_scf += ( $curPos - $curBeg_scf );
                            }

                            if ( $strand eq "+" ) { $curBeg_lr = $mapLen_sf; }
                            if ( $strand eq "-" ) {
                                $curEnd_lr = ( $lrLen - $mapLen_sf );
                            }

                            last;
                        }
                    }
                }
            }

            if ( $mBeg_scf < $mEnd_sf_scf && $mEnd_sf_scf < $mEnd_scf ) {
                my $curPos = $mBeg_scf;
                $curEnd_scf = $mEnd_sf_scf;

                my $mapLen_sf = 0;
                my $map_beg   = 0;
                if ( my @mapPat = $cigar =~ /(\d+)([M|I|D|N|S|H|P|=|X])/g ) {
                    for ( my $pi = 0 ; $pi <= $#mapPat ; $pi += 2 ) {
                        my $len  = $mapPat[$pi];
                        my $type = $mapPat[ $pi + 1 ];

                        if ( $type eq "M" ) { $map_beg = 1; }
                        if (
                            $map_beg
                            && (   $type eq "M"
                                || $type eq "D"
                                || $type eq "N"
                                || $type eq "="
                                || $type eq "X" )
                          )
                        {
                            $curPos += $len;
                        }
                        if (   $type eq "M"
                            || $type eq "I"
                            || $type eq "S"
                            || $type eq "H"
                            || $type eq "="
                            || $type eq "X" )
                        {
                            $mapLen_sf += $len;
                        }

                        if ( $curPos >= $curEnd_scf ) {
                            if ( $type ne "N" ) {
                                $mapLen_sf -= ( $curPos - $curEnd_scf );
                            }
                            else {
                                $curEnd_scf = ( $curPos - $len );
                            }

                            if ( $strand eq "+" ) { $curEnd_lr = $mapLen_sf; }
                            if ( $strand eq "-" ) {
                                $curBeg_lr = ( $lrLen - $mapLen_sf );
                            }

                            last;
                        }
                    }
                }
            }
            if ( $curBeg_lr == $curEnd_lr ) { next; }    # Alignment gap

            my $mBeg_sf = $curBeg_scf - $mBeg_sf_scf;
            my $mEnd_sf = $curEnd_scf - $mBeg_sf_scf;

            if ( exists $mInfo_sf{$lr} ) {
                my $curid = scalar @{ $mInfo_sf{$lr} };
                push(
                    @{ $mInfo_sf{$lr} },
                    {
                        "sf"     => $sf,
                        "sf_beg" => $mBeg_sf,
                        "sf_end" => $mEnd_sf,
                        "lr_beg" => $curBeg_lr,
                        "lr_end" => $curEnd_lr,
                        "dir"    => $strand,
                        "scfid"  => $scf
                    }
                );
            }
            else {
                $mInfo_sf{$lr} = [
                    {
                        "sf"     => $sf,
                        "sf_beg" => $mBeg_sf,
                        "sf_end" => $mEnd_sf,
                        "lr_beg" => $curBeg_lr,
                        "lr_end" => $curEnd_lr,
                        "dir"    => $strand,
                        "id"     => 0,
                        "scfid"  => $scf
                    }
                ];
            }
        }
    }
}

# 3.Detect synteny pairs from each long reads
foreach my $lr ( keys %mInfo_sf ) {
    my %fltout_scf_link  = ();
    my %group_aln        = ();

    my %cluster_mDir     = ();
    my %cluster_position = ();
    my $gid              = 0;

    my @mInfo_sf     = @{ $mInfo_sf{$lr} };
    my @mInfo_lr_scf = @{ $mInfo_scf{$lr} };

    my %rel_sfs;
    for ( my $i = 0 ; $i <= $#mInfo_sf ; ++$i ) {
        my $cursf = $mInfo_sf[$i]->{"sf"};
        $rel_sfs{$cursf} = 1;
    }
    if ( ( scalar keys %rel_sfs ) <= 1 ) { next; }

    # 3.1 Check overlaps among mapping positions in scaffolds
    @mInfo_lr_scf = sort { $$a{"beg"} <=> $$b{"beg"} || $$a{"end"} <=> $$b{"end"} } @mInfo_lr_scf;
    for ( my $i = 0 ; $i < $#mInfo_lr_scf ; ++$i ) {
        my $scf_i = $mInfo_lr_scf[$i]->{"scf"};
        my $beg_i = $mInfo_lr_scf[$i]->{"beg"};
        my $end_i = $mInfo_lr_scf[$i]->{"end"};

        for ( my $j = $i + 1 ; $j <= $#mInfo_lr_scf ; ++$j ) {
            my $scf_j = $mInfo_lr_scf[$j]->{"scf"};
            my $beg_j = $mInfo_lr_scf[$j]->{"beg"};
            my $end_j = $mInfo_lr_scf[$j]->{"end"};

            if ( $end_i <= $beg_j ) { last; }
            my $o_len = ( $end_i - $beg_j < $end_j - $beg_j ) ? $end_i - $beg_j : $end_j - $beg_j;

            if ( $scf_i ne $scf_j && $o_len >= $scf_ol_co ) {
                my $scfid_i = $mInfo_lr_scf[$i]->{"id"};
                my $scfid_j = $mInfo_lr_scf[$j]->{"id"};
                $fltout_scf_link{$scfid_i}{$scfid_j} = 1;

                print STDERR "[! scaffold] $lr\t$scf_i\t$scf_j\t$o_len\n";
                print STDERR "\t$scf_i: $beg_i\t$end_i\n";
                print STDERR "\t$scf_j: $beg_j\t$end_j\n";
            }
        }
    }

    @mInfo_sf = sort { $$a{"lr_beg"} <=> $$b{"lr_beg"} || $$a{"lr_end"} <=> $$b{"lr_end"} } @mInfo_sf;

    # 3.2 Cluster the alignments
    for ( my $i = 0 ; $i <= $#mInfo_sf ; ++$i ) {
        my $newSf     = $mInfo_sf[$i]->{"sf"};
        my $newDir    = $mInfo_sf[$i]->{"dir"};
        my $newBeg_sf = $mInfo_sf[$i]->{"sf_beg"};
        my $newEnd_sf = $mInfo_sf[$i]->{"sf_end"};

        my @groups = keys %group_aln;
        if ( $#groups >= 0 ) {
            my $save = 0;

            foreach my $g (@groups) {
                my $saveSf     = $group_aln{$g}->[-1]->{"sf"};
                my $saveBeg_sf = $group_aln{$g}->[-1]->{"sf_beg"};
                my $saveEnd_sf = $group_aln{$g}->[-1]->{"sf_end"};
                my $saveDir    = $group_aln{$g}->[-1]->{"dir"};

                if ( $newSf eq $saveSf && $newDir eq $saveDir && ( ( $newDir eq "+" && $newBeg_sf >= $saveBeg_sf ) || ( $newDir eq "-" && $newEnd_sf <= $saveEnd_sf ) ) )
                {
                    push( @{ $group_aln{$g} }, $mInfo_sf[$i] );
                    $save = 1;
                }
            }

            if ($save) { next; }
        }
        $group_aln{"g$gid"} = [ $mInfo_sf[$i] ];
        ++$gid;
    }

    my %unwanted    = ();
    my %check_dup_g = ();
    foreach my $g ( keys %group_aln ) {
        my @mInfo_sf_g = @{ $group_aln{$g} };
        my $sf_g       = $mInfo_sf_g[0]->{"sf"};

        if ( exists $check_dup_g{$sf_g} ) {
            my $dup_g = $check_dup_g{$sf_g};
            $unwanted{$dup_g} = 1;
            $unwanted{$g}     = 1;
        }
        else {
            $check_dup_g{$sf_g} = $g;
        }
    }

    if ( scalar( keys %unwanted ) >= 1 ) {
        my %mul_sfs = ();
        foreach my $g ( keys %unwanted ) {
            $mul_sfs{ $group_aln{$g}[0]->{'sf'} } = 1;
        }
        my $mul_sf_str = join( ", ", keys %mul_sfs );
        print STDERR "[! Multiple alignment groups in a SF] $lr\t$mul_sf_str\n";

        delete @group_aln{ keys %unwanted };
    }

    my @g_cnt = keys %group_aln;
    if ( $#g_cnt <= 0 ) { next; }

    # 3.3 Calculate mapping/overlap positions for groups in the read
    foreach my $g ( keys %group_aln ) {
        my @sort_lrPos = ();
        my @sort_sfPos = ();
        my $sf_g;
        my $dir_g;
        my $mBeg_sf_min_g = 0;
        my $mBeg_lr_min_g = 0;
        my $mEnd_sf_max_g = 0;
        my $mEnd_lr_max_g = 0;
        my $alnLen_sf_g   = 0;
        my $alnLen_lr_g   = 0;
        my $oBeg_lr_g;
        my $oEnd_lr_g;

        my @mInfo_sf_g = @{ $group_aln{$g} };

        $sf_g          = $mInfo_sf_g[0]->{"sf"};
        $dir_g         = $mInfo_sf_g[0]->{"dir"};
        $mBeg_sf_min_g = $sizeInfo{"sf"}{$sf_g};
        $mBeg_lr_min_g = $sizeInfo{"lr"}{$lr};

        for ( my $i = 0 ; $i <= $#mInfo_sf_g ; ++$i ) {
            my $lr_beg_i = $mInfo_sf_g[$i]->{"lr_beg"};
            my $lr_end_i = $mInfo_sf_g[$i]->{"lr_end"};

            my $sf_beg_i = $mInfo_sf_g[$i]->{"sf_beg"};
            my $sf_end_i = $mInfo_sf_g[$i]->{"sf_end"};
			
			print STDERR "\t*$lr\t$g\t$lr_beg_i\t$lr_end_i\t$sf_g\t$sf_beg_i\t$sf_end_i\n";

            if ( $mBeg_sf_min_g > $sf_beg_i ) {
                $mBeg_sf_min_g = $sf_beg_i;
            }
            if ( $mBeg_lr_min_g > $lr_beg_i ) {
                $mBeg_lr_min_g = $lr_beg_i;
            }
            if ( $mEnd_sf_max_g < $sf_end_i ) {
                $mEnd_sf_max_g = $sf_end_i;
            }
            if ( $mEnd_lr_max_g < $lr_end_i ) {
                $mEnd_lr_max_g = $lr_end_i;
            }
        }

        my $fOvlLen_lr_g = $mBeg_lr_min_g;
        my $bOvlLen_lr_g = $sizeInfo{"lr"}{$lr} - $mEnd_lr_max_g;

        my $fOvlLen_sf_g = 0;
        my $bOvlLen_sf_g = 0;
        if ( $dir_g eq "+" ) {
            $fOvlLen_sf_g = $mBeg_sf_min_g;
            $bOvlLen_sf_g = $sizeInfo{"sf"}{$sf_g} - $mEnd_sf_max_g;
        }
        else {
            $fOvlLen_sf_g = $sizeInfo{"sf"}{$sf_g} - $mEnd_sf_max_g;
            $bOvlLen_sf_g = $mBeg_sf_min_g;
        }

        my $fOvlLen_g =
          ( $fOvlLen_sf_g < $fOvlLen_lr_g ) ? $fOvlLen_sf_g : $fOvlLen_lr_g;
        my $bOvlLen_g =
          ( $bOvlLen_sf_g < $bOvlLen_lr_g ) ? $bOvlLen_sf_g : $bOvlLen_lr_g;

        $oBeg_lr_g = $mBeg_lr_min_g - $fOvlLen_g;
        $oEnd_lr_g = $mEnd_lr_max_g + $bOvlLen_g;

        $cluster_position{$g}{"sf"} = $sf_g;
        $cluster_position{$g}{"ob"} = $oBeg_lr_g;
        $cluster_position{$g}{"oe"} = $oEnd_lr_g;
        $cluster_position{$g}{"mb"} = $mBeg_lr_min_g;
        $cluster_position{$g}{"me"} = $mEnd_lr_max_g;
        $cluster_mDir{$g}           = $dir_g;
    }

    # 3.4 Detect synteny pairs
    my @g_list = sort {
        $cluster_position{$a}{"mb"} <=> $cluster_position{$b}{"mb"}
          || $cluster_position{$a}{"me"} <=> $cluster_position{$b}{"me"}
    } keys %cluster_position;

    my @overlap = ( $g_list[0] );
    my @cluster = ( \@overlap );

    for ( my $i = 1 ; $i <= $#g_list ; ++$i ) {
        my $scf_fg = $sfInfo{ $cluster_position{ $g_list[$i] }{"sf"} }{"scf"};
        my $mBeg_lr_fg = $cluster_position{ $g_list[$i] }{"mb"};
        my $mEnd_lr_fg = $cluster_position{ $g_list[$i] }{"me"};

        my $save = 0;
        foreach my $pre_g (@overlap) {
            my $scf_sg     = $sfInfo{ $cluster_position{$pre_g}{"sf"} }{"scf"};
            my $mBeg_lr_sg = $cluster_position{$pre_g}{"mb"};
            my $mEnd_lr_sg = $cluster_position{$pre_g}{"me"};

            if (
                $scf_fg eq $scf_sg
                && (
                    (
                           $mBeg_lr_fg <= $mBeg_lr_sg
                        && $mEnd_lr_sg <= $mEnd_lr_fg
                    )
                    || (   $mBeg_lr_sg <= $mBeg_lr_fg
                        && $mEnd_lr_fg <= $mEnd_lr_sg )
                )
              )
            {
                push( @{ $cluster[-1] }, $g_list[$i] );
                $save = 1;
                last;
            }
        }

        if ( !$save ) {
            my @new_ol = ( $g_list[$i] );
            push( @cluster, \@new_ol );
        }
    }

    for ( my $i = 0 ; $i <= $#cluster ; ++$i ) {
        @{ $cluster[$i] } = sort {
            ( $sfInfo{ $cluster_position{$a}{"sf"} }{"beg"} )
              <=> ( $sfInfo{ $cluster_position{$b}{"sf"} }{"beg"} )
        } @{ $cluster[$i] };
    }

    my %head = ();
    for ( my $i = 0 ; $i <= $#cluster ; ++$i ) {
		if( scalar keys %head ){ last; }
    	for ( my $j = 0 ; $j <= $#{ $cluster[$i] } ; ++$j ) {
			&Find_sf_link( $lr, "h", $i, $j, \@cluster, \%cluster_position, \%cluster_mDir, \%sfInfo, \%fltout_scf_link, \%head );
			if( scalar keys %head ) { last; }
		}
	}

	@cluster = reverse @cluster;
    for ( my $i = 0 ; $i <= $#cluster ; ++$i ) {
        @{ $cluster[$i] } = reverse @{ $cluster[$i] };
    }

    my %tail = ();
    for ( my $i = 0 ; $i <= $#cluster ; ++$i ) {
		if( scalar keys %tail ){ last; }
    	for ( my $j = 0 ; $j <= $#{ $cluster[$i] } ; ++$j ) {
    		&Find_sf_link( $lr, "t", $i, $j, \@cluster, \%cluster_position, \%cluster_mDir, \%sfInfo, \%fltout_scf_link, \%tail );
			if( scalar keys %tail ){ last; }
		}
	}

    foreach my $sf1 ( keys %head ) {
        foreach my $sf2 ( keys %{ $head{$sf1} } ) {
			my $scf1 = $sfInfo{abs ($sf1)}{"scf"};
			my $scf2 = $sfInfo{abs ($sf2)}{"scf"};

            if ( exists $tail{$sf1}{$sf2} && $scf1 ne $scf2 ) {
                if ( exists $counts{$sf1}{$sf2} ) {
                    ++$counts{$sf1}{$sf2};
                }
                else {
                    $counts{$sf1}{$sf2} = 1;
                }
            }
        }
    }
}

open (WF, ">$out_f");
# 4. Print results
foreach my $fsf ( sort { $a cmp $b } keys %counts ) {
	my $fcrd; my $fdir; my $fbid;

	if($fsf =~ /([-]*)(\w+)/){
		my $dir;
		( $fcrd, $dir ) = split( /\s+/, $blocks{$2} );
		$fdir = ( $dir eq $1 )? "+" : "-";
	}

    my $sfs_linked_fsf = $counts{$fsf};
    foreach my $ssf ( sort { $a cmp $b } keys %$sfs_linked_fsf ) {
		my $scrd; my $sdir; my $sbid;

		if($ssf =~ /([-]*)(\w+)/){
			my $dir;
			($scrd, $dir) = split( /\s+/, $blocks{$2});
			$sdir = ( $dir eq $1 )? "+" : "-";
		}
        my $cnt   = $counts{$fsf}{$ssf};

        print WF "$fcrd $fdir\t$scrd $sdir\t$cnt\t$fsf\t$ssf\n";
    }
}
close WF;

#-- Funtions 
sub Find_sf_link {
    my $_lr    = shift;
    my $_hot   = shift;
    my $_anc_i = shift;
    my $_anc_j = shift;

    my @_cluster          = @{ shift() };
    my %_cluster_pos      = %{ shift() };
    my %_cluster_dir      = %{ shift() };
    my %_sf_info          = %{ shift() };
    my %_unwanted_scflink = %{ shift() };
    my $_results          = shift();

    for ( my $i = $_anc_i ; $i <= $#_cluster ; ++$i ) {
        for ( my $j = $_anc_j ; $j <= $#{ $_cluster[$i] } ; ++$j ) {
            if ( $i == $_anc_i && $j == $_anc_j ) { next; }
            my $sf_fg;
            my $sf_sg;
            my $dir_fg;
            my $dir_sg;
            my $cBeg_lr_fg;
            my $cEnd_lr_fg;
            my $cBeg_lr_sg;
            my $cEnd_lr_sg;

            my $g_overlap          = 0;
            my $score_pair         = 0;
            my $mBeg_lr_fg         = 0;
            my $mBeg_lr_sg         = 0;
            my $mEnd_lr_fg         = 0;
            my $mEnd_lr_sg         = 0;
            my $map_overlap_region = 0;

            $sf_fg = $_cluster_pos{ $_cluster[$_anc_i][$_anc_j] }{"sf"};

            $cBeg_lr_fg = $_cluster_pos{ $_cluster[$_anc_i][$_anc_j] }{"ob"};
            $cEnd_lr_fg = $_cluster_pos{ $_cluster[$_anc_i][$_anc_j] }{"oe"};

            $mBeg_lr_fg = $_cluster_pos{ $_cluster[$_anc_i][$_anc_j] }{"mb"};
            $mEnd_lr_fg = $_cluster_pos{ $_cluster[$_anc_i][$_anc_j] }{"me"};

            if ( $_anc_i != $i ) {
                $sf_sg = $_cluster_pos{ $_cluster[$i][$j] }{"sf"};

                $cBeg_lr_sg = $_cluster_pos{ $_cluster[$i][$j] }{"ob"};
                $cEnd_lr_sg = $_cluster_pos{ $_cluster[$i][$j] }{"oe"};

                $mBeg_lr_sg = $_cluster_pos{ $_cluster[$i][$j] }{"mb"};
                $mEnd_lr_sg = $_cluster_pos{ $_cluster[$i][$j] }{"me"};

                if ( $_hot eq "h" ) {
                    $dir_fg = ( $_cluster_dir{ $_cluster[$_anc_i][$_anc_j] } eq $_sf_info{$sf_fg}{"dir"} ) ? "" : "-";
                    $dir_sg = ( $_cluster_dir{ $_cluster[$i][$j] } eq $_sf_info{$sf_sg}{"dir"} ) ? "" : "-";
                }
                else {
                    $dir_fg = ( $_cluster_dir{ $_cluster[$_anc_i][$_anc_j] } eq $_sf_info{$sf_fg}{"dir"} ) ? "-" : "";
                    $dir_sg = ( $_cluster_dir{ $_cluster[$i][$j] } eq $_sf_info{$sf_sg}{"dir"} ) ? "-" : "";
                }
            }
            else {
                $sf_sg = $_cluster_pos{ $_cluster[$i][$j] }{"sf"};

                $cBeg_lr_sg = $_cluster_pos{ $_cluster[$i][$j] }{"ob"};
                $cEnd_lr_sg = $_cluster_pos{ $_cluster[$i][$j] }{"oe"};

                $mBeg_lr_sg = $_cluster_pos{ $_cluster[$i][$j] }{"mb"};
                $mEnd_lr_sg = $_cluster_pos{ $_cluster[$i][$j] }{"me"};

                if ( $_hot eq "h" ) {
                    $dir_fg = ( $_sf_info{$sf_fg}{"dir"} eq "+" ) ? "" : "-";
                    $dir_sg = ( $_sf_info{$sf_sg}{"dir"} eq "+" ) ? "" : "-";
                }
                else {
                    $dir_fg = ( $_sf_info{$sf_fg}{"dir"} eq "+" ) ? "-" : "";
                    $dir_sg = ( $_sf_info{$sf_sg}{"dir"} eq "+" ) ? "-" : "";
                }

                $map_overlap_region = 1;
            }

            my $scf_fg = $_sf_info{$sf_fg}{"scf"};
            my $scf_sg = $_sf_info{$sf_sg}{"scf"};

            if (   exists $_unwanted_scflink{$scf_fg}{$scf_sg}
                || exists $_unwanted_scflink{$scf_sg}{$scf_fg} )
            {
                print STDERR "[! Scaffold link] $_lr\t$sf_fg\t$sf_sg\n";
                print STDERR "$sf_fg: $scf_fg\n";
                print STDERR "$sf_sg: $scf_sg\n";
                next;
            }

            if ( $cBeg_lr_fg < $cBeg_lr_sg ) {
                if ( $cBeg_lr_sg < $cEnd_lr_fg ) {
                    if ( $cEnd_lr_fg < $cEnd_lr_sg ) {
                        $g_overlap = $cEnd_lr_fg - $cBeg_lr_sg;
                    }
                    else {
                        $g_overlap = $cEnd_lr_sg - $cBeg_lr_sg;
                    }
                }
            }
            else {
                if ( $cBeg_lr_fg < $cEnd_lr_sg ) {
                    if ( $cEnd_lr_sg < $cEnd_lr_fg ) {
                        $g_overlap = $cEnd_lr_sg - $cBeg_lr_fg;
                    }
                    else {
                        $g_overlap = $cEnd_lr_fg - $cBeg_lr_fg;
                    }
                }
            }

            if ( ! $map_overlap_region && ( $g_overlap / ( $cEnd_lr_fg - $cBeg_lr_fg ) >= $ol_co || $g_overlap / ( $cEnd_lr_sg - $cBeg_lr_sg ) >= $ol_co ))
            {
                print STDERR "[! overlap region] $_lr\t$sf_fg\t$sf_sg\n";
                print STDERR "\t$sf_fg: $scf_fg $cBeg_lr_fg-$mBeg_lr_fg-$mEnd_lr_fg-$cEnd_lr_fg\n";
                print STDERR "\t$sf_sg: $scf_sg $cBeg_lr_sg-$mBeg_lr_sg-$mEnd_lr_sg-$cEnd_lr_sg\n";
                next;
            }

            my $fdir = $dir_fg;
            my $fsf  = $sf_fg;
            my $sdir = $dir_sg;
            my $ssf  = $sf_sg;

            if ( ( $sf_fg cmp $sf_sg ) > 0 ) {
                $fdir = ( $dir_sg eq "" ) ? "-" : "";
                $fsf  = $sf_sg;
                $sdir = ( $dir_fg eq "" ) ? "-" : "";
                $ssf  = $sf_fg;
            }

            print STDERR "[" . uc($_hot) . ". Find link] $_lr\t$fdir$fsf\t$sdir$ssf\n";
            print STDERR "\t$sf_fg: $scf_fg $cBeg_lr_fg-$mBeg_lr_fg-$mEnd_lr_fg-$cEnd_lr_fg\n";
            print STDERR "\t$sf_sg: $scf_sg $cBeg_lr_sg-$mBeg_lr_sg-$mEnd_lr_sg-$cEnd_lr_sg\n";

            if ( exists $_results->{"$fdir$fsf"}{"$sdir$ssf"} ) {
                ++$_results->{"$fdir$fsf"}{"$sdir$ssf"};
            }
            else {
                $_results->{"$fdir$fsf"}{"$sdir$ssf"} = 1;
            }

            &Find_sf_link( $_lr, $_hot, $i, $j, \@_cluster, \%_cluster_pos, \%_cluster_dir, \%_sf_info, \%_unwanted_scflink, $_results );
            return 1;
        }
    }
	return 0;
}
