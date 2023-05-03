#!/usr/bin/perl

use strict;
use warnings;

my $rspc = shift;   # reference species
my $tspc = shift;   # target species
my $f = shift;      # Genomes.Order

my ($frspc, $ftspc) = (0, 0);
my %hs_ref = ();
my $chr = "";
my $total_scf = 0;
my $inter_scf = 0;
open(F, $f); 
while(<F>) {
    chomp;
    if (length($_) == 0) { next; }
    if ($_ =~ /^# (\S+)/) {
        $chr = $1;
        next; 
    }

    if ($_ =~ /^>$rspc/) {
        $frspc = 1;
        $ftspc = 0;
    } elsif ($_ =~ /^>$tspc/) {
        $frspc = 0;
        $ftspc = 1;
    } elsif ($_ =~ /^>/) {
        last;
    } else {
        if ($frspc == 1) {
            my @ar = split(/\s+/);
            pop(@ar);
            for (my $i = 0; $i <= $#ar; $i++) {
                $hs_ref{abs($ar[$i])} = $chr;
            }
        } 
        
        if ($ftspc == 1) {
            $total_scf++;
            my @ar = split(/\s+/);
            pop(@ar);
            my %hs_tspc = ();
            for (my $i = 0; $i <= $#ar; $i++) {
                my $bid = abs($ar[$i]);
                my $rchr = $hs_ref{$bid};
                $hs_tspc{$rchr} = 1;
            }

            my @rchrs = keys %hs_tspc;
            if (scalar(@rchrs) > 1) {
                $inter_scf++;
                print "$chr\t@rchrs\n";
            }
        }
    }
}
close(F);

printf "$inter_scf\t$total_scf\t%.2f\n", $inter_scf/$total_scf;
