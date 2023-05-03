#!/usr/bin/perl

use strict;
use warnings;

my $f = shift;

my ($spc, $seq) = ("", "");
open(F, $f);
while(<F>) {
	chomp;
	if (length($_) == 0 || $_ =~ /^#/) { next; }
	if ($_ =~ /^>(\S+)/) {
		if ($spc ne "") {
			print "$spc\t" . length($seq) . "\n";
		}
		($spc, $seq) = ($1, "");
	} else {
		$seq .= $_;
	}
}
close(F);
		
if ($spc ne "") {
	print "$spc\t" . length($seq) . "\n";
}
