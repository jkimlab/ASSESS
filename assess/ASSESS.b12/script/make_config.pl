#!/usr/bin/perl

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../lib/perl";
use Bio::TreeIO;

my $refspc = shift;
my $tarspc = shift;
my $res = shift;
my $cn_dir = shift;
my $tree_f = shift;

print ">netdir\n";
print "$cn_dir\n\n";

print ">chaindir\n";
print "$cn_dir\n\n";

print ">resolution\n";
print "$res\n\n";

print ">species\n";
print "$refspc 0\n";
print "$tarspc 1\n";

my $in = new Bio::TreeIO(-format => "newick", -file => "$tree_f");
my $tree = $in->next_tree;
my $tarnode = $tree->find_node(-id => "$tarspc@");
my $refnode = $tree->find_node(-id => "$refspc");
my $lcanode = $tree->find_node(-id => "lca#");

my %hs_used = ();
foreach my $node ($lcanode->get_all_Descendents) {
	if (!$node->is_Leaf) { next; }
	if ($node->id eq $tarnode->id) { next; }
	if ($node->id eq $refnode->id) { next; }

	print $node->id, " 1\n";
	$hs_used{"$node->id"} = 1;
}

foreach my $node ($tree->get_nodes) {
	if (!$node->is_Leaf) { next; }
	if ($node->id eq $tarnode->id) { next; }
	if ($node->id eq $refnode->id) { next; }
	if (defined($hs_used{"$node->id"})) { next; }	

	print $node->id, " 2\n";
}
