#!/usr/bin/perl

use strict;
use warnings;

# with equal weight
my $type = shift;
my $cgfile = shift;
my $srfile = shift; 
my $linkrfile = shift;
my $lrfile = shift;
my $rp = shift;

my %hs_pairs = ();

my %cs_scores = ();
my %cs_labels = ();

if($cgfile ne "NA"){
	read_scores($cgfile, \%cs_scores, \%cs_labels, \%hs_pairs);
}

my %sr_scores = ();
my %sr_labels = ();
if ($type =~ /SR/) {
    read_scores($srfile, \%sr_scores, \%sr_labels, \%hs_pairs);
}

my %linkr_scores = ();
my %linkr_labels = ();
if ($type =~ /LnR/) {
    read_scores($linkrfile, \%linkr_scores, \%linkr_labels, \%hs_pairs);
}

my %lr_scores = ();
my %lr_labels = ();
if ($type =~ /LR/) {
    read_scores($lrfile, \%lr_scores, \%lr_labels, \%hs_pairs);
}

my @ar_types = split(/:/, $type);
my $total = scalar(@ar_types);

my %used_pairs = ();
foreach my $key (sort keys %hs_pairs) {
    if (defined($used_pairs{$key})) { next; }
    $used_pairs{$key} = 1;
    my ($bid1, $bid2) = split(/\s+/, $key);
    my $rkey = -1*$bid2 . " " . -1*$bid1;
    $used_pairs{$rkey} = 1;

    my $sum = 0;
    my $cscore = $cs_scores{$key};
    my $clabel = $cs_labels{$key};
    
    my $sr_score = $sr_scores{$key};
    my $sr_label = $sr_labels{$key};
    
    my $linkr_score = $linkr_scores{$key};
    my $linkr_label = $linkr_labels{$key};

    my $lr_score = $lr_scores{$key};
    my $lr_label = $lr_labels{$key};

    if (defined($cscore)) { $sum += $cscore; }
    if (defined($sr_score)) { $sum += $sr_score; }
    if (defined($linkr_score)) { $sum += $linkr_score; }
    if (defined($lr_score)) { $sum += $lr_score; }
    
    my $final_score = $sum / $total;

    if (!defined($sr_label)) { $sr_label = "RN"; }
    if (!defined($linkr_label)) { $linkr_label = "RN"; }
    if (!defined($lr_label)) { $lr_label = "RN"; }

    if ($rp && ($sr_label eq "RS" || $sr_label eq "RM" || $sr_label eq "RW" ||
        $lr_label eq "RS" || $lr_label eq "RM" || $lr_label eq "RW" ||
		$linkr_label eq "RS" || $linkr_label eq "RM" || $linkr_label eq "RW")) { $final_score += 1.0; }
    
    print "$key\t$final_score\n"; 
}

########################################################################
sub read_scores {
    my $f = shift;
    my $rhs_scores = shift;
    my $rhs_labels = shift;
    my $rhs_pairs = shift;
   
    open(F, $f);
    while(<F>) {
        chomp;
        if ($_ =~ /^#/) { next; }

        my ($bid1, $bid2, $score, $label) = split(/\s+/);
        my $key = "$bid1 $bid2";
        my $rkey = -1*$bid2 . " " . -1*$bid1;
        $$rhs_scores{$key} = $score;
        $$rhs_scores{$rkey} = $score;
        $$rhs_labels{$key} = $label;
        $$rhs_labels{$rkey} = $label;

        $$rhs_pairs{$key} = 1;
        $$rhs_pairs{$rkey} = 1;
    }
    close(F);
}
