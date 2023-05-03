#!/usr/bin/perl

use strict;
use warnings;

# input
my $refspc = shift;	# reference species
my $maf_f = shift;	# rec_chrs.refspc.segments.refined.txt

open(F, $maf_f);
my @lines = <F>;
close(F);

my %raca_chrs = ();
my %spc_chrs = ();
my %spc_dirs = ();
my %raca_dirs = ();
my %raca2spc = ();
my %spc2raca = ();
for (my $i = 0; $i < $#lines; $i+=4) {
	my $raca_line = $lines[$i+1];	
	my $spc_line = $lines[$i+2];	

	$raca_line =~ /\S+\.(\S+):(\S+)\-(\S+) (\S)/;
	my ($rchr, $rstart, $rend, $rdir) = ($1, $2, $3, $4);
	$raca_chrs{$rchr}{$rstart} = $rend;
	$raca_dirs{"$rchr:$rstart"} = $rdir;

	$spc_line =~ /\S+\.(\S+):(\S+)\-(\S+) (\S)/;
	my ($schr, $sstart, $send, $sdir) = ($1, $2, $3, $4);
	$spc_chrs{$schr}{$sstart} = $send;
	$spc_dirs{"$schr:$sstart"} = $sdir;
	
	$raca2spc{"$rchr:$rstart"} = "$schr:$sstart";
	$spc2raca{"$schr:$sstart"} = "$rchr:$rstart";
}

my %raca_key2id = (); 
make_key2id(\%raca_chrs, \%raca_key2id);
my %spc_key2id = (); 
make_key2id(\%spc_chrs, \%spc_key2id);

print "#RACA centric configuration\n";
print_conf(\%raca_chrs, \%raca2spc, \%spc_key2id, \%raca_dirs, \%spc_dirs, 1);
print "\n";
print "#$refspc centric configuration\n";
print_conf(\%spc_chrs, \%spc2raca, \%raca_key2id, \%spc_dirs, \%raca_dirs, 0);

sub print_conf {
	my $ref_main_chrs = shift;
	my $ref_main2tar = shift;
	my $ref_tarkey2id = shift;
	my $ref_main_dirs = shift;
	my $ref_tar_dirs = shift;
	my $sortflag = shift;

	my @chrs = sort keys %$ref_main_chrs;
	if ($sortflag == 1) {
		@chrs = sort {$a<=>$b} keys %$ref_main_chrs;
	}

	my $num_inter = 0;
	my $num_intra = 0;
    my $num_intrax = 0;
	foreach my $chr (@chrs) {
		my $rhs = $$ref_main_chrs{$chr};
		my $oline = "$chr\t";
		my %tarchrs = ();
		my $prev_tchr = "";
		foreach my $start (sort {$a<=>$b} keys %$rhs) {
			my $mainkey = "$chr:$start";
			my $tarkey = $$ref_main2tar{$mainkey};
			my $maindir = $$ref_main_dirs{$mainkey};
			my $tardir = $$ref_tar_dirs{$tarkey};
			my $tarid = $$ref_tarkey2id{$tarkey};
            
            my $end = $$rhs{$start};
            my $len = $end - $start;

			my ($tchr, $tstart) = split(/:/, $tarkey);
			$tarchrs{$tchr} = 1;
			if ($prev_tchr ne "" && $prev_tchr ne $tchr) {
				$num_inter++;
			} elsif ($prev_tchr eq $tchr) {
				$num_intra++;
                if ($tchr =~ /chrX/ || $tchr =~ /chrx/) { $num_intrax++; }
			}

			if ($maindir eq $tardir) {
				$oline .= "$tarid $len ";
			} else {
				$oline .= "-$tarid $len ";
			}

			$prev_tchr = $tchr;
		}

		if (scalar(keys %tarchrs) > 1) {
			print "INTER\t$oline\n";
		} else {
			print "INTRA\t$oline\n";
		}
	}
	print "Total inter: $num_inter\n";
	print "Total intra: $num_intra\n";
	print "Total intrax: $num_intrax\n";
}

sub make_key2id {
	my $ref_chrs = shift;
	my $ref_out = shift;

	foreach my $chr (sort keys %$ref_chrs) {
		my $rhs = $$ref_chrs{$chr};
		my @starts = sort {$a<=>$b} keys %$rhs;
		if (scalar(@starts) == 1) {
			my $key = "$chr:$starts[0]";
			$$ref_out{$key} = $chr;
		} else {
			my $cnt = 0;
			foreach my $start (@starts) {
				my $key = "$chr:$start";
				my $new_chr = $chr . "_$cnt";
				$$ref_out{$key} = $new_chr;
				$cnt++;
			}
		}
	}
}


