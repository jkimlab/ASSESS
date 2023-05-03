#!/usr/bin/perl

use strict;
use warnings;

my $refine_bidout_f = shift;
my $fill_f = shift;

my %new_bid = ();
my @anchor = ();
open (RF, $fill_f);
while (<RF>){
	chomp;

	if($_ =~ /^#/){
		$_ =~ s/#//;
		@anchor = split(/\s+/);
	}else{
		my @bid = split(/\s+/);
		
		if(scalar @bid >= 3){
			shift @bid;
			pop @bid;
			$new_bid{$anchor[0]}{$anchor[1]} = \@bid;
		}
	}
}
close RF;

my %update_chr = ();
my $chr_id = "";
my %end = ();
open(RF, $refine_bidout_f);
while (<RF>){
	chomp;

	if($_ =~ /^#/){
		$chr_id = substr($_, 1);
	}else{
		my @ori_bid = split(/\s+/);
		my @upd_bid = ( $ori_bid[0] );

		for(my $i=1; $i<$#ori_bid; ++$i){
			if(exists $new_bid{$ori_bid[$i-1]}{$ori_bid[$i]}){
				push (@upd_bid, @{$new_bid{$ori_bid[$i-1]}{$ori_bid[$i]}});
			}
			push (@upd_bid, $ori_bid[$i]);
		}

		$update_chr{$chr_id} = \@upd_bid;
	}
}
close RF;

my %proc = ();
foreach my $chr (sort keys %update_chr){
	if(exists $proc{$chr}){ next; }
	
	my $ext_forward= 1;
	while ($ext_forward){
		my $fbid = $update_chr{$chr}[0];
		my $rfbid;
		if($fbid =~ /([-]*)(\S+)/){
			$rfbid = ($1 ne "-")? "-".$2 : $2;
		}

		if(! exists $new_bid{$rfbid}{"\$"}){
			$ext_forward = 0;
		}else{
			$ext_forward = 1;
		
			my @bid = @{$new_bid{$rfbid}{"\$"}};
			my @rbid = Get_reverse (\@bid);
			
			push(@rbid, @{$update_chr{$chr}});
			$update_chr{$chr} = \@rbid;
			
			# Check the raca chromosomes can be connected 
			$fbid = $update_chr{$chr}[0];
			if($fbid =~ /([-]*)(\S+)/){
				$rfbid = ($1 ne "-")? "-".$2 : $2;
			}

			foreach my $cur_chr (sort keys %update_chr){
				if(exists $proc{$cur_chr}){ next; }

				my @cur_bid = @{$update_chr{$cur_chr}};
				my @rcur_bid = Get_reverse (\@cur_bid);

				my $cur_lbid = $cur_bid[-1];
				if($cur_lbid eq $fbid){
					$proc{$cur_chr} = 1;
					shift @{$update_chr{$chr}};
				
					@{$update_chr{$chr}} = (@cur_bid, @{$update_chr{$chr}});
					last;
				}
			
				my $rcur_lbid = $rcur_bid[-1];
				if($rcur_lbid eq $fbid){
					$proc{$cur_chr} = 1;
					shift @{$update_chr{$chr}};
					
					@{$update_chr{$chr}} = (@rcur_bid, @{$update_chr{$chr}});
					last;
				}
			}
		}
	}
	
	my $ext_backward= 1;
	while ($ext_backward){
		my $lbid = ${$update_chr{$chr}}[-1];
		
		if(! exists $new_bid{$lbid}{"\$"}){
			$ext_backward = 0;
		}else{
			$ext_backward = 1;

			my @bid = @{$new_bid{$lbid}{"\$"}};
			@{$update_chr{$chr}} = (@{$update_chr{$chr}}, @bid);
			
			# Check the raca chromosomes can be connected 
			$lbid = $update_chr{$chr}[-1];
			my $rlbid;
			if($lbid =~ /([-]*)(\S+)/){
				$rlbid = ($1 ne "-")? "-".$2 : $2;
			}
			
			foreach my $cur_chr (sort keys %update_chr){
				if(exists $proc{$cur_chr}){ next; }
				
				my @cur_bid = @{$update_chr{$cur_chr}};
				my @rcur_bid = Get_reverse (\@cur_bid);
				
				my $cur_fbid = $cur_bid[0];
				if($cur_fbid eq $lbid){
					$proc{$cur_chr} = 1;
					pop @{$update_chr{$chr}};
				
					@{$update_chr{$chr}} = (@{$update_chr{$chr}}, @cur_bid);
					last;
				}
				
				my $rcur_fbid = $rcur_bid[0];
				if($rcur_fbid eq $lbid){
					$proc{$cur_chr} = 1;
					pop @{$update_chr{$chr}};
					
					@{$update_chr{$chr}} = (@{$update_chr{$chr}}, @rcur_bid);
					last;
				}
			}
		}
	}
}

my $new_chr = 1;
foreach my $chr (sort {$a <=> $b} keys %update_chr){
	if($proc{$chr}){ next; } 
	print "#$new_chr\n@{$update_chr{$chr}}\t\$\n";
	
	$new_chr++;
}
	
sub Get_reverse {
	my @_bid = @{$_[0]};

	my @_rbid = ();
	foreach my $id ( reverse @_bid ){
		my $rid;
		if($id =~ /([-]*)(\S+)/){
			$rid = ($1 ne "-")? "-".$2 : $2;
		}
		push(@_rbid, $rid);
	}
	
	return @_rbid;
}
