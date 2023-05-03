#!/usr/bin/per/
use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';
use File::Basename;
use Bio::TreeIO;
use IO::String;

#my $in_params = $ARGV[0];
#my $in_cpu = $ARGV[1];
#my $out_dir = $ARGV[-1];

## parameters
my $lib_param_f;
my $run_param_f;
my $reference_manager;
my $help;
## setting default parameters
my $in_cpu = 1;
my $resolution = 1000;

## get arguments
my $options = GetOptions(
				"lib=s" => \$lib_param_f,
				"params=s" => \$run_param_f,
				"p|process=i" => \$in_cpu,
				"h|help" => \$help
);

## checking mandatory arguments
if(!defined($run_param_f)){
		print STDERR "# running parameter file is mandatory, please check the -params option \n";
		exit();
}
if(!defined($lib_param_f)){
		print STDERR "# lib parameter file is mandatory, please check the -lib option \n";
		exit();
}
if(defined($help)){
		print STDERR "# Usage: making_cmd.pl -params running_option_file -lib input_option_file > run_cmd.sh\n";
		exit();
}
print '#!/bin/bash'."\n";

## global variables
my %hs_lib_param = ();
my %hs_sequencing_read = ();
my %hs_run_param = ();
my $param_name;
my $initasm_in_container = "/assess_wd/final_out/initial_assembly";
my $mapping_in_container = "/assess_wd/final_out/mapping";
my $initasm_user_out = "";
my %hs_binding_path = ();
my $binding_path = "";
my $tmp_binding_path;

## reading input library files 
open(FLPA, "$lib_param_f");
while(<FLPA>){
		chomp;
		if($_ eq ""){ next; }
		if($_ =~ /^#/){ next; } ## commentary
		if($_ =~ /^>(\S+)/){
				$param_name = $1;
		}
		else{
				my @arr_input = split(/\s+/,$_);
				if($param_name eq "LIB"){
						if(scalar(@arr_input) != 3){ next; }
						# $arr_input[0] : SR or LR or HIC or 10X
						# $arr_input[1] : prefix of the library
						# $arr_input[2] : file path
						$hs_lib_param{$arr_input[0]}{$arr_input[1]} = $arr_input[2];
				}
				elsif($param_name eq "SR_INS"){
						if(scalar(@arr_input) != 3){ next; }
						$hs_lib_param{"SR_INS"}{$arr_input[0]} = "$arr_input[1]:$arr_input[2]";
				}
				else{
						$hs_lib_param{$param_name}{$arr_input[0]} = $arr_input[1];
				}
		}
}
close(FLPA);
## checking list of sequence library files
print STDERR "# ---+ recording lib parameters (not checking the input files are in proper position)\n";
sub CHECK_PARAM_DEFINED{ 
		my $param_for_check = shift(@_);
		if(defined($param_for_check)){ print STDERR "# yes\n"; }
		else{ print STDERR "# no use\n"; }
		print STDERR "\n";
}
print STDERR "# \t---+ reference manager ... ";
CHECK_PARAM_DEFINED($hs_lib_param{"REF"}{"REF_DIR"});
print STDERR "# \t---+ initial assembly ... ";
CHECK_PARAM_DEFINED($hs_lib_param{"ASSEMBLY"}{"SEQFILE"});
print STDERR "# \t---+ short read mapping directory... ";
CHECK_PARAM_DEFINED($hs_lib_param{"MAPPING"}{"SR_MAPPING_DIR"});
print STDERR "# \t---+ long read mapping directory... ";
CHECK_PARAM_DEFINED($hs_lib_param{"MAPPING"}{"LR_MAPPING_DIR"});
print STDERR "# \t---+ hi-c read mapping directory... ";
CHECK_PARAM_DEFINED($hs_lib_param{"MAPPING"}{"HIC_MAPPING_DIR"});
print STDERR "# \t---+ 10x linked read mapping directory... ";
CHECK_PARAM_DEFINED($hs_lib_param{"MAPPING"}{"10X_MAPPING_DIR"});
foreach my $seq_type ("SR", "HIC", "10X"){
		print STDERR "# \t---+ $seq_type sequencing files...\n";
		foreach my $lib_num (sort keys(%{$hs_lib_param{$seq_type}})){
				print "# \t\t---+ $lib_num... ";
				CHECK_PARAM_DEFINED($hs_lib_param{$seq_type}{$lib_num});
				print "\n";
		}
}
print "\n";
print STDERR "# \t---+ LR sequencing files... ";
foreach my $lib_num (sort keys(%{$hs_lib_param{LR}})){
		print "# \t\t---+ $lib_num... ";
		CHECK_PARAM_DEFINED($hs_lib_param{LR}{$lib_num});
		print "\n";
}
print "\n";

## command printing 
print 'OUTDIR=$(realpath $1)'."\n";
my $output_in_container = "/assess_wd/final_out";
$tmp_binding_path = "-v \$OUTDIR:$output_in_container";
if(!exists($hs_binding_path{$tmp_binding_path})){
		$binding_path = "-v \$OUTDIR:$output_in_container";
}
$param_name = ();

## input ASSESS running parameters

print STDERR "# ---+ creating params.txt for running ASSESS\n";
print STDERR "# \t---+ data directory for running ASSESS: \$OUTDIR/data\n";
print STDERR "# \t---+ params.txt file is reset now\n";
print 'mkdir -p $OUTDIR/data'."\n"; # !! /assess_wd/final_out/data is created !!
print 'rm -f $OUTDIR/data/params.txt'."\n"; # !! the params file already existed in /assess_wd/final_out/data is removed !!
open(FRPA, "$run_param_f");
while(<FRPA>){
		chomp;
		if($_ =~ /^#/){ next; } ## commentary
		my ($indicator, $value) = split(/=/,$_);
		if($indicator eq "TREE"){
				print "echo \"$value\" > \$OUTDIR/data/tree.nwk\n";
				print "echo \"$_\" >> \$OUTDIR/data/params.txt\n";
		}
		else{ print "echo \"$_\" >> \$OUTDIR/data/params.txt\n"; }
		if($indicator =~ /^DIVTIME_(\S+)/){
				$hs_run_param{"$1"} = lc($value);
				next;
		}
		$hs_run_param{$indicator} = $value;

}
close(FRPA);
#print STDERR "# ---+ params.txt for running ASSESS has been created\n";
## Main works (assembly, masking, read mapping to assembly, whole-genome assembly for references, ASSESS)
### 1. Setting target assembly 
print STDERR "# ---+ MAIN 1. initial assembly\n";
if(defined($hs_lib_param{'ASSEMBLY'}{'SEQFILE'})){ ## initial assembly (O)
		print STDERR "# \t---+ initial assembly is already defined (skipped)\n";
		my @file_path = split(/\//,$hs_lib_param{'ASSEMBLY'}{'SEQFILE'});
		my $file_name = pop(@file_path);
		$initasm_user_out = join("/",@file_path);
		$tmp_binding_path .= " -v $initasm_user_out:$initasm_in_container";
		if(!exists($hs_binding_path{$tmp_binding_path})){
				$binding_path .= " -v $initasm_user_out:$initasm_in_container";
		}
		$hs_lib_param{'ASSEMBLY'}{'SEQFILE'} = "$initasm_in_container/$file_name"; ## initial assembly file is updated.
}
else{
		print STDERR "# \t---+ initial assembly will be generated\n";
		$hs_run_param{ASSEMBLER} = lc($hs_run_param{ASSEMBLER});
		print STDERR "# \t---+ initial assembly needs to be created.\n# \t---+ assembly will be performed $hs_run_param{ASSEMBLER}\n";
		$initasm_user_out = "\$OUTDIR/initial_assembly";
		$tmp_binding_path .= " -v $initasm_user_out:$initasm_in_container";
		if(!exists($hs_binding_path{$tmp_binding_path})){
				$binding_path .= " -v $initasm_user_out:$initasm_in_container";
		}
		print STDERR "# \t---+ initial assembly directory: $initasm_user_out\n";
		my $tmp_lib_options = "";
		foreach my $seq_type ("SR", "LR"){
				foreach my $num_lib (sort keys(%{$hs_lib_param{$seq_type}})){
						if(exists($hs_lib_param{$seq_type}{$num_lib})){
								my @file_path = split(/\//,$hs_lib_param{$seq_type}{$num_lib});
								my $file_name = pop(@file_path);
								my $original_lib_dir = join("/",@file_path);
								$tmp_binding_path .= " -v $original_lib_dir:/assess_wd/lib/$seq_type$num_lib";
								if(!exists($hs_binding_path{$tmp_binding_path})){
										$binding_path .= " -v $original_lib_dir:/assess_wd/lib/$seq_type$num_lib";
								}
								$tmp_lib_options .= " $seq_type:$num_lib:/assess_wd/lib/$seq_type$num_lib/$file_name";
								if($seq_type eq "SR"){
										my ($lib_id, $fr) = split(/-/,$num_lib);
										if(exists($hs_lib_param{"SR_INS"}{$lib_id})){ $tmp_lib_options .= ":$hs_lib_param{SR_INS}{$lib_id}"; }
										else{ $tmp_lib_options .= ":500:50"; }
								}
								else{
										$tmp_lib_options .= ":$hs_lib_param{LR_TYPE}{$num_lib}";
								}
						}
				}
		}
		if(defined($hs_lib_param{GENOMEINFO}{GENOMESIZE})){
				$tmp_lib_options .= " -g $hs_lib_param{GENOMEINFO}{GENOMESIZE}";
		}
		my $initasm_cmd = "docker run --rm $binding_path assessw:latest do_assembly_$hs_run_param{ASSEMBLER}.pl $initasm_in_container";
		print "$initasm_cmd"."$tmp_lib_options\n";
		$hs_lib_param{'ASSEMBLY'}{'SEQFILE'} = "$initasm_in_container/scaf.fa";
}

### 2. masking target assembly

print STDERR "# \t---+ masking initial assembly\n";
print "docker run --rm $binding_path assessw:latest check_scaf_masking.pl $hs_lib_param{'ASSEMBLY'}{'SEQFILE'} $in_cpu $hs_run_param{'TARSPC'} $hs_run_param{'TARSPC'} $output_in_container/data >> \$OUTDIR/data/params.txt\n";
$hs_lib_param{'ASSEMBLY'}{'SEQFILE'} = "$output_in_container/data/$hs_run_param{TARSPC}.fa";

### 3. mapping all reads to target assembly
#### 3.1 Short read mapping
my $mapping_bind = "";
if(defined($hs_lib_param{"MAPPING"}{"SR_MAPPING_DIR"})){
		$binding_path .= " -v $hs_lib_param{MAPPING}{SR_MAPPING_DIR}:$mapping_in_container/mapping_SR";
		print "echo \"SR_MAPPING_DIR=$mapping_in_container/mapping_SR\" >> \$OUTDIR/data/params.txt\n";
}
else{
		print "mkdir -p \$OUTDIR/data/mapping/mapping_SR\n";
		$mapping_bind = " -v \$OUTDIR/data/mapping/mapping_SR:$mapping_in_container/mapping_SR";
}
my $flag = 0;
foreach my $num_lib (sort keys(%{$hs_lib_param{SR}})){
		if($flag == 0){
				print "echo \"SR_MAPPING_DIR=$mapping_in_container/mapping_SR\" >> \$OUTDIR/data/params.txt\n";
				$flag = 1;
		}
				
		my ($lib_num, $frag) = split(/-/,$num_lib);
		if($frag == 1){
				my @file_path = split(/\//,$hs_lib_param{SR}{$num_lib}); ## eg. /mss1/MSIP/SITRI/clustering/evaluations/yeast/data/S288C_simulated1.fq.gz
				my $file_name = pop(@file_path); ## eg. S288C_simulated1.fq.gz
				my $original_lib_dir = join("/",@file_path); ## eg. /mss1/MSIP/SITRI/clustering/evaluations/yeast/data
				my $tmp_bind = " -v $original_lib_dir:/assess_wd/lib/SR$num_lib"; ## eg. /mss1/MSIP/SITRI/clustering/evaluations/yeast/data == /assess_wd/lib/SR
				my $additional_pair = $hs_lib_param{SR}{"$lib_num-2"};
				@file_path = split(/\//,$additional_pair); ## eg. /mss1/MSIP/SITRI/clustering/evaluations/yeast/data/S288C_simulated1.fq.gz
				$additional_pair = pop(@file_path); ## eg. S288C_simulated1.fq.gz
				$original_lib_dir = join("/",@file_path); ## eg. /mss1/MSIP/SITRI/clustering/evaluations/yeast/data
				$tmp_bind .= " -v $original_lib_dir:/assess_wd/lib/SR$lib_num-2"; ## eg. /mss1/MSIP/SITRI/clustering/evaluations/yeast/data == /assess_wd/lib/SR
				print "docker run --rm $binding_path $tmp_bind assessw:latest do_mapping_bwa.pl SR-$lib_num $output_in_container/data/$hs_run_param{TARSPC}.fa  $mapping_in_container/mapping_SR /assess_wd/lib/SR$num_lib/$file_name /assess_wd/lib/SR$lib_num-2/$additional_pair $in_cpu $output_in_container/data\n";
		}
}

#### 3.2 Long read mapping
$mapping_bind = "";
if(defined($hs_lib_param{"MAPPING"}{"LR_MAPPING_DIR"})){
		$binding_path .= " -v $hs_lib_param{MAPPING}{LR_MAPPING_DIR}:$mapping_in_container/mapping_LR";
		print "echo \"LR_MAPPING_DIR=$mapping_in_container/mapping_LR\" >> \$OUTDIR/data/params.txt\n";
}
else{
		print "mkdir -p \$OUTDIR/data/mapping/mapping_LR\n";
		$mapping_bind = " -v \$OUTDIR/data/mapping/mapping_LR:$mapping_in_container/mapping_LR";
}
$flag = 0;
foreach my $num_lib (sort keys(%{$hs_lib_param{LR}})){
		if($flag == 0){
				print "echo \"LR_MAPPING_DIR=$mapping_in_container/mapping_LR\" >> \$OUTDIR/data/params.txt\n";
				$flag = 1;
		}
		my @file_path = split(/\//,$hs_lib_param{LR}{$num_lib}); ## eg. /mss1/MSIP/SITRI/clustering/evaluations/yeast/data/S288C_simulated1.fq.gz
		my $file_name = pop(@file_path); ## eg. S288C_simulated1.fq.gz
		my $original_lib_dir = join("/",@file_path); ## eg. /mss1/MSIP/SITRI/clustering/evaluations/yeast/data
		my $tmp_bind = " -v $original_lib_dir:/assess_wd/lib/LR$num_lib"; ## eg. /mss1/MSIP/SITRI/clustering/evaluations/yeast/data == /assess_wd/lib/SR
		print "docker run --rm $binding_path $tmp_bind assessw:latest do_mapping_minimap2.pl LR-$num_lib $output_in_container/data/$hs_run_param{TARSPC}.fa $mapping_in_container/mapping_LR /assess_wd/lib/LR$num_lib/$file_name $hs_lib_param{LR_TYPE}{$num_lib} $in_cpu $output_in_container/data\n";
}

#### 3.3 Hi-C read mapping
$mapping_bind = "";
if(defined($hs_lib_param{"MAPPING"}{"HIC_MAPPING_DIR"})){
		$binding_path .= " -v $hs_lib_param{MAPPING}{HIC_MAPPING_DIR}:$mapping_in_container/mapping_HIC";
		print "echo \"HIC_MAPPING_DIR=$mapping_in_container/mapping_HIC/bam\" >> \$OUTDIR/data/params.txt\n";
}
else{
		print "mkdir -p \$OUTDIR/data/mapping/mapping_HIC\n";
		$mapping_bind = " -v \$OUTDIR/data/mapping/mapping_HIC:$mapping_in_container/mapping_HIC";
}
$flag = 0;
foreach my $num_lib (sort keys(%{$hs_lib_param{HIC}})){
		if($flag == 0){
				print "echo \"HIC_MAPPING_DIR=$mapping_in_container/mapping_HIC/bam\" >> \$OUTDIR/data/params.txt\n";
				$flag = 1;
		}
		my ($lib_num, $frag) = split(/-/,$num_lib);
		if($frag == 1){
				my @file_path = split(/\//,$hs_lib_param{HIC}{$num_lib}); 
				my $file_name = pop(@file_path); 
				my $original_lib_dir = join("/",@file_path); 
				my $tmp_bind = " -v $original_lib_dir:/assess_wd/lib/HIC$num_lib"; 
				my $additional_pair = $hs_lib_param{HIC}{"$lib_num-2"};
				@file_path = split(/\//,$additional_pair); 
				$additional_pair = pop(@file_path); 
				$original_lib_dir = join("/",@file_path); 
				$tmp_bind .= " -v $original_lib_dir:/assess_wd/lib/HIC$lib_num-2"; 
				print "docker run --rm $binding_path $tmp_bind assessw:latest do_mapping_hic.pl HIC-$lib_num $output_in_container/data/$hs_run_param{TARSPC}.fa  $mapping_in_container/mapping_HIC /assess_wd/lib/HIC$num_lib/$file_name /assess_wd/lib/HIC$lib_num-2/$additional_pair $in_cpu $output_in_container/data\n";
		}
}
#### 3.4 10X read mapping
$mapping_bind = "";
if(defined($hs_lib_param{"MAPPING"}{"10X_MAPPING_DIR"})){
		$binding_path .= " -v $hs_lib_param{MAPPING}{'10X_MAPPING_DIR'}:$mapping_in_container/mapping_10X";
		print "echo \"LINKR_MAPPING_DIR=$mapping_in_container/mapping_10X/bam\" >> \$OUTDIR/data/params.txt\n";
}
else{
		print "mkdir -p \$OUTDIR/data/mapping/mapping_10X\n";
		$mapping_bind = " -v \$OUTDIR/data/mapping/mapping_10X:$mapping_in_container/mapping_10X";
}
$flag = 0;
foreach my $num_lib (sort keys(%{$hs_lib_param{'10X'}})){
		if($flag == 0){
				print "echo \"LINKR_MAPPING_DIR=$mapping_in_container/mapping_10X/bam\" >> \$OUTDIR/data/params.txt\n";
				$flag = 1;
		}
		my ($lib_num, $frag) = split(/-/,$num_lib);
		if($frag == 1){
				my @file_path = split(/\//,$hs_lib_param{'10X'}{$num_lib}); 
				my $file_name = pop(@file_path); 
				my $original_lib_dir = join("/",@file_path); 
				my $tmp_bind = " -v $original_lib_dir:/assess_wd/lib/10X$num_lib"; 
				my $additional_pair = $hs_lib_param{'10X'}{"$lib_num-2"};
				@file_path = split(/\//,$additional_pair); 
				$additional_pair = pop(@file_path); 
				$original_lib_dir = join("/",@file_path); 
				$tmp_bind .= " -v $original_lib_dir:/assess_wd/lib/10X$lib_num-2"; 
				print "docker run --rm $binding_path $tmp_bind assessw:latest do_mapping_10x.pl LINK$lib_num $output_in_container/data/$hs_run_param{TARSPC}.fa  $mapping_in_container/mapping_10X /assess_wd/lib/10X$num_lib/$file_name /assess_wd/lib/10X$lib_num-2/$additional_pair $in_cpu $output_in_container/data\n";
		}
}


### 4. Reference preparation
my $ref_in_container = "";
if(defined($hs_lib_param{'REF'}{'REF_DIR'})){ # !! /assess_wd/REF is created by linking with $reference_manager
		$reference_manager = $hs_lib_param{'REF'}{'REF_DIR'};
		$tmp_binding_path .= " -v $reference_manager:/assess_wd/REF";
		if(!exists($hs_binding_path{$tmp_binding_path})){
				$binding_path .= " -v $reference_manager:/assess_wd/REF";
		}
		print "mkdir -p $reference_manager/fasta $reference_manager/list $reference_manager/chainNet $reference_manager/temp\n";
		print "echo \"CHAINNET_DIR=/assess_wd/REF/chainNet\" >> \$OUTDIR/data/params.txt\n";
		$ref_in_container = "/assess_wd/REF";
}
else{ # !! /assess_wd/final_out/REF is newly created 
		print "mkdir -p \$OUTDIR/REF/fasta \$OUTDIR/REF/list \$OUTDIR/REF/chainNet \$OUTDIR/REF/temp\n";
		print "echo \"CHAINNET_DIR=$output_in_container/REF/chainNet\" >> \$OUTDIR/data/params.txt\n";
		$ref_in_container = "$output_in_container/REF";
		$reference_manager = "\$OUTDIR/REF";
}

my $io = IO::String->new($hs_run_param{TREE});
my $treeio = Bio::TreeIO->new(-fh => $io,
                              -format => 'newick');
my $tree = $treeio->next_tree;
my $rootnode = $tree->get_root_node;
# reference genome sequence preparation ==> need to be remake
foreach my $node ( $rootnode->get_all_Descendents() ) { 
    if( $node->is_Leaf ) { 
      my $n_id = $node->id;
				if($n_id =~ /\@/){
						next;
				}
				print "docker run --rm $binding_path assessw:latest genome_preparation.pl $n_id $ref_in_container/fasta\n";
    }   
}

########## whole-genome alignment preparation
foreach my $node ( $rootnode->get_all_Descendents() ) { 
    if( $node->is_Leaf ) { 
      	my $n_id = $node->id;
				if($n_id =~ /\@/){
						next;
				}
				if($n_id eq $hs_run_param{'REFSPC'}){ next; }
				if(!-f "$ref_in_container/list/$hs_run_param{'REFSPC'}/$n_id.done"){
						print "docker run --rm $binding_path assessw:latest buildChainNet_forkManager.pl -p $in_cpu --divtime $hs_run_param{$n_id} -o $ref_in_container/chainNet/$hs_run_param{'REFSPC'}/$n_id/raw_chainNet -chainNet_dir $ref_in_container/chainNet/$hs_run_param{'REFSPC'}/$n_id $ref_in_container/fasta/$hs_run_param{'REFSPC'}.fa.masked.gz $ref_in_container/fasta/$n_id.fa.masked.gz\n"; 
#						print "docker run --rm $binding_path assessw:latest touch $reference_manager/list/$hs_run_param{'REFSPC'}/$n_id.done\n";
				}
    }   
}
my $temp_tar = "D"."$hs_run_param{TARSPC}";
print "docker run --rm $binding_path assessw:latest buildChainNet_forkManager.pl -p $in_cpu --divtime $hs_run_param{$hs_run_param{'TARSPC'}} -o $ref_in_container/chainNet/$hs_run_param{'REFSPC'}/$hs_run_param{TARSPC}/raw_chainNet -chainNet_dir $ref_in_container/chainNet/$hs_run_param{'REFSPC'}/$hs_run_param{TARSPC} $ref_in_container/fasta/$hs_run_param{REFSPC}.fa.masked.gz $output_in_container/data/$hs_run_param{TARSPC}.fa\n";
#}
#=cut
## Run ASSESS
print "docker run --rm $binding_path assessw:latest /ASSESS/ASSESS.b12/ASSESS -p /assess_wd/final_out/data/params.txt -r $hs_run_param{'RESOLUTION'} -o /assess_wd/final_out/assess_out\n";



########## Finishing
if(defined($reference_manager)){
		print "docker run --rm -v $reference_manager:$ref_in_container assessw:latest chmod 777 -R $ref_in_container\n";
}
print "docker run --rm -v \$OUTDIR:$output_in_container assessw:latest chmod 777 -R /assess_wd/final_out\n";
