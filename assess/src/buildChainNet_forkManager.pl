#!/usr/bin/env perl

use strict;
use warnings;
use Cwd;
use Cwd 'abs_path';
use File::Basename;
use Parallel::ForkManager;
use Getopt::Long qw(:config no_ignore_case);
use FindBin '$Bin';

my $lastz_E;
my $lastz_H;
my $lastz_K;
my $lastz_L;
my $lastz_M;
my $lastz_O;
my $lastz_Q;
my $lastz_T;
my $lastz_Y;
my $chain_minScore;
my $chain_linearGap;
my $divtime;    # near, medium (human vs mouse), or far (human vs chicken), custom
my $target_drop = 1; # 1 or 0
my $out_dir = getcwd;
my $help;
my $core = 1;
my $lastz_bit = 1;
my $chainNet_dir = "";

GetOptions (
	"core|p=i" => \$core,
        "divtime|D=s" => \$divtime,
        "E=i" => \$lastz_E,
        "K=i" => \$lastz_K,
        "L=i" => \$lastz_L,
        "M=i" => \$lastz_M,
        "O=i" => \$lastz_O,
        "T=i" => \$lastz_T,
        "Y=i" => \$lastz_Y,
        "Q=s" => \$lastz_Q,
        "minScore|S=i" => \$chain_minScore,
        "linearGap|G=s" => \$chain_linearGap,
				"chainNet_dir=s" => \$chainNet_dir,
        "NonTargetDrop|ND:i" => \$target_drop,
        "outDir|o=s" => \$out_dir,
        "help|h" => \$help,
        "lastz_32:i" => \$lastz_bit,
);

if(!$divtime||$help||$#ARGV != 1){
	print "buildChainNet_forkManager.pl --divtime [divtime] -p [core_num] -o [out_dir] [Ref(.fasta)] [Tar(.fasta)]\n";
        exit(1);
}

if($divtime eq "custom"){
        if(!$lastz_E || !$lastz_K || !$lastz_L || !$lastz_M || !$lastz_O || !$lastz_T || !$lastz_Y || !$lastz_Q || !$chain_minScore || !$chain_linearGap){
                exit(1);
        }
}

print STDERR "# ---# Whole-genome alignment start\n";
if((-d "$chainNet_dir/chain") && (-d "$chainNet_dir/net")){
		print STDERR "# \t---# chainNet files are already prepared (skipped)\n";
		exit();
}
my $pm = new Parallel::ForkManager($core);
my $time = `date`;
chomp($time);
print STDERR "# \t---# [START] $time\n";
print STDERR "# \t---# setting paramters\n";
my $ref_fa = shift;
my $tar_fa = shift;

my $lastz = "lastz";

$ref_fa = abs_path($ref_fa);
$tar_fa = abs_path($tar_fa);
system("mkdir -p $out_dir");
$out_dir = abs_path($out_dir);


my $chainPar = "";
my $lastzPar = "";
if ($divtime eq "near") {
        $chainPar = "-minScore=5000 -linearGap=$Bin/params/human_chimp.linearGap";
        $lastzPar = "E=150 H=2000 K=4500 L=2200 M=254 O=600 Q=$Bin/params/human_chimp.v2.q T=2 Y=15000";
} elsif ($divtime eq "medium") {
        $chainPar = "-minScore=3000 -linearGap=medium";
        $lastzPar = "E=30 H=2000 K=3000 L=2200 M=50 O=400 Q=$Bin/params/human_mouse.v2.q T=1 Y=9400";
} elsif ($divtime eq "medium2") {
        $chainPar = "-minScore=3000 -linearGap=medium";
        $lastzPar = "E=30 H=2000 K=3000 L=2200 M=254 O=400 Q=$Bin/params/human_mouse.v2.q T=1 Y=9400";
} elsif ($divtime eq "medium3") {
        $chainPar = "-minScore=3000 -linearGap=medium";
        $lastzPar = "E=30 H=2000 K=3000 L=3000 M=50 O=400 Q=$Bin/params/human_mouse.v2.q T=1 Y=9400";
} elsif ($divtime eq "far") {
        $chainPar = "-minScore=5000 -linearGap=loose";
        $lastzPar = "E=30 H=2000 K=2200 L=6000 M=50 O=400 Q=$Bin/params/HoxD55.q T=2 Y=3400";
} elsif ($divtime eq "fish") {
        $chainPar = "-minScore=3000 -linearGap=medium";
        $lastzPar = "E=30 H=2000 K=3000 L=3000 M=50 O=400 Q=$Bin/params/HoxD55.q T=2 Y=9400";
} elsif ($divtime eq "custom") {
        $chainPar = "-minScore=$chain_minScore -linearGap=$chain_linearGap";
        $lastzPar = "E=$lastz_E H=$lastz_H K=$lastz_K L=$lastz_L M=$lastz_M O=$lastz_O Q=$lastz_Q T=$lastz_T Y=$lastz_Y";
} else {
        print STDERR "Unrecognized a divergence time parameter!!\n";
        exit(1);
}

# split fasta sequences
print STDERR "# \t---# Split fasta files...\n";
`mkdir -p $out_dir/ref`;
`mkdir -p $out_dir/tar`;
my @split_cmds = ("faSplit byName $ref_fa $out_dir/ref/","faSplit byName $tar_fa $out_dir/tar/");
foreach my $c (0..$#split_cmds)
{
	my $pid = $pm->start($split_cmds[$c]) and next;
	system("$split_cmds[$c]");
	$pm->finish($c);
}
$pm -> wait_all_children;

# create 2bit files
print STDERR "# \t---# Create 2bit files...\n";
my $bref_name = basename($ref_fa,".numID.fa",".fa",".fna",".fasta");
my $ref_2bit = "$out_dir/ref/$bref_name.2bit";
`faToTwoBit $ref_fa $ref_2bit`;
my $btar_name = basename($tar_fa,".numID.fa",".fa",".fna",".fasta");
my $tar_2bit = "$out_dir/tar/$btar_name.2bit";
`faToTwoBit $tar_fa $tar_2bit`;

print STDERR "# \t---# Create size files...\n";
my $ref_size = "$out_dir/ref/$bref_name.sizes";
my $tar_size = "$out_dir/tar/$btar_name.sizes";

my @size_cmds = ("faSize $ref_fa -detailed > $ref_size","faSize $tar_fa -detailed > $tar_size");
foreach my $c (0..$#size_cmds){
	my $pid = $pm->start($size_cmds[$c]) and next;
	system("$size_cmds[$c]");
	$pm->finish($c);
}
$pm -> wait_all_children;

my $dcnt = 0;
if($target_drop == 1){
	print STDERR "# \t---# Removing small (<1Kbp) target files...\n";
	open(F,"$ref_size");
	while(<F>) {
		chomp;
		my ($name, $size) = split(/\s+/);
		if (-f "$out_dir/ref/$name.fa" && $size < 1000) {
			`rm -f $out_dir/ref/$name.fa`; 
			$dcnt++;
		}
	}
	close(F);
	print STDERR "\t$dcnt target files were removed.\n";
}
$dcnt = 0;

print STDERR "# \t---# Merging small (<10Mbp) query files...\n";
my $merge_f = "$out_dir/tar/merged.fa";
open(F,"$tar_size");
while(<F>) {
	chomp;
	my ($name, $size) = split(/\s+/);
	if (-f "$out_dir/tar/$name.fa" && $size < 10000000) {
		`cat $out_dir/tar/$name.fa >> $merge_f`;
		`rm -f $out_dir/tar/$name.fa`; 
		$dcnt++;
	}
}
close(F);
print STDERR "\t$dcnt query files were merged to $merge_f.\n";

my @ref_fas = <$out_dir/ref/*.fa>;
my @tar_fas = <$out_dir/tar/*.fa>;

my $num_ref = scalar(@ref_fas);
my $num_tar = scalar(@tar_fas);

my $total_jobs = $num_ref * $num_tar;

my $lav_dir = "$out_dir/lav";
`mkdir -p $lav_dir`;
my $psl_dir = "$out_dir/psl";
`mkdir -p $psl_dir`;
my $log_dir = "$out_dir/log";
`mkdir -p $log_dir`;
my $done_dir = "$out_dir/done";
`mkdir -p $done_dir`;
my $chain_dir = "$out_dir/chain";
`mkdir -p $chain_dir`;

# create a job file
my $job_cnt = 1;
foreach my $ref_file (@ref_fas) {
	my $ref_name = basename($ref_file,".fa");
	foreach my $tar_file (@tar_fas) {
		my $tar_name = basename($tar_file,".fa");

		`mkdir -p $lav_dir/$ref_name`;
		`mkdir -p $psl_dir/$ref_name`;
		`mkdir -p $chain_dir/$ref_name`;
		`mkdir -p $done_dir/$ref_name`;
		my $lav_f = "$lav_dir/$ref_name/$ref_name.$tar_name.lav";
		my $psl_f = "$psl_dir/$ref_name/$ref_name.$tar_name.psl";
		my $chain_f = "$chain_dir/$ref_name/$ref_name.$tar_name.chain";
	
		`mkdir -p $out_dir/jobs`;
		my $job_f = "$out_dir/jobs/$job_cnt.job"; 
		open(O,">$job_f");
		print O "#!/usr/bin/perl\n";
		print O "use strict;\n";
		print O "use warnings;\n";
		print O "system(\"$lastz $ref_file $tar_file $lastzPar > $lav_f\");\n";
		print O "system(\"lavToPsl $lav_f $psl_f\");\n";
		print O "system(\"axtChain -psl -verbose=0 $chainPar $psl_f $ref_2bit $tar_2bit stdout | chainAntiRepeat $ref_2bit $tar_2bit stdin $chain_f\");\n"; 
		print O "system(\"echo $ref_name vs $tar_name | cat > $done_dir/$ref_name/$job_cnt.done\");\n";
		close(O);
		$job_cnt++;
	}
}

# wait for all job completion
my @fork_jobs = <$out_dir/jobs/*.job>;
my @align_cmds;
foreach my $job_f (@fork_jobs)
{
	push(@align_cmds, "perl $job_f");
}

foreach my $c (0..$#align_cmds)
{
	my $pid = $pm -> start($align_cmds[$c]) and next;
	my $job_num = $c+1;
	print STDERR "# \t\t---# Lastz ($job_num/$total_jobs)\n";
	system("$align_cmds[$c]");
	sleep(1);
	$pm -> finish($c);
}
print STDERR "# \t---# Wait for all job completion...\n";
$pm -> wait_all_children;

print STDERR "# \t---# Chain merge sort...\n";
`chainMergeSort $chain_dir/*/*.chain > $out_dir/all.chain`;

print STDERR "# \t---# Chain net...\n";
`chainNet $out_dir/all.chain -minSpace=1 $ref_size $tar_size stdout /dev/null | netSyntenic stdin $out_dir/all.net`;

my $chainNet_dir_1 = "";
if($chainNet_dir ne ""){
print STDERR "8. Split Chain net...\n";
#		$chainNet_dir_1 = "$/$bref_name/$btar_name";
		`splitChainNet.pl --chain $out_dir/all.chain --net $out_dir/all.net -o $chainNet_dir > $out_dir/log.chainnetsplit.txt 2>&1`;
}

print STDERR "All done.\n";
if(-d $chainNet_dir_1){
    `mkdir -p /raca_wd/reference/list/$bref_name`;
#    `touch /raca_wd/reference/list/$bref_name/$btar_name.done`;
	`rm -rf $out_dir`;
}
$time = `date`;
chomp($time);
print STDERR "# ---# [ END ] $time\n";
print STDERR "# ---# whole-genome alignment is finished (chainNet output: $chainNet_dir/$bref_name/$btar_name)\n";
