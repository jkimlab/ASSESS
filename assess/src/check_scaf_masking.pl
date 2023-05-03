#!/usr/bin/perl
use strict;
use warnings;

my $in_scaf = $ARGV[0];
my $in_cpu = $ARGV[1];
my $tar_name = $ARGV[2];
my $in_target_species = $ARGV[3];
my $out_dir = $ARGV[-1];

print STDERR "# ---# Masking start (repeatmasker)\n";
my $time = `date`;
chomp($time);
print STDERR "# \t---# [START] $time\n";

print STDERR "# \t---# reading sequence size\n";
my $masked_frac = `faSize -tab $in_scaf | grep "fracMasked" -`;
chomp($masked_frac);
my ($t, $f) = split(/\s+/,$masked_frac);
$f = $f+0;

print STDERR "# \t---# converting unrecognized characters\n";
`mkdir -p $out_dir/tmp`;
`sed 's/\\./_/g' $in_scaf > $out_dir/tmp/masked_assembly.rmsk1.fa`;
`sed 's/\\s/_/g' $out_dir/tmp/masked_assembly.rmsk1.fa > $out_dir/tmp/$tar_name.fa`;

if($f == 0){
		print STDERR "# \t---# sequence is not masked\n";
		print STDERR "# \t---# running repeatmasker\n";
		`RepeatMasker -dir $out_dir/tmp -s -a -pa $in_cpu -species $in_target_species $out_dir/tmp/$tar_name.fa`;
		`faToTwoBit $out_dir/tmp/$tar_name.fa $out_dir/tmp/$tar_name.2bit`;
		`twoBitMask $out_dir/tmp/$tar_name.2bit -add $out_dir/tmp/$tar_name.fa.out $out_dir/tmp/masked_assembly.rmsk.2bit`;
		`twoBitToFa $out_dir/tmp/masked_assembly.rmsk.2bit $out_dir/$tar_name.fa`;
		print "SEQFILE=$out_dir/$tar_name.fa\n";
}
else{
		print STDERR "# \t---# sequence is already masked (skipped)\n";
		`cp $in_scaf $out_dir/$tar_name.fa`;
		print "SEQFILE=$out_dir/$tar_name.fa\n";
}
$time = `date`;
chomp($time);
print STDERR "# \t---# [ END ] $time\n";
print STDERR "# ---# Masking finished (masking out: $out_dir/$tar_name.fa)\n";
