#!/usr/bin/perl -w

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# MANUAL FOR SNPmergeranalysisPipeline.pl

=pod

=head1 NAME


=head1 SYNOPSIS

SNPmergeranalysisPipeline.pl --in1 fastqfolder/ --in2 outputfolder/ [--help] [--manual]

=head1 DESCRIPTION


=head1 COPYRIGHT

Copyright 2016 MOA.  

Please acknowledge author and affiliation in published work arising from this script's 
usage.
=cut

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - U S E R  P O D  V A R I A B L E S - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

use strict;
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Getopt::Long;
use File::Basename;
use Pod::Usage;

#ARGUMENTS
my($in1,$in2,$help,$manual);

GetOptions (	
				"1|in|list=s"		=>	\$in1,
				"2|out|output=s"		=>	\$in2,
                                "h|help"       		=>      \$help,
                                "man|manual"		=>      \$manual);

# VALIDATE ARGS
pod2usage( -verbose => 2 )  if ($manual);
pod2usage( -verbose => 1 )  if ($help);
pod2usage( -msg  => "ERROR!  Required argument -1 (sam/bam files list) not found/made. \n", -exitval => 2, -verbose => 1)  if (! $in1 );

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - M A I N  W O R K F L O W - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#open status directory
my $status = '/home/modupe/status/list-'.`date +%Y_%m_%d_%T`; chomp $status; $status = $status.'.txt';
open (STATUS, ">$status");

#LOGFILES
my $std_out = '/home/modupe/.LOG/SNPmerger-'.`date +%m-%d-%y_%T`; chomp $std_out; $std_out = $std_out.'.log';
my $std_err = '/home/modupe/.LOG/SNPmerger-'.`date +%m-%d-%y_%T`; chomp $std_err; $std_err = $std_err.'.err';

open(STDOUT, '>', "$std_out") or die "Log file doesn't exist";
open(STDERR, '>', "$std_err") or die "Error file doesn't exist";

#OPEN LIST
open(IN, "<$in1") or die "$in1 doesn't exist";

#DESTINATION VARIABLES
my $PICARDDIR="/home/modupe/.software/picard-tools-1.119";
my $REF="/home/modupe/DATA/genome/TophatChickengenome.fa";
my $TMP_DIR="/home/modupe/DATA/temp";
my $gtf_file="/home/modupe/DATA/gtf/Galgal4_76.gtf";
my $index="/home/modupe/DATA/TophatIndex/Tophat";
my $GATKDIR="/home/modupe/.software/GenomeAnalysisTK-3.3-0";
my $snpEFF="/home/modupe/.software/snpEff";
my $filteringScript="/home/modupe/SCRIPTS/VARIANT_analysis/01-filteringDP5.pl";
my $VEP="/home/modupe/.software/ensembl-tools-release-78/scripts/variant_effect_predictor/variant_effect_predictor.pl";
my $snpEFFgff="/home/modupe/DATA/snpEFFgff";
my $Database="/home/modupe/SCRIPTS/INSERT-VARIANTS/specificDatabaseVariants.pl";

#WORKING....
my $picardsyntax="java -jar $PICARDDIR/MergeSamFiles.jar TMP_DIR=$TMP_DIR ASSUME_SORTED=true ";

while (<IN>){
	chomp;
	$picardsyntax .= "INPUT=$_ ";
}
if($in2){
	$picardsyntax .= "OUTPUT=$in2.bam";
} else {
	$picardsyntax .= "OUTPUT=/home/modupe/newmergedfile.bam";
}

my $time = `date`;
print STATUS "running MergeSAMfiles\t$time\n==\n";
`time $picardsyntax`;

if($in2){
	`samtools sort $in2.bam $in2.sorted`;
	`samtools index $in2.sorted.bam`; 
}else {
	`samtools sort /home/modupe/newmergedfile.bam /home/modupe/newmergedfile.sorted`;
	`samtools index /home/modupe/newmergedfile.sorted.bam`;
}
close(IN);
close (STATUS);
