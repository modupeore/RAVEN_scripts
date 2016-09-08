#!/usr/bin/perl -w

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# MANUAL FOR SNPmergeranalysisPipeline.pl

=pod

=head1 NAME

$0 -- Variant detection and functional annotation pipeline for Tophat bam file

=head1 SYNOPSIS

SNPmergeranalysisPipeline.pl --in1 fastqfolder/ --in2 outputfolder/ [--help] [--manual]


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
my($in1,$folder,$help,$manual);

GetOptions (	
				"1|in|list=s"		=>	\$in1,
				"2|out|output=s"		=>	\$folder,
                                "h|help"       		=>      \$help,
                                "man|manual"		=>      \$manual);

# VALIDATE ARGS
pod2usage( -verbose => 2 )  if ($manual);
pod2usage( -verbose => 1 )  if ($help);
pod2usage( -msg  => "ERROR!  Required argument -1 (sam/bam files list) not found/made. \n", -exitval => 2, -verbose => 1)  if (! $in1 );
pod2usage( -msg  => "ERROR!  Required argument -2 (folder) not found. \n", -exitval => 2, -verbose => 1)  if (! $folder );

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
open(IN, "<$in1") or die "$in1 doesn't exist";
`mkdir $folder`;
#sort samfiles
my $number = 0;
while (<IN>){
	chomp;
	$number++;
	my $samfile = $folder."/S".$number."sorted";
	`samtools sort $_ $samfile`;
}
close (IN);
print "\n$number\n";

#index samfiles
opendir(FOLDER,$folder) or die "Folder \"$folder\" doesn't exist\n";
my @Directory = readdir(FOLDER);
close(FOLDER);
foreach my $sam (@Directory){
	if ($sam =~ /^S.*bam$/){
		my $newsam = $folder."/".$sam;
		print "samtools index $newsam\n";
		`samtools index $newsam`;
	}
}

#merge samfiles
my $finalsam = $folder."/finalmerged.bam";
my $mergesyntax = "samtools merge $finalsam ";
foreach my $merge (@Directory){
	if ($merge =~ /^S.*bam$/){
		$mergesyntax .= "$folder/$merge ";
	}
}
`time $mergesyntax`;

close (STATUS);
