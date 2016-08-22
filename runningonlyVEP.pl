#!/usr/bin/perl -w

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# MANUAL FOR SNPanalysisPipeline.pl

=pod

=head1 NAME

$0 -- Variant detection and functional annotation pipeline for RNA-Seq fastq files

=head1 SYNOPSIS

SNPanalysisPipeline.pl --in1 fastqfolder/ --in2 outputfolder/ [--help] [--manual]

=head1 DESCRIPTION

Accepts a folder containing one or more fastq files.
This handles folders containing either zipped or not zipped files.
It runs through Tophat, PicardTools, GATK, Filtering steps,and gene annotation using SNPeff and VEP.

=head2 FASTQ FORMAT

Default: Sequence name is all characters between beginning '@' and first space or '/'.  Also first 
character after space or '/' must be a 1 or 2 to indicate pair relationship.  
Compatible with both original Illumina (Casava version 1.7 and previous) and newer
(Casava >= 1.8) formatted headers:
  @B<HWUSI-EAS100R:6:73:941:1973#0>/I<1>
 
=head2 OUTPUT FORMAT

OUTPUTS are a vcf file, filtered_DP5 vcf files, SNPeff_GTF, SNPeff_GTF-GFF, SNPef-GFF and VEP annotated files.

=head1 OPTIONS

=over 3

=item B<-1, --in1, --folder>=FOLDER

Specify Fastq files folder (e.g. the zipped files folder).  (Required) 

=item B<-2, --in2, --output>=FOLDER

Specify desired output file folder.  (Required)

=item B<-h, --help>

Displays the usage message.  (Optional) 

=item B<-man, --manual>

Displays full manual.  (Optional) 

=back

=head1 DEPENDENCIES

Requires the following Perl libraries (all standard in most Perl installs).
   IO::Compress::Gzip
   IO::Uncompress::Gunzip
   Getopt::Long
   File::Basename
   Pod::Usage

=head1 AUTHOR

Written by Modupe Adetunji, 
Schmidt's Lab, University of Delaware.

=head1 REPORTING BUGS

Report bugs to amodupe@udel.edu

=head1 COPYRIGHT

Copyright 2015 MOA.  

Please acknowledge author and affiliation in published work arising from this script's 
usage.
=cut

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - U S E R  P O D  V A R I A B L E S - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# CODE FOR SNPanalysisPipeline.pl

use strict;
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Getopt::Long;
use File::Basename;
use Pod::Usage;

#ARGUMENTS
my($folder,$help,$manual);

GetOptions (	
				"1|in1|folder=s"	=>	\$folder,
                "h|help"       		=>      \$help,
                "man|manual"		=>      \$manual);

# VALIDATE ARGS
pod2usage( -verbose => 2 )  if ($manual);
pod2usage( -verbose => 1 )  if ($help);
pod2usage( -msg  => "ERROR!  Required argument -1 (folder) not found. \n", -exitval => 2, -verbose => 1)  if (! $folder );

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - M A I N  W O R K F L O W - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#open status directory
my $status = '/home/modupe/status/'.`date +%Y_%m_%d_%T`; chomp $status; $status = $status.'.txt';
open (STATUS, ">$status");

# multiple files for tophat & picardtools on RAVEN

#VARIABLES
my ($libname, $Mfile, $Vfile, $libraryno, $filename);
my ($time, $picard, $tophat, $tophatfolder, $markduplicates, $splittrim, $gatksnps, $DATA);
my $number = 0;

#FILE VARIABLES to change

#-------------


#DESTINATION VARIABLES

my $VEP="/home/modupe/.software/ensembl-tools-release-78/scripts/variant_effect_predictor/variant_effect_predictor.pl";
my $snpEFFgff="/home/modupe/DATA/snpEFFgff";
#--------------

###WORKING
#OPENING DIRECTORY
opendir(DIR,$folder) or die "Folder \"$folder\" doesn't exist\n";
my @Directory = readdir(DIR);
close(DIR);

#processing each fastq file
foreach my $FILE (@Directory){
	if ($FILE =~ /^library-(\d*)$/){ # parsing FASTQ name 
	print "$FILE\n";
		$libraryno = $1; #library_id
		$libname = "library-$libraryno";
		$number++;
		
		$Vfile = "$folder/$libname";
		
		#running VEP
		$time = `date`;
		print STATUS "running VEP\t$time\n==\n";
		my $vepsyntax="perl $VEP -i $Vfile/".$libname."_DP5.vcf --fork 24 --species Chicken  --cache --everything on --output_file $Vfile/".$libname."_VEP.txt";
		`time $vepsyntax`;
	}
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - -T H E  E N D - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

$time = `date`;
print STATUS "================DONE!!!!================\n\t$time\n";
close (STATUS);
exit;


