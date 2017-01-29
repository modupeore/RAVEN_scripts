#!/usr/bin/perl -w

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# MANUAL FOR SNPeveryanalysisPipeline.pl

=pod

=head1 NAME

$0 -- Variant detection and functional annotation pipeline for RNA-Seq fastq files

=head1 SYNOPSIS

SNPeveryanalysisPipeline.pl --in1 fastqfolder/ --in2 outputfolder/ [--help] [--manual]

=head1 DESCRIPTION

Accepts a folder containing one or more fastq files.
This handles folders containing either zipped or not zipped files.
It runs through Tophat, Cufflinks, PicardTools, GATK, Filtering steps,and gene annotation using SNPeff and VEP.

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
my($in1,$in2,$help,$manual);

GetOptions (	
				"1|in1|folder=s"	=>	\$in1,
				"2|in2|output=s"	=>	\$in2,
                                "h|help"       		=>      \$help,
                                "man|manual"		=>      \$manual);

# VALIDATE ARGS
pod2usage( -verbose => 2 )  if ($manual);
pod2usage( -verbose => 1 )  if ($help);
pod2usage( -msg  => "ERROR!  Required argument -1 (fastq folder) not found. \n", -exitval => 2, -verbose => 1)  if (! $in1 );
pod2usage( -msg  => "ERROR!  Required argument -2 (output folder) not found/made. \n", -exitval => 2, -verbose => 1)  if (! $in2 );

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - M A I N  W O R K F L O W - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#open status directory
my $status = '/home/modupe/status/'.`date +%Y_%m_%d_%T`; chomp $status; $status = $status.'.txt';
open (STATUS, ">$status");

# multiple files for tophat & picardtools on RAVEN

#VARIABLES
my ($filename, $time, $picard, $tophat, $cufflinks, $cufflinksfolder, $tophatfolder,$storedfilename, $DATA);
my $number = 0;

#FILE VARIABLES to change
my $fastqfolder=$in1;
my $destfolder=$in2;
my $folder="$destfolder/MAPPING";

system "mkdir $folder";
$time = `time`;
print STATUS "opened $folder\n$time\n===\n\n";
#-------------


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
my $snpDAT="/home/modupe/.software/SNPdat_package_v1.0.5";
#--------------

###WORKING
#OPENING DIRECTORY
opendir(DIR,$fastqfolder) or die "Folder \"$fastqfolder\" doesn't exist\n";
my @Directory = readdir(DIR);
close(DIR);

#processing each fastq file

foreach my $FILE (@Directory){
	if ($FILE =~ /.*\.gz$/){
		$number++;
		#parsing NAME of file
		$DATA = $fastqfolder."/".$FILE;
		$filename = fileparse($FILE, qr/\.[^.]*(\.gz)?$/);

$time = `date`;
print STATUS "\n\nWorking on library $filename\t$time\n==\n";

		system "mkdir $folder/$filename";
		#make tophat & cufflinks directory
		$tophatfolder = "$folder/$filename/tophat_out";
		$cufflinksfolder = "$folder/$filename/cufflinks_out";
		system "mkdir $tophatfolder";
		system "mkdir $cufflinksfolder";

		##TOPHAT
		$tophat="tophat -p 24 --rg-id Label --rg-sample Label --rg-library Label --no-coverage-search -o $tophatfolder $index $DATA";
                
        ##CUFFLINKS
        $cufflinks="cufflinks -p 24 -g $gtf_file -o $cufflinksfolder $tophatfolder/accepted_hits.bam";

		##PICARD
		$picard="java -jar $PICARDDIR/SortSam.jar TMP_DIR=$TMP_DIR INPUT=$tophatfolder/accepted_hits.bam OUTPUT=$tophatfolder/accepted_hits.sorted.bam SO=coordinate";


$time = `date`;
print STATUS "running TOPHAT\t$time\n==\n";
		`time $tophat`;
$time = `date`;
print STATUS "running CUFFLINKS\t$time\n==\n";
		`time $cufflinks`;
$time = `date`;
print STATUS "running PICARD\t$time\n==\n";
		`time $picard`;
		
		$storedfilename .= "$tophatfolder/accepted_hits.sorted.bamrrr";
	}
}

#merging output help.
my @storage = split("rrr", $storedfilename);
my $picardsyntax="java -jar $PICARDDIR/MergeSamFiles.jar TMP_DIR=$TMP_DIR ASSUME_SORTED=true ";
foreach my $figure (@storage){
	$picardsyntax .= "INPUT=$figure ";
	}
$picardsyntax .= "OUTPUT=$folder/newbwasorted.bam";

$time = `date`;
print STATUS "running MergeSAMfiles\t$time\n==\n";
`time $picardsyntax`;

#markduplicates
$time = `date`;
print STATUS "running MarkDuplicates\t$time\n==\n";
`time java -jar $PICARDDIR/MarkDuplicates.jar INPUT=$folder/newbwasorted.bam OUTPUT=$folder/newbwasorted_mdup.bam M=$folder/newbwassorted_mdup.metrics CREATE_INDEX=true`;

#splitandtrim
$time = `date`;
print STATUS "running Split and Trim\t$time\n==\n";
`time java -jar $GATKDIR/GenomeAnalysisTK.jar -T SplitNCigarReads -R $REF -I $folder/newbwasorted_mdup.bam -o $folder/newbwasplit.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS`;

##GATK
#split haplotype
$time = `date`;
print STATUS "running Haplotype Caller\t$time\n==\n";
`time java -jar $GATKDIR/GenomeAnalysisTK.jar -T HaplotypeCaller -U ALLOW_N_CIGAR_READS -R $REF -I $folder/newbwasplit.bam -o $destfolder/all_Unified.vcf`;

#perl to select DP > 5 & get header information
$time = `date`;
print STATUS "running Haplotype Caller\t$time\n==\n";
`time perl $filteringScript $destfolder  $destfolder/all_Unified.vcf`;


##ANNOTATION types 
my ($vepsyntax, $snpEFFsyntax, $snpdatsyntax);

#running snpEFF-GTF
$time = `date`;
print STATUS "running snpEFF-GTF\t$time\n==\n";
$snpEFFsyntax ="java -Xmx4g -jar $snpEFF/snpEff.jar eff -c $snpEFF/snpEff.config -v -i vcf -o gatk -s $destfolder/all_DP5_GTF.html Galgal4.76 $destfolder/all_Unified_DP5.vcf > $destfolder/all_DP5_GTF.vcf";
`time $snpEFFsyntax`; 
print STATUS "running snpEFF-GTF with GFF as the interval\t$time\n==\n";
`time java -Xmx4g -jar $snpEFF/snpEff.jar -c $snpEFF/snpEff.config -v -i vcf -o gatk -interval /home/modupe/DATA/gff/Galgal_NCBI.gff -s $destfolder/all_DP5_GTF-interval.html Galgal4.76 $destfolder/all_Unified_DP5.vcf > $destfolder/all_DP5_GTF-interval.vcf`;

#running snpEFF-GFF
$time = `date`;
print STATUS "running snpEFF-GFF\t$time\n==\n";
`cp -rf $snpEFFgff/* ./`;
`time java -jar $snpEFF/snpEff.jar eff -v -i vcf -o gatk -s $destfolder/all_DP5_GFF.html Galgal4_76 $destfolder/all_Unified_DP5.vcf > $destfolder/all_DP5_GFF.vcf`; 
print STATUS "removing snpeffGFF folder\t$time\n==\n";
opendir(GFFDIR,$snpEFFgff) or die "Folder \"$snpEFFgff\" doesn't exist\n";
my @SnpGFFdir = readdir(GFFDIR);
foreach (@SnpGFFdir){`rm -rf ./$_`;}
close(GFFDIR);

#running VEP
$time = `date`;
print STATUS "running VEP\t$time\n==\n";
$vepsyntax="perl $VEP -i $destfolder/all_Unified_DP5.vcf --fork 24 --species Chicken  --cache --everything on --output_file $destfolder/all_DP5_VEP.txt";
`time $vepsyntax`;



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - T H E  E N D - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

$time = `date`;
print STATUS "================DONE!!!!================\n\t$time\n";
close (STATUS);
exit;


