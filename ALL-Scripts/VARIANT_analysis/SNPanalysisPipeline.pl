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
It also imports the results into the database : CHICKENSNPS
=head2 FASTQ FORMAT

Default: Sequence name is all characters between beginning '@' and first space or '/'.  Also first 
character after space or '/' must be a 1 or 2 to indicate pair relationship.  
Compatible with both original Illumina (Casava version 1.7 and previous) and newer
(Casava >= 1.8) formatted headers:
  @B<HWUSI-EAS100R:6:73:941:1973#0>/I<1>
 
=head2 OUTPUT FORMAT

OUTPUTS are a vcf file, filtered_DP5 vcf files, SNPeff_GTF, SNPef-GFF and VEP annotated files.

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
my($fastqfolder,$destfolder,$help,$manual);

GetOptions (	
				"1|in1|folder=s"	=>	\$fastqfolder,
				"2|in2|output=s"	=>	\$destfolder,
                                "h|help"       		=>      \$help,
                                "man|manual"		=>      \$manual);

# VALIDATE ARGS
pod2usage( -verbose => 2 )  if ($manual);
pod2usage( -verbose => 1 )  if ($help);
pod2usage( -msg  => "ERROR!  Required argument -1 (fastq folder) not found. \n", -exitval => 2, -verbose => 1)  if (! $fastqfolder );
pod2usage( -msg  => "ERROR!  Required argument -2 (output folder) not found/made. \n", -exitval => 2, -verbose => 1)  if (! $destfolder );

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
my $Mfolder="$destfolder/MAPPING";
my $Vfolder="$destfolder/VARIANTS";

my $filelocation = `ls $Mfolder`; # making sure the FASTQ file isn't reprocessed
unless (length $filelocation > 1){
	system "mkdir $Mfolder"; system "mkdir $Vfolder";
	$time = `date`;
	print STATUS "created $Mfolder and $Vfolder at \t$time\n===\n";
}
else {
	$time = `date`;
	print STATUS " Folders $Mfolder and $Vfolder already exist at \t$time\n===\n";
}
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
my $Database="/home/modupe/SCRIPTS/INSERT-VARIANTS/specificDatabaseVariants.pl";
#--------------

###WORKING
#OPENING DIRECTORY
opendir(DIR,$fastqfolder) or die "Folder \"$fastqfolder\" doesn't exist\n";
my @Directory = readdir(DIR);
close(DIR);
my $check = 0;

#processing each fastq file
foreach my $FILE (@Directory){
	if ($FILE =~ /^(s*(\d*)).*\.f.*q\.gz$/){ # parsing FASTQ name
		$check = 1;
		$libraryno = $2; #library_id
		my $raw_reads = `find $fastqfolder -name '$1*gz'`; 
		my @allraw_reads = split("\n", $raw_reads);
		$DATA = '';
		if ($#allraw_reads > 0){
			foreach my $i (0..$#allraw_reads){
				$DATA .= "$allraw_reads[$i] ";
			}
		}
		else {
			$DATA = $raw_reads;
		}
    }
    if ($check == 1){
		$libname = "library-$libraryno";
		$number++;
		
		$Mfile = "$Mfolder/$libname";
		$Vfile = "$Vfolder/$libname";
		
        $filename = `ls $Mfile`; # making sure the FASTQ file isn't reprocessed
        my @locations = split("\n", $filename);
        unless (length $filename < 1){
        	next;
        }
        		
		system "mkdir $Mfile"; system "mkdir $Vfile";
		#make tophat directory
		$tophatfolder = "$Mfile/tophat_out";
		system "mkdir $tophatfolder";

		##TOPHAT
		$tophat="tophat -p 24 --rg-id Label --rg-sample Label --rg-library Label --no-coverage-search -o $tophatfolder $index $DATA";

		##PICARD
		$picard="java -jar $PICARDDIR/SortSam.jar TMP_DIR=$TMP_DIR INPUT=$tophatfolder/accepted_hits.bam OUTPUT=$tophatfolder/$libname.bam SO=coordinate";

		#MARKDUPLICATES
		$markduplicates = "java -jar $PICARDDIR/MarkDuplicates.jar INPUT=$tophatfolder/$libname.bam OUTPUT=$Mfile/".$libname."_mdup.bam M=$Mfile/".$libname."_mdup.metrics CREATE_INDEX=true";

		#SPLIT&TRIM
		$splittrim="java -jar $GATKDIR/GenomeAnalysisTK.jar -T SplitNCigarReads -R $REF -I $Mfile/".$libname."_mdup.bam -o $Mfile/".$libname."_split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS";

		#GATK
		$gatksnps = "java -jar $GATKDIR/GenomeAnalysisTK.jar -T HaplotypeCaller -U ALLOW_N_CIGAR_READS -R $REF -I $Mfile/".$libname."_split.bam -o $Vfile/$libname.vcf";

		$time = `date`;
		print STATUS "\nWorking on library $libraryno\t$time\n==\n";
		
		$time = `date`;
		print STATUS "running TOPHAT\t$time\n==\n";
		`time $tophat`;

		$time = `date`;
		print STATUS "running PICARD\t$time\n==\n";
		`time $picard`;
		
		$time = `date`;
		print STATUS "running MarkDuplicates\t$time\n==\n";
		`time $markduplicates`;

		$time = `date`;
		print STATUS "running Split and Trim\t$time\n==\n";
		$time = `date`;
		`time $splittrim`;

		$time = `date`;
		print STATUS "running Haplotype Caller\t$time\n==\n";
		`time $gatksnps`;
		
		#perl to select DP > 5 & get header information
		$time = `date`;
		print STATUS "running filtering DP>5\t$time\n==\n";
		`time perl $filteringScript $Vfile  $Vfile/$libname.vcf`;

		##ANNOTATION types 
		my ($vepsyntax, $snpEFFGTF, $snpeffGFF);

		#running snpEFF-GTF
		$time = `date`;
		print STATUS "running snpEFF-GTF\t$time\n==\n";
		$snpEFFGTF ="java -Xmx4g -jar $snpEFF/snpEff.jar eff -c $snpEFF/snpEff.config -v -i vcf -o gatk -s $Vfile/".$libname."_GTF.html Galgal4.76 $Vfile/".$libname."_DP5.vcf > $Vfile/".$libname."_GTF.vcf";
		`time $snpEFFGTF`; 
		
		#running snpEFF-GFF
		$time = `date`;
		print STATUS "running snpEFF-GFF\t$time\n==\n";
		`cp -rf $snpEFFgff/* ./`;
		$snpeffGFF="java -jar $snpEFF/snpEff.jar eff -v -i vcf -o gatk -s $Vfile/".$libname."_GFF.html Galgal4_76 $Vfile/".$libname."_DP5.vcf > $Vfile/".$libname."_GFF.vcf"; 
		`time $snpeffGFF`;
		
		print STATUS "removing snpeffGFF folder\t$time\n==\n";
		opendir(GFFDIR,$snpEFFgff) or die "Folder \"$snpEFFgff\" doesn't exist\n";
		my @SnpGFFdir = readdir(GFFDIR);
		foreach (@SnpGFFdir){`rm -rf ./$_`;}
		close(GFFDIR);

		#running VEP
		$time = `date`;
		print STATUS "running VEP\t$time\n==\n";
		$vepsyntax="perl $VEP -i $Vfile/".$libname."_DP5.vcf --fork 24 --species Chicken  --cache --everything on --output_file $Vfile/".$libname."_VEP.txt";
		`time $vepsyntax`;

		$check = 0;
		$time = `date`;
		print STATUS "running into Database : CHICKENSNPs\t$time\n==\n";
		`time perl $Database -1 $Vfile -2 $libname`;
		
		print STATUS "Successful for library $libraryno\n==\n";
	}
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - -T H E  E N D - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

$time = `date`;
print STATUS "================DONE!!!!================\n\t$time\n";
close (STATUS);
exit;


