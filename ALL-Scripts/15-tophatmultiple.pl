#!usr/bin/perl
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use File::Basename;

# multiple files for tophat & picardtools

#VARIABLES
my ($filename, $bwa, $picard, $tophat, $cufflinks, $storedfilename, $tophatfolder, $DATA );
my $number = 0;

#FILE VARIABLES to change
my $fastqfolder="/home/amodupe/CARL/PERMANENT/sequences/ROSS_BREAST_RICK"; #PERMANENT/sequences/ROSS_BREAST_RICK";
my $folder="/home/amodupe/CARL/TOPHATBreastMuscle/ROSS_B_R";

#-------------

#DESTINATION VARIABLES
my $PICARDDIR="/usr/local/picard-tools-1.67";
my $REF="/home/amodupe/CARL/TophatChickengenome.fa";
my $TMP_DIR="/home/amodupe/CARL/temp";
my $gtf_file="/home/amodupe/CARL/Chicken_genome/Galgal.gtf";
my $index="/home/amodupe/CARL/Chicken/Tophat";
my $EXPORTPATH='export PATH=${PATH}:/usr/local/tophat-2.0.11:/usr/local/cufflinks-2.2.0';
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
		
		system "mkdir $folder/$filename";
		#make tophat & cufflinks directory
		$tophatfolder = "$folder/$filename/tophat_out";
		$cufflinksfolder = "$folder/$filename/cufflinks_out";
		system "mkdir $tophatfolder";
		system "mkdir $cufflinksfolder";

		##TOPHAT
		$tophat="tophat -p 24 --rg-id Label --rg-sample Label --rg-library Label -o $tophatfolder $index $DATA";
                ##CUFFLINKS
                $cufflinks="cufflinks -p 24 -g $gtf_file -o $cufflinksfolder $tophatfolder/accepted_hits.bam";

		##BWA
#		$bwa="bwa mem -t 8 -M -R '\@RG\\tID:Label\\tLB:Label\\tSM:Label' $REF $fastqfolder/$FILE > $folder/bwa_files/$filename.sam";

		##PICARD
		$picard="java -jar $PICARDDIR/SortSam.jar TMP_DIR=$TMP_DIR INPUT=$tophatfolder/accepted_hits.bam OUTPUT=$tophatfolder/accepted_hits.sorted.bam SO=coordinate";

#PERMANENT
my $shell = "#!/bin/sh
#PBS -N TUXEDO-$number
#PBS -S /bin/bash
#PBS -V
#PBS -l walltime=300:00:00,cput=300:00:00,nodes=1
#PBS -q long
#PBS -d ";

		# Print to shell script
		open(OUT,">$number-tuxedo.sh");
		print OUT $shell;	
		print OUT "$folder\n";
		print OUT "\n$EXPORTPATH\n\n";
		print OUT "\#TOPHAT\ntime $tophat\n\n";
		print OUT "\n\#CUFFLINK\ntime $cufflinks\n\n";
		print OUT "\n\#PICARD\ntime $picard\n\n";
		close(OUT);
		
		print `qsub $number-tuxedo.sh`;
		system "rm -rf $number-tuxedo.sh";  ##removing the shellscript
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

print "===SUBMITTED==\n";

open(syntaxtorun, ">$folder/finalsyntax.txt");
print syntaxtorun "$picardsyntax\n";
close (syntaxtorun);


