#!usr/bin/perl
use File::Basename;

# multiple files for BWA & samtools or picardtools

#VARIABLES
my ($filename, $bwa, $picard, $storedfilename );
my $number = 0;

#FILE VARIABLES to change
my $fastqfolder="/home/amodupe/CARL/FINAL";
my $folder="/home/amodupe/CARL/FINAL";

#-------------

#DESTINATION VARIABLES
my $PICARDDIR="/usr/local/picard-tools-1.67";
my $REF="/home/amodupe/CARL/Galgal4_genome.fa";
my $TMP_DIR="/home/amodupe/CARL/temp";

#--------------

###WORKING
#OPENING DIRECTORY
opendir(DIR,$fastqfolder) or die "Folder \"$fastqfolder\" doesn't exist\n";
my @Directory = readdir(DIR);
close(DIR);

#make bwa_file directory
system "mkdir $folder/bwa_files";

foreach my $FILE (@Directory){
	if ($FILE =~ /.*\.gz$/){
		$number++;
		$filename = fileparse($FILE, qr/(\.f.*\.gz)?$/);
		##BWA
		$bwa="bwa mem -t 8 -M -R '\@RG\\tID:Label\\tLB:Label\\tSM:Label' $REF $fastqfolder/$FILE > $folder/bwa_files/$filename.sam";
		##PICARD
		$picard="java -jar $PICARDDIR/SortSam.jar TMP_DIR=$TMP_DIR INPUT=$folder/bwa_files/$filename.sam OUTPUT=$folder/bwa_files/$filename"."_sorted.bam SO=coordinate";

#PERMANENT
my $shell = "#!/bin/sh
#PBS -N BWA-$number
#PBS -S /bin/bash
#PBS -V
#PBS -l walltime=300:00:00,cput=300:00:00,nodes=1
#PBS -q long
#PBS -d ";

		# Print to shell script
		open(OUT,">$number-bwa.sh");
		print OUT $shell;	
		print OUT "$folder\n";
		print OUT "time $bwa\n";
		print OUT "time $picard\n";
		close(OUT);

		print `qsub $number-bwa.sh`;
		system "rm -rf $number-bwa.sh";  ##removing the shellscript
		$storedfilename .= "$folder/bwa_files/$filename"."_sorted.bamrrr";
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


