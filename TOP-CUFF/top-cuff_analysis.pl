#!/usr/bin/perl -w

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# MANUAL FOR top-cuff_analysis.pl

=pod

=head1 NAME

$0 -- Tophat and cufflinks automated pipeline for an entire folder containing the zipped fastq files

=head1 SYNOPSIS

top-cuff_analysis.pl --in fastq.gz --num 24 --gtf genes.gtf --index genome.index 
		[--in2 OUTPUT_FOLDER] [--help] [--manual]

=head1 DESCRIPTION

Accepts a folder containing the zipped fastq files.
This handles only a folder containing zipped files.

=head2 FASTQ FORMAT

Default: Sequence name is all characters between beginning '@' and first space or '/'.  Also first 
character after space or '/' must be a 1 or 2 to indicate pair relationship.  
Compatible with both original Illumina (Casava version 1.7 and previous) and newer
(Casava >= 1.8) formatted headers:
  @B<HWUSI-EAS100R:6:73:941:1973#0>/I<1>
 
=head1 OPTIONS

=over 3

=item B<-1, --in, --folder>=FOLDER

Fastq file folder (e.g. the zipped files folder).  (Required) 

=item B<-2, --num>=NUMBER

Specify number of processors that will be used (e.g 24). (Required)

=item B<-3, --gtf>=FILE

Specify the GTF file (e.g genes.gtf). (Required)

=item B<-4, --index>=FILE

Specify the Reference genome index file. (Required)

=item B<-5, --out, --output>=FOLDER

Specify the output directory (e.g output folder. (Optional))

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
Center for Bioinformatics and Computational Biology Core Facility, University of Delaware.

=head1 REPORTING BUGS

Report bugs to amodupe@udel.edu

=head1 COPYRIGHT

Copyright 2013 MOA.  
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.  
This is free software: you are free to change and redistribute it.  
There is NO WARRANTY, to the extent permitted by law.  

Please acknowledge author and affiliation in published work arising from this script's 
usage <http://bioinformatics.udel.edu/Core/Acknowledge>.

=cut

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - U S E R  V A R I A B L E S- - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# CODE FOR top-cuff_analysis.pl

use strict;
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Getopt::Long;
use File::Basename;
use Pod::Usage;

#ARGUMENTS
my($in1,$in2,$num_cores,$gtf_file,$index,$help,$manual);

GetOptions (	
				"1|in|folder=s"	=>	\$in1,
				"2|num=i"	=>	\$num_cores,
				"3|gtf=s"	=>	\$gtf_file,
				"4|index=s"	=>	\$index,
				"5|out|output=s"=>	\$in2,
                                "h|help"        =>      \$help,
                                "man|manual"	=>      \$manual);

# VALIDATE ARGS
pod2usage( -verbose => 2 )  if ($manual);
pod2usage( -verbose => 1 )  if ($help);
pod2usage( -msg  => "ERROR!  Required argument -1 (input folder) not found. \n", -exitval => 2, -verbose => 1)  if (! $in1 );
pod2usage( -msg  => "ERROR!  Required argument -2 (number of processors) not found. \n", -exitval => 2, -verbose => 1)  if (! $num_cores );
pod2usage( -msg  => "ERROR!  Required argument -3 (GTF file) not found. \n", -exitval => 2, -verbose => 1)  if (! $gtf_file );
pod2usage( -msg  => "ERROR!  Required argument -4 (reference genome index file) not found. \n", -exitval => 2, -verbose => 1)  if (! $index );
if (!$in2) {$in2 = "OUTPUT_FOLDER"};

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - G L O B A L  V A R I A B L E S- - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#OPENING FOLDER
opendir(DIR,$in1);
my @Folder = readdir(DIR);
close(DIR);

#MAKE GENERAL OUTPUT FOLDER
system "mkdir $in2";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - M A I N - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ACCESSING EACH FILE IN FOLDER
foreach my $filename (@Folder){
	if ($filename =~ /\.gz$/){
	        my $filepath = "$in1/$filename";
        	#UNZIP FILE
		my $DATA=fileparse($filepath, qr/\.[^.]*(\.gz)?$/);
		gunzip $filepath => $DATA;

        	#OUTPUTFILENAME
		$filename=~ m/^s*(\d*).*.gz/;
		my $out_file = "library_$1"; 
		#PRINT STATUS.
		print "\n\n\t\t=== WORKING ON \"$out_file\" ===\n\n";

		#PARSE OUTPUT FILEBASE
        	my $out_dir = "$in2/$out_file";
        	system "mkdir $out_dir";
        	my $tophat = "$out_dir/tophat_out";
        	my $cufflinks = "$out_dir/cufflinks_out";
        	
		#TOPHAT
        	system "tophat -p $num_cores -G $gtf_file -o $tophat $index $DATA";

        	#CUFFLINKS
        	system "cufflinks -p $num_cores -g $gtf_file -o $cufflinks $tophat/accepted_hits.bam";
		
		#REMOVE UNZIP FILE
		system "rm -f $DATA";
	}
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
print "\n\n*********DONE*********\n\n";
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - S U B R O U T I N E S - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
