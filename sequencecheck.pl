#!/usr/bin/perl -w

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# MANUAL FOR creating a seperate fasta and quality score files from fastq file.pl

=pod

=head1 NAME

$0 -- to convert the fastq to fasta file and keep the quality scores in a seperate file as well

=head1 SYNOPSIS

convertfastqtofasta_qualityscoresseperate.pl --in1 file.fastq[.gz] [--zip] [--help] [--manual]

=head1 DESCRIPTION

Accepts one fastq files. Convert to fasta and quality score file.
The sequences are extracted to an output fasta and quality score file.
This handles only paired end & single reads properly.

=head2 FASTQ FORMAT

Default: Sequence name is all characters between beginning '@' and first space or '/'.  Also first 
character after space or '/' must be a 1 or 2 to indicate pair relationship.  
Compatible with both original Illumina (Casava version 1.7 and previous) and newer
(Casava >= 1.8) formatted headers:
  @B<HWUSI-EAS100R:6:73:941:1973#0>/I<1>
  @B<EAS139:136:FC706VJ:2:2104:15343:197393> I<2>:Y:18:ATCACG

=head1 OPTIONS

=over 3

=item B<-1, -a, -in, --in1>=FILE

One input file must be specified (e.g. parent paired reads).  (Required) 

=item B<-z, --zip>

Specify that input is gzipped fastq (output files will also be gzipped).  Slower, but space-saving.  (Optional)

=item B<-h, --help>

Displays the usage message.  (Optional) 

=item B<-m, --manual>

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

Copyright 2012 MOA.  
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.  
This is free software: you are free to change and redistribute it.  
There is NO WARRANTY, to the extent permitted by law.  

Please acknowledge author and affiliation in published work arising from this script's 
usage <http://bioinformatics.udel.edu/Core/Acknowledge>.

=cut

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - U S E R  V A R I A B L E S- - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# CODE FOR convertfastqtofasta_qualityscoresseperate.p

use strict;
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Getopt::Long;
use File::Basename;
use Pod::Usage;

#ARGUMENTS
my($in1,$zip,$help,$manual);

GetOptions (
                                "1|a|in|in1=s"       =>      \$in1,
                                "z|zip"         =>      \$zip,
                                "h|help"        =>      \$help,
                                "m|manual"      =>      \$manual);

# VALIDATE ARGS
pod2usage(-verbose  => 2)  if ($manual);
pod2usage( -verbose => 1 )  if ($help);
pod2usage( -msg  => "ERROR!  Required argument -1 is not found.\n", -exitval => 2, -verbose => 1)  if (! $in1 );

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - G L O B A L  V A R I A B L E S- - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# PARSE OUTPUT FILEBASE
my $out = fileparse($in1, qr/\.[^.]*(\.gz)?$/);

# FILE HANDLES
my ($DATA,$OUT1, $OUT2);

# OPEN INPUT FILE(s)(in1 &/ in2)
if($zip) {
   $DATA = IO::Uncompress::Gunzip->new( $in1 ) or die "IO::Uncompress::Gunzip failed: $GunzipError\n";
}
else {
   open ($DATA,$in1) || die $!;
}


my $sequence = 0;
my $quality = 0;
my $spec = 0;
my ($firstsyntax, $first, $second, $secondsyntax);
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - M A I N - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
##

while (<$DATA>) {
chomp;
   if (/^@(\S+.*)$/ && $sequence == 0 && $quality == 0){
      $sequence = 1;
   }
   elsif ($sequence == 1 && $quality == 0) {
   	  $firstsyntax = $_;
      $first = length($_);
      $sequence = 2;
   }
   elsif (/^\+/ && $sequence == 2 && $quality == 0) {
      $quality = 1;
   }
   elsif ($quality == 1 && $sequence == 2) {
      $secondsyntax = $_;
      $second = length($_);
      $quality = 0;
      $sequence = 0;
      $spec = 1;
   }
   if ($sequence == 0 && $quality == 0 && $spec == 1){
   		$spec = 0;
   		unless ($first==$second){
   		 print "$firstsyntax\t$secondsyntax\n";
   		}
   }
}
close $DATA;
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - S U B R O U T I N E S - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -

