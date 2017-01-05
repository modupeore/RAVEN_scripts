#!/usr/bin/perl -w

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# MANUAL FOR fastq-random-sample.pl

=pod

=head1 NAME

$0 -- randomly sampled fastq sequenced data based on defined percentage or number

=head1 SYNOPSIS

fastq-random-sample.pl --in1 fwd.fastq[.gz] [--in2 rev.fastq[.gz]] --number 
		[--zip] [--percentage]	[--batches] [--help] [--manual]

=head1 DESCRIPTION

Accepts one/two fastq files. Based on user selection either a defined percentage or number.
The sequences will be randomly sampled into the output fastq.'
This handles single end and paired end reads properly.

=head2 FASTQ FORMAT

Default: Sequence name is all characters between beginning '@' and first space or '/'.  Also first 
character after space or '/' must be a 1 or 2 to indicate pair relationship.  
Compatible with both original Illumina (Casava version 1.7 and previous) and newer
(Casava >= 1.8) formatted headers:
  @B<HWUSI-EAS100R:6:73:941:1973#0>/I<1>
  @B<EAS139:136:FC706VJ:2:2104:15343:197393> I<2>:Y:18:ATCACG
 
=head1 OPTIONS

=over 3

=item B<-1, --i, --in, --in1>=FILE

One input file must be specified (e.g. forward reads).  (Required) 

=item B<-2, --in2>=FILE

Paired fastq file of the first input file (e.g. reverse reads).  (Optional) 

=item B<-n, --number>

Specify the number/percentage of sequences that will be randomly sampled. (Required)

=item B<-z, --zip>

Specify that input is gzipped fastq (output files will also be gzipped).  Slower, but space-saving.  (Optional)

=item B<-p, --percentage>

Specify that the number specified is in percentage. (Optional)

=item B<-b, --batches>

Specify the number of batches of the sampled reads you want. (Optional)

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

# CODE FOR fastq-random-sample.pl

use strict;
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Getopt::Long;
use File::Basename;
use Pod::Usage;

#ARGUMENTS
my($in1,$in2,$zip,$number,$percentage,$help,$manual,$decision);
my $index = 1;

GetOptions (	
				"1|i|in|in1=s"	=>	\$in1,
				"2|in2=s"	=>	\$in2,
				"z|zip" 	=>	\$zip,
				"n|number=i" 	=>	\$number,
				"p|percentage" 	=>	\$percentage,
				"h|help"	=>	\$help,
				"m|manual"	=>	\$manual,
				"b|batches=i"	=>	\$decision);

# VALIDATE ARGS
pod2usage(-verbose  => 2)  if ($manual);
pod2usage( -verbose => 1 )  if ($help);
pod2usage( -msg  => "ERROR!  Required arguments -1 and/or -2, and -n are not found. \n", -exitval => 2, -verbose => 1)  if (! $in1 || ! $number );
if (! $decision){$decision=1};

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - G L O B A L  V A R I A B L E S- - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

while ($index <=$decision){

# PARSE OUTPUT FILEBASE
my ($out1, $out2);
$out1=fileparse($in1, qr/\.[^.]*(\.gz)?$/);
if ($in2){
   $out2=fileparse($in2, qr/\.[^.]*(\.gz)?$/);
}
# FILE HANDLES
my ($DATA1,$DATA2,$OUT1,$OUT2);

# OPEN INPUT FILE(s)(in1 &/ in2)
if($zip) {
   $DATA1 = IO::Uncompress::Gunzip->new( $in1 ) or die "IO::Uncompress::Gunzip failed: $GunzipError\n";
   if($in2) {
      $DATA2 = IO::Uncompress::Gunzip->new( $in2 ) or die "IO::Uncompress::Gunzip failed: $GunzipError\n";
   }
}
else {
   open ($DATA1,$in1) || die $!;
   if($in2) {
   open ($DATA2,$in2) || die $!;
   }
}

#OPEN OUTPUT FILE(s)
if($zip){
   $OUT1 = IO::Compress::Gzip->new( "$out1-$index.selected.fastq.gz" ) or die "IO::Compress::Gzip failed: $GzipError\n";
   if($in2){
      $OUT2 = IO::Compress::Gzip->new( "$out2-$index.selected.fastq.gz" ) or die "IO::Compress::Gzip failed: $GzipError\n";
   }
}
else {
   open ($OUT1,"> $out1-$index.selected.fastq") or die $!;
   if($in2){
      open ($OUT2,"> $out2-$index.selected.fastq") or die $!;
   }
}

# REG EXP FOR FASTQ HEADERS
my $hdr_ptrn;
# Illumina 1-1.7:	@HWUSI-EAS100R:6:73:941:1973#0/1
# Illumina 1.8+:	@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG
$hdr_ptrn='^\@(\S+)[\/ ]*[12]*';

#HASH TABLES & VARIABLES
my (%INhash, %IN1hashdetails, %IN2hashdetails); my $uid=0;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - M A I N - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
##
## PROCESS in1 for single reads
#
print "Reading forward strand";
while(<$DATA1>) {
  if (/$hdr_ptrn/){
      $INhash{$uid} = $1;
      $IN1hashdetails{$1}=$_;
      $IN1hashdetails{$1}.=<$DATA1>.<$DATA1>.<$DATA1>;
      $uid++;
   }
   else {
      die "\nERROR! File format error in $in1 near line ".$..".\n$_\n";
   }
}
close $DATA1;
print " Done\n";
##
## PROCESS in2 for paired reads (optional)
#
if ($in2){
print "Reading reverse strand";
   undef %INhash;
   $uid=0;
   while(<$DATA2>){
      if (/$hdr_ptrn/){
         if ($IN1hashdetails{$1}){
            $INhash{$uid} = $1;
            $_.=<$DATA2>.<$DATA2>.<$DATA2>;
	    $IN2hashdetails{$1}= $_;
            $uid++;
         }
      }
   }
   close $DATA2;
print "  Done\n";
   
   #OUTPUT PAIRED READS
   #PERCENTAGE/NUMBER
print "Processing output for paired reads";
   if ($percentage){$number = sprintf "%.0f", ($number/100)*$uid;} else {$number = $number;}
   if ($number>$uid){$number = $uid;}

   foreach my $no (1..$number){
      my $check = 1;
      my $s_header = '';
      my $s_uid = '';
      while ($check == 1){
         $s_uid = int(rand($uid));
         if (exists $INhash{$s_uid}){
            $s_header = $INhash{$s_uid};
            delete $INhash{$s_uid};
            $check = 0;
         }
      }
   print $OUT1 "$IN1hashdetails{$s_header}";
   print $OUT2 "$IN2hashdetails{$s_header}";
   }
   close $OUT1;
   close $OUT2;
print "  Done\n";
} #end if $in2
else{
   #OUTPUT FOR SINGLE READS
   #PERCENTAGE/NUMBER
print "Processing output for single reads";
   if ($percentage){$number = sprintf "%.0f", ($number/100)*$uid;} else {$number = $number;}
   if ($number>$uid){$number = $uid;}

   foreach my $no (1..$number){
      my $check = 1;
      my $s_header = '';
      my $s_uid = '';
      while ($check == 1){
         $s_uid = int(rand($uid));
         if (exists $INhash{$s_uid}){
            $s_header = $INhash{$s_uid};
            delete $INhash{$s_uid};
            $check = 0;
         }
      }
   print $OUT1 "$IN1hashdetails{$s_header}";
   }
   close $OUT1;
}
$index++;
print "  Done\n";
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
print "\n\n*********DONE*********\n\n";
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - S U B R O U T I N E S - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -

