#!/usr/bin/perl -w

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# MANUAL FOR sorting a fasta file to be GATK compartible.

=pod

=head1 NAME

$0 -- to convert the reference fasta file to be GATK compartible with the .fai index.
=head1 SYNOPSIS

GATKSORT.pl --fa file.fasta[.gz] [--help] [--manual]

=head1 DESCRIPTION

Accepts one reference fasta file. Convert to the GATK comparible fasta file with its .fai index file
The output are sorted in the respective fasta (.fna) and index (.fna.fai) index files.

=head1 OPTIONS

=over 3

=item B<-1, -fa, --fasta, --in1>=FILE

One input file must be specified (e.g. fasta.fa ).  (Required) 

=item B<-h, --help>

Displays the usage message.  (Optional) 

=item B<-m, --manual>

Displays full manual.  (Optional) 

=back

=head1 DEPENDENCIES

Requires the following Perl libraries (all standard in most Perl installs).
   IO::Uncompress::Gunzip
   Getopt::Long
   File::Basename
   Pod::Usage

=head1 AUTHOR

Written by Modupe Adetunji, 
Bioinformatics Research Assistant, University of Delaware.

=head1 REPORTING BUGS

Report bugs to amodupe@udel.edu

=head1 COPYRIGHT

Copyright 2016 MOA.  
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.  
This is free software: you are free to change and redistribute it.  
There is NO WARRANTY, to the extent permitted by law.  

Please acknowledge author and affiliation in published work arising from this script's usage

=cut

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - U S E R  V A R I A B L E S- - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# CODE FOR confqtofaandqual.pl

use strict;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Getopt::Long;
use File::Basename;
use Pod::Usage;

#ARGUMENTS
my($in1,$zip,$help,$manual);

GetOptions (
                                "1|fa|fasta|in1=s"       =>      \$in1,
                                "h|help"        =>      \$help,
                                "m|manual"      =>      \$manual);

# VALIDATE ARGS
pod2usage(-verbose  => 2)  if ($manual);
pod2usage( -verbose => 1 )  if ($help);
pod2usage( -msg  => "ERROR!  Required argument -1 is not found.\n", -exitval => 2, -verbose => 1)  if (! $in1);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - G L O B A L  V A R I A B L E S- - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
$/ = "\>";

# PARSE OUTPUT FILEBASE
my $out = fileparse($in1, qr/\.[^.]*(\.gz)?$/);
if ($in1 =~ /\.gz$/){$zip=1;}
# FILE HANDLES
my ($DATA,$OUT1,$OUT2, %SEQ, %SEQnum, %SEQheader);

# OPEN INPUT FILE(s)(in1 &/ in2)
if($zip) {
   $DATA = IO::Uncompress::Gunzip->new( $in1 ) or die "IO::Uncompress::Gunzip failed: $GunzipError\n";
}
else {
   open ($DATA,$in1) or die $!;
}

#OPEN OUTPUT FILE(s)
open($OUT1, "> $out.fna") or die $!;
open($OUT2, "> $out.fna.fai") or die $!;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - M A I N - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
##
#  if (/^@(\S+[\/ ][12].*)$/){
my @fastqfile = <$DATA>;
shift(@fastqfile);
foreach my $entry (@fastqfile){
  my @pieces = split(/\n/, $entry);
  $pieces[0] = (split(' ',$pieces[0]))[0];
  my $seq = '';
  foreach my $num (1.. $#pieces-1){
    $seq .= $pieces[$num];
  }
  $seq .= substr($pieces[$#pieces],0,-1);
  $SEQ{$pieces[0]} = $seq;
  $SEQnum{$pieces[0]} = length($seq);
  $SEQheader{$pieces[0]} = length($pieces[0]);
}
my ($check, $start, $newstart, $last);
foreach my $header (sort keys %SEQ){
  if (length($header) >= 1) {
    print $OUT1 ">$header\n$SEQ{$header}\n";
    unless ($check){
      $start = $SEQheader{$header}+2;
      $last = $SEQnum{$header}+1;
      print $OUT2 "$header\t$SEQnum{$header}\t$start\t$SEQnum{$header}\t$last\n";
      $check = "yes";
    }
    else {
      $newstart = $SEQheader{$header}+2+$last+$start;
      $start = $newstart;
      $last = $SEQnum{$header}+1;
      print $OUT2 "$header\t$SEQnum{$header}\t$start\t$SEQnum{$header}\t$last\n";
      $check = "yes";
    }
  }
}
close $DATA; close $OUT1; close $OUT2;
$/ = "\n";
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -


