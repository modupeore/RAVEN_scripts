#!/usr/bin/perl -w

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# MANUAL FOR extractlibrarylist.pl

=pod

=head1 NAME

$0 -- Extracting libraries present in the mysql database : PENGUIN & CHICKENSNPS

=head1 SYNOPSIS

extractingallfromDB.pl [--help] [--manual]

 
=head1 OPTIONS

=over 3

=item B<-h, --help>

Displays the usage message.  (Optional) 

=item B<-man, --manual>

Displays full manual.  (Optional) 

=back

=head1 DEPENDENCIES

Requires the following Perl libraries (all standard in most Perl installs).
   DBI
   DBD::mysql
   Getopt::Long
   Pod::Usage

=head1 AUTHOR

Written by Modupe Adetunji, 
Animal and Food Sciences Department, University of Delaware.

=head1 REPORTING BUGS

Report bugs to amodupe@udel.edu

=head1 COPYRIGHT

Copyright 2015 MOA.  

=cut

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - U S E R  V A R I A B L E S- - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# CODE FOR extractingallfromDB.pl

use strict;
use DBI;
use DBD::mysql;
use Getopt::Long;
use Pod::Usage;

#ARGUMENTS
my($help,$manual);
GetOptions(
                                "h|help"        =>      \$help,
                                "man|manual"	=>      \$manual );

# VALIDATE ARGS
pod2usage( -verbose => 2 )  if ($manual);
pod2usage( -verbose => 1 )  if ($help);

# DATABASE ATTRIBUTES
my $dsn = 'dbi:mysql:PENGUIN';
my $dcn = 'dbi:mysql:CHICKENSNPS';
my $user = 'modupe';
my $passwd = 'penguin123';

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - G L O B A L  V A R I A B L E S- - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# DATABASE VARIABLES
my ($dbh, $dch, $sth, $syntax,@row,@row1,@row2,@row3,@row4,@row5,@row6,@row7);

#OPENING FILES
#PENGUIN
# open (TAB,">/home/modupe/Results.txt");
# open (ACT,">/home/modupe/Actual.txt");
# open (GEN,">/home/modupe/Genes.txt");
# open (ISO,">/home/modupe/Isoforms.txt");
# open (IDG,">/home/modupe/IdsGenes.txt");
#CHICKENSNPS
open (FIL,">/home/modupe/VFiles.txt");
open (RES,">/home/modupe/VResults.txt");
open (SUM,">/home/modupe/VSummary.txt");

#CONNECT TO THE DATABASE
print "\n\n\tCONNECTING TO THE DATABASE : $dsn and $dcn\n\n";
$dbh = DBI->connect($dsn, $user, $passwd) or die "Connection Error: $DBI::errstr\n";
$dch = DBI->connect($dcn, $user, $passwd) or die "Connection Error: $DBI::errstr\n";

# VARIABLES
my $number=0;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - M A I N  W O R K F L O W - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# SUMMARY OF RESULTS

# $syntax ="SELECT `COLUMN_NAME` FROM `INFORMATION_SCHEMA`.`COLUMNS` WHERE `TABLE_SCHEMA`=\'PENGUIN\' AND `TABLE_NAME`=\'RESULTS_SUMMARY\'";
# $sth = $dbh->prepare($syntax);
# $sth->execute or die "SQL Error: $DBI::errstr\n";
# while ($row = $sth->fetchrow_array() ) {
# 	print OUT "$row\t";
# }
# print OUT "\n";

#SELECTING RESULTS SUMMARY TABLE & printing out results
# $syntax = "select * from RESULTS_SUMMARY";
# $sth = $dbh->prepare($syntax);
# $sth->execute or die "SQL Error: $DBI::errstr\n";
# while (@row = $sth->fetchrow_array() ) {
# 	foreach my $list(0..$#row-1){print TAB "$row[$list]\t";}print TAB "$row[$#row]\n";
# }
# close TAB;

# $syntax = "select * from ACTUAL_FILES";
# $sth = $dbh->prepare($syntax);
# $sth->execute or die "SQL Error: $DBI::errstr\n";
# while (@row1 = $sth->fetchrow_array() ) {
# 	foreach my $list(0..$#row1-1){print ACT "$row1[$list]\t";}print ACT "$row1[$#row1]\n";
# }
# close ACT;

# $syntax = "select * from GENES_FPKM";
# $sth = $dbh->prepare($syntax);
# $sth->execute or die "SQL Error: $DBI::errstr\n";
# while (@row2 = $sth->fetchrow_array() ) {
# 	foreach my $list(0..$#row2-1){print GEN "$row2[$list]\t";}print GEN "$row2[$#row2]\n";
# }
# close GEN;

# $syntax = "select * from ISOFORMS_FPKM";
# $sth = $dbh->prepare($syntax);
# $sth->execute or die "SQL Error: $DBI::errstr\n";
# print "first";
# while (@row3 = $sth->fetchrow_array() ) {
# 	#print "yes";	
# 	foreach my $list(0..$#row3-1){print ISO "$row3[$list]\t";}print ISO "$row3[$#row3]\n";
# }
# print "\n\nsecond\n\n";
# close ISO;

# $syntax = "select * from GENE_ID";
# $sth = $dbh->prepare($syntax);
# $sth->execute or die "SQL Error: $DBI::errstr\n";
# while (@row4 = $sth->fetchrow_array() ) {
# 	foreach my $list(0..$#row4-1){print IDG "$row4[$list]\t";}print IDG "$row4[$#row4]\n";
# }
# close IDG;

$syntax = "select * from VARIANT_FILES";
$sth = $dch->prepare($syntax);
$sth->execute or die "SQL Error: $DBI::errstr\n";
while (@row5 = $sth->fetchrow_array() ) {
	foreach my $list(0..$#row5-1){print FIL "$row5[$list]\t";}print FIL "$row5[$#row5]\n";
}
close FIL;

$syntax = "select * from VARIANT_SUMMARY";
$sth = $dch->prepare($syntax);
$sth->execute or die "SQL Error: $DBI::errstr\n";
while (@row6 = $sth->fetchrow_array() ) {	
	foreach my $list(0..$#row6-1){print SUM "$row6[$list]\t";}print SUM "$row6[$#row6]\n";
}
close SUM;

$syntax = "select * from VARIANT_RESULTS";
$sth = $dch->prepare($syntax);
$sth->execute or die "SQL Error: $DBI::errstr\n";
while (@row7 = $sth->fetchrow_array() ) {	
	foreach my $list(0..$#row7-1){print RES "$row7[$list]\t";}print RES "$row7[$#row7]\n";
}
close RES;


# DISCONNECT FROM THE DATABASE
print "\n\tDISCONNECTING FROM THE DATABASE : $dsn\n\n";
$sth->finish;
$dbh->disconnect();
$dch->disconnect();

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
print "\n\n*********DONE*********\n\n";
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - -T H E  E N D - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
exit;
