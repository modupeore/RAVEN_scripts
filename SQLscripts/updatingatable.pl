#!/usr/bin/perl -w

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# MANUAL FOR insert_results.pl

=pod

=head1 NAME

$0 -- Insert results from tophat and cufflinks analysis into the mysql database : PENGUIN

=head1 SYNOPSIS

summarydb.pl [--help] [--manual]

 
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

# CODE FOR insert_results.pl

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
my $user = 'modupe';
my $passwd = 'penguin123';

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - G L O B A L  V A R I A B L E S- - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# DATABASE VARIABLES
my ($dbh, $sth, $syntax, $row, @row);
my $DE = "DE";
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - M A I N - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# CONNECT TO THE DATABASE
print "\n\n\tCONNECTING TO THE DATABASE : $dsn\n\n";
$dbh = DBI->connect($dsn, $user, $passwd) or die "Connection Error: $DBI::errstr\n";

# SUMMARY OF RESULTS
print "\n\tEXECUTING SELECT STATEMENT ON THE DATABASE TABLES \n";
print "Summary of Results gotten from the database \"$dsn\"\n\n";

#RESULTS_SUMMARY
my $table = "GHI";
$syntax = "select distinct(library_id) from $table";
$sth = $dbh->prepare($syntax);
$sth->execute or die "SQL Error: $DBI::errstr\n";
while ($row = $sth->fetchrow_array() ) {
	my $newrow = $row.$DE;
	my $syntax2 = "update $table SET library_id = ? where library_id = ?";
	my $sth2 = $dbh->prepare($syntax2);
	$sth2->bind_param(1,$newrow);
	$sth2->bind_param(2,$row);
	$sth2->execute ();
}
#print "Number of rows in \"RESULTS_SUMMARY\" table \t:\t @row\n";}
#ACTUAL_FILES
#$syntax = "select count(*) from ACTUAL_FILES";
#$sth = $dbh->prepare($syntax);
#$sth->execute or die "SQL Error: $DBI::errstr\n";
#while (@row = $sth->fetchrow_array() ) {print "Number of rows in \"ACTUAL_FILES\" table   \t:\t @row\n";}
#GENES_FPKM
#$syntax = "select count(distinct library_id) from GENES_FPKM";
#$sth = $dbh->prepare($syntax);
#$sth->execute or die "SQL Error: $DBI::errstr\n";
#while (@row = $sth->fetchrow_array() ) {print "Number of unique libraries in \"GENES_FPKM\" table \t\t:\t @row\n";}
#ISOFORMS_FPKM
#$syntax = "select count(distinct library_id) from ISOFORMS_FPKM";
#$sth = $dbh->prepare($syntax);
#$sth->execute or die "SQL Error: $DBI::errstr\n";
#while (@row = $sth->fetchrow_array() ) {print "Number of unique libraries in \"ISOFORMS_FPKM\" table \t:\t @row\n";}

# DISCONNECT FROM THE DATABASE
print "\n\tDISCONNECTING FROM THE DATABASE : $dsn\n\n";
$dbh->disconnect();

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
print "\n\n*********DONE*********\n\n";
# - - - - - - - - - - - - - - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
exit;
