#!/usr/bin/perl

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# MANUAL FOR DBcleaner.pl

=pod

=head1 NAME

$0 -- Clean incomplete results from the mysql database : PENGUIN

=head1 SYNOPSIS

cleaningDB.pl [--help] [--manual]

=head1 DESCRIPTION

Cleans the database of incomplete inserts.
 
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
my($in1,$help,$manual);
GetOptions (	
                                "h|help"        =>      \$help,
                                "man|manual"	=>      \$manual );

# VALIDATE ARGS
pod2usage( -verbose => 2 )  if ($manual);
pod2usage( -verbose => 1 )  if ($help);

#making sure the input file is parsable
#my @temp = split('',$in1); $in1 = undef; my $checking = pop(@temp); push (@temp, $checking); unless($checking eq "/"){ push (@temp,"/")}; foreach(@temp){$in1 .= $_};

# DATABASE ATTRIBUTES
my $dsn = 'dbi:mysql:PENGUIN';
my $user = 'modupe';
my $passwd = 'penguin123';

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - G L O B A L  V A R I A B L E S- - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# FOLDER VARIABLES
my ($top_folder, $cuff_folder);

# RESULTS_HASH
my %Hashresults; my $number=0;

# DATABASE VARIABLES
my ($dbh, $sth, $syntax, $row, @row);
my $code="DE"; 
#PARSING VARIABLES
my @parse; my $len;

# TABLE VARIABLES
my ($lib_id, $total, $mapped, $unmapped, $deletions, $insertions, $junctions, $genes, $isoforms, $prep, $date); #RESULTS_SUMMARY TABLE
my ($raw_reads, $fastqc, $accepted, $unmapped_bam, $deletions_bed, $insertions_bed, $junctions_bed, $skipped_gtf, $transcripts_gtf, $isoforms_fpkm, $genes_fpkm, $run_log); #ACTUAL_FILES TABLE
my ($track, $class, $ref_id, $gene, $gene_name, $tss, $locus, $chrom_no, $chrom_start, $chrom_stop, $length, $coverage, $fpkm, $fpkm_low, $fpkm_high, $fpkm_stat); # GENES_FPKM & ISOFORMS_FPKM TABLE


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - M A I N - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# CONNECT TO THE DATABASE
print "\n\n\tCONNECTING TO THE DATABASE : $dsn\n\n";
$dbh = DBI->connect($dsn, $user, $passwd) or die "Connection Error: $DBI::errstr\n";

#CHECKING TO MAKE SURE NOT "done" FILES ARE REMOVED
$syntax = "select library_id from RESULTS_SUMMARY where status is NULL";
$sth = $dbh->prepare($syntax);
$sth->execute or die "SQL Error: $DBI::errstr\n";
my $coltorem; my $count=0; my @columntoremove;
while ($row = $sth->fetchrow_array() ) {
	$count++; 
	$coltorem .= "AAA$row";
} 
@columntoremove = split("AAA", $coltorem);
for (my $i= 1; $i < $count+1; $i++){
 	my $colzz =  $columntoremove[$i]; 
	print "\nLibrary_$colzz wasn't completed!!!!\n\tBeing removed!!!\n\n";
	#DELETE FROM ACTUAL FILES
        $sth = $dbh->prepare("delete from ACTUAL_FILES where library_id = ?"); 
	$sth->execute( $colzz ); 
        #DELETE FROM GENES_FPKM
        $sth = $dbh->prepare("delete from GENES_FPKM where library_id = ?"); 
	$sth->execute( $colzz );
        #DELETE FROM ISOFORMS_FPKM
        $sth = $dbh->prepare("delete from ISOFORMS_FPKM where library_id = ?"); 
	$sth->execute( $colzz );
        #DELETE FROM RESULTS_SUMMARY
        $sth = $dbh->prepare("delete from RESULTS_SUMMARY where library_id = ?"); 
	$sth->execute( $colzz );
}

# SUMMARY OF RESULTS
print "\n\tEXECUTING SELECT STATEMENT ON THE DATABASE TABLES \n";
print "Summary of Results gotten from \"$in1\" folder in the database \"$dsn\"\n\n";

#RESULTS_SUMMARY
$syntax = "select count(*) from RESULTS_SUMMARY";
$sth = $dbh->prepare($syntax);
$sth->execute or die "SQL Error: $DBI::errstr\n";
while (@row = $sth->fetchrow_array() ) {print "Number of rows in \"RESULTS_SUMMARY\" table \t:\t @row\n";}
#ACTUAL_FILES
$syntax = "select count(*) from ACTUAL_FILES";
$sth = $dbh->prepare($syntax);
$sth->execute or die "SQL Error: $DBI::errstr\n";
while (@row = $sth->fetchrow_array() ) {print "Number of rows in \"ACTUAL_FILES\" table   \t:\t @row\n";}
#GENES_FPKM
$syntax = "select count(*) from GENES_FPKM";
$sth = $dbh->prepare($syntax);
$sth->execute or die "SQL Error: $DBI::errstr\n";
while (@row = $sth->fetchrow_array() ) {print "Number of rows in \"GENES_FPKM\" table \t\t:\t @row\n";}
#ISOFORMS_FPKM
$syntax = "select count(*) from ISOFORMS_FPKM";
$sth = $dbh->prepare($syntax);
$sth->execute or die "SQL Error: $DBI::errstr\n";
while (@row = $sth->fetchrow_array() ) {print "Number of rows in \"ISOFORMS_FPKM\" table \t:\t @row\n";}

# DISCONNECT FROM THE DATABASE
print "\n\tDISCONNECTING FROM THE DATABASE : $dsn\n\n";
$dbh->disconnect();

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
print "\n\n*********DONE*********\n\n";
# - - - - - - - - - - - - - - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
exit;
