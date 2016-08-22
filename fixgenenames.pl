#!/usr/bin/perl

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# MANUAL FOR fixgenenames.pl

=pod

=head1 NAME

$0 -- Insert results from tophat and cufflinks analysis into the mysql database : PENGUIN

=head1 SYNOPSIS

insert_results.pl --in FOLDER [--help] [--manual]

=head1 DESCRIPTION

##
 
=head1 OPTIONS

=over 3

=item B<-i, -1, -a, --in, --folder>=FOLDER

Results folder (e.g. the folder with results to insert into the database).  (Required)

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

# DATABASE ATTRIBUTES
my $dsn = 'dbi:mysql:PENGUIN';
my $user = 'modupe';
my $passwd = 'penguin123';

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - G L O B A L  V A R I A B L E S- - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# DATABASE VARIABLES
my ($dbh, $sth, $syntax, @row);
#my $code="DE"; 

my ($gene_new_name, $library_id, $number, $chrom1, $chrom2, $chrom3);


# CONNECT TO THE DATABASE
print "\n\n\tCONNECTING TO THE DATABASE : $dsn\n\n";
$dbh = DBI->connect($dsn, $user, $passwd) or die "Connection Error: $DBI::errstr\n";

#SQL syntax
print "\n\n\tWorking on the \"at fault\" columns\n\n";
$syntax = "select library_id, gene_short_name, chrom_no, chrom_start, chrom_stop from GENES_FPKM where gene_short_name like \"%,%\" and library_id between 1 and 500";
$sth = $dbh->prepare($syntax);
$sth->execute or die "SQL Error: $DBI::errstr\n";
while (@row = $sth->fetchrow_array() ) { 
	my @first = split('DE', $row[0]);
	$library_id = $first[0];
	$gene_new_name = $row[1];
	$chrom1 = $row[2];
	$chrom2 = $row[3];
	$chrom3 = $row[4];
	if (length($gene_new_name)==20){ 
		my $outfile = "/home/modupe/genesfpkm/$library_id"."genes.fpkm_tracking";
		open(GENES, "<$outfile") or die "Can't open file $outfile\n";
		print "Working on library $library_id and fixing $gene_new_name\n";
		while (<GENES>){
			chomp;
 			my ($track, $class, $ref_id, $gene, $gene_name, $tss, $locus, $length, $coverage, $fpkm, $fpkm_low, $fpkm_high, $fpkm_stat ) = split /\t/;
 			$locus =~ /^(.+)\:(.+)\-(.+)$/;
             		my $chrom_no = $1; my $chrom_start = $2; my $chrom_stop = $3;
             		if ($gene_new_name eq $gene_name && $chrom_no == $chrom1 && $chrom_start == $chrom2 && $chrom_stop == $chrom3){
				print "No need all is well\n";next;
			}
			elsif ($chrom_no == $chrom1 && $chrom_start == $chrom2 && $chrom_stop == $chrom3 && $gene_new_name ne $gene_name){
				my $syntax2 = "update GENES_FPKM set gene_short_name = \"$gene_name\" where library_id = $library_id and chrom_no = \"$chrom1\" and chrom_start = $chrom2 and chrom_stop = $chrom3"; print "$syntax2\n";
 				my $sth2 = $dbh->prepare($syntax2);
     				$sth2->execute ();
				$sth2->finish;
     			}
		}
	}        	
	$number++;	
}
print "\n$number\n";
# DISCONNECT FROM THE DATABASE
print "\n\tDISCONNECTING FROM THE DATABASE : $dsn\n\n";
$sth->finish;
$dbh->disconnect();

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
print "\n\n*********DONE*********\n\n";
# - - - - - - - - - - - - - - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
exit;

