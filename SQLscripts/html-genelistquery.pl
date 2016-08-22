#!/usr/bin/perl

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# MANUAL FOR genelistquery.pl

=pod

=head1 NAME

$0 -- Extract genes information from specified library_id in the mysql database : PENGUIN

=head1 SYNOPSIS

genelistquery.pl [--help] [--manual]

 
=head1 OPTIONS

=over 3

=item B<-1, -a, in>

List of libraries numbers separated by comma (e.g 123,975).  (Required)

=item B<-2, --b, -out>

Output1 .  (Optional)

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

# CODE FOR genelistquery.pl

use strict;
use DBI;
use DBD::mysql;
use Getopt::Long;
use Pod::Usage;

#ARGUMENTS
my($specifics,$output1,$output2,$help,$manual);
GetOptions(
				"1|a|in|in1|list=s"	=>	\$specifics,
				"2|b|out1|output1=s"        =>      \$output1,
				"3|c|out2|output2=s"        =>      \$output2,
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
my ($dbh, $sth, $syntax, @row);

# CONNECT TO THE DATABASE
#open (OUT,">/home/modupe/genelist.txt");
$dbh = DBI->connect($dsn, $user, $passwd) or die "Connection Error: $DBI::errstr\n";


#TABLE COLUMNS
my ($gene_id, $chrom_no, $chrom_start, $chrom_stop, $fpkm);
my ($fpkm_conf_low, $fpkm_conf_high, $library_id);

#EMPTY SPACES COUNTER
my $no = 0;
#HASH TABLES
my (%CHROM, %FPKM, %GENES);

# OPENING OUTPUT FILE
open (OUT, ">$output1");
open (OUTCHROM, ">$output2");


#SPECIFYING LIBRARIES OF INTEREST
#my $specifics = "123,973,974,976,977,979";
my @headers = split("\,", $specifics);

# HEADER print out
print OUT "GENE\tCHROM POSITION\t";
print OUTCHROM "CHROM POSITION\tGENE\t";
foreach my $name (0..$#headers-1){
	print OUT "library_$headers[$name]\t";
	print OUTCHROM "library_$headers[$name]\t";
}
print OUT "library_$headers[$#headers]\n";
print OUTCHROM "library_$headers[$#headers]\n";
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - M A I N  W O R K F L O W - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#TABLE COLUMNS
#gene_short_name, chrom_no, chrom_start, chrom_stop, fpkm, library_id
$syntax = "select gene_short_name, chrom_no, chrom_start, chrom_stop, fpkm, 
	library_id from GENES_FPKM where library_id in 
	($specifics) ORDER BY gene_id desc;";
$sth = $dbh->prepare($syntax);
$sth->execute or die "SQL Error: $DBI::errstr\n";

while (@row = $sth->fetchrow_array() ) { 
	$gene_id = $row[0];
	if ($row[1] =~ /^chr(.*)$/) {$chrom_no = $1;} 
	else {$chrom_no = $row[1]};
	$chrom_start = $row[2];
	$chrom_stop = $row[3];
	$fpkm = $row[4];
	my @initial = split('DE',$row[5]);
	$library_id = $initial[0];
	my $realrow = $chrom_no."_".$chrom_start."_".$chrom_stop;
	$FPKM{$gene_id}{$realrow}{$library_id} = $fpkm;
	$CHROM{$gene_id} = $realrow;
	$GENES{$realrow} = $gene_id;
}
my ($post, $lib);
foreach $post (sort keys %CHROM){
	if ($post =~ /^\-.*$/ || $post =~ /^\d.*$/){next;}
	else{
		print OUT "$post\t$CHROM{$post}\t";
		foreach $lib (0..$#headers-1){ 
			if (exists $FPKM{$post}{$CHROM{$post}}{$headers[$lib]}){
				print OUT "$FPKM{$post}{$CHROM{$post}}{$headers[$lib]}\t";
			}
			else {
				print OUT "null\t";
			}
		}
		if (exists $FPKM{$post}{$CHROM{$post}}{$headers[$#headers]}){
			print OUT "$FPKM{$post}{$CHROM{$post}}{$headers[$#headers]}\n";
		}
		else {
			print OUT "null\n";
		} 
	}
}
my ($position, $ident);
foreach $position (sort keys %GENES){
	print OUTCHROM "$position\t$GENES{$position}\t";
	foreach $ident (0..$#headers-1){ 
		if (exists $FPKM{$GENES{$position}}{$position}{$headers[$ident]}){
			print OUTCHROM "$FPKM{$GENES{$position}}{$position}{$headers[$ident]}\t";
		}
		else {print OUTCHROM "null\t";}
	}
	if (exists $FPKM{$GENES{$position}}{$position}{$headers[$#headers]}){
		print OUTCHROM "$FPKM{$GENES{$position}}{$position}{$headers[$#headers]}\n";
	}
	else {print OUTCHROM "null\n";} 
}
# DISCONNECT FROM THE DATABASE
close(OUT);
close(OUTCHROM);
$dbh->disconnect();

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - -T H E  E N D - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
exit;
