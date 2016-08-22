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

Output directory.  (Optional)

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
my($specifics,$dir,$help,$manual);
GetOptions(
                                "1|a|in|in1|list=s"	=>	\$specifics,
                              	"2|b|out|dir=s"        =>      \$dir,
			 	"h|help"        =>      \$help,
                                "man|manual"	=>      \$manual );

# VALIDATE ARGS
pod2usage( -verbose => 2 )  if ($manual);
pod2usage( -verbose => 1 )  if ($help);
pod2usage( -msg  => "ERROR!  Required argument -1 (list of libraries separated by comma) is not included. \n", 
-exitval => 2, -verbose => 1)  if (! $specifics);


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
print "\n\n\tCONNECTING TO THE DATABASE : $dsn\n\n";
$dbh = DBI->connect($dsn, $user, $passwd) or die "Connection Error: $DBI::errstr\n";


#TABLE COLUMNS
my ($gene_id, $chrom_no, $chrom_start, $chrom_stop, $fpkm);
my ($fpkm_conf_low, $fpkm_conf_high, $library_id);

#EMPTY SPACES COUNTER
my $no = 0;
#HASH TABLES
my (%CHROM, %FPKM, %GENES);

# OPENING OUTPUT FILE

my $output = "gene-list.txt";
my $output2 = "chr-list.txt";
if ($dir){
	`mkdir $dir`; 
	open (OUT, ">$dir/$output");
	open (OUTCHROM, ">$dir/$output2");
}
else { 
	open (OUT, ">$output");
	open (OUTCHROM, ">$output2");
}

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
#gene_id, chrom_no, chrom_start, chrom_stop, fpkm, library_id
print "\n\n\tPROCESSING libraries specified : $specifics \n\n";
$syntax = "select gene_id, chrom_no, chrom_start, chrom_stop, fpkm, 
	library_id from GENES_FPKM where library_id in 
	($specifics) ORDER BY gene_id desc;";
$sth = $dbh->prepare($syntax);
$sth->execute or die "SQL Error: $DBI::errstr\n";

while (@row = $sth->fetchrow_array() ) { 
	if (length($row[0]) < 1){$no++; $gene_id = "+$no";}
	else {$gene_id = $row[0];}
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
	$GENES{$realrow}{$gene_id}{$library_id} = $row[0];
}

print "\n\n\tPRODUCING GENE LIST TABLE\n\n";
my ($post, $lib, $unique);
foreach $post (sort keys %CHROM){
	foreach $unique (sort keys %{ $FPKM{$post} }){
		if ($post =~ /^\+\d*/){next;}
		else{
			print OUT "$post\t$CHROM{$post}\t";
			foreach $lib (0..$#headers-1){ 
				if (exists $FPKM{$post}{$unique}{$headers[$lib]}){
					print OUT "$FPKM{$post}{$unique}{$headers[$lib]}\t";
				}
				else {print OUT "null\t";}
			}
			if (exists $FPKM{$post}{$unique}{$headers[$#headers]}){
				print OUT "$FPKM{$post}{$unique}{$headers[$#headers]}\n";
			}
			else {print OUT "null\n";} 
		}
	}
}

print "\n\n\tPRODUCING CHROMSOME LIST TABLE\n\n";
my ($position, $ident, $unique2, $unique3);
foreach $position (sort {$a cmp $b} (keys %GENES)){
	foreach $unique2 (keys %{ $GENES{$position} }){
		foreach $unique3 (keys %{ $GENES{$position}{$unique2} } ) {
			print OUTCHROM "$position\t$GENES{$position}{$unique2}{$unique3}\t";
			foreach $ident (0..$#headers-1){ 
				if (exists $FPKM{$unique2}{$position}{$headers[$ident]}){
					print OUTCHROM "$FPKM{$unique2}{$position}{$headers[$ident]}\t";
				}
				else {print OUTCHROM "null\t";}
			}
			if (exists $FPKM{$unique2}{$position}{$headers[$#headers]}){
				print OUTCHROM "$FPKM{$unique2}{$position}{$headers[$#headers]}\n";
			}
			else {print OUTCHROM "null\n";} 
		}
	}
}
print "\n\n\tKUDOS ALL DONE\n\n";
# DISCONNECT FROM THE DATABASE
close(OUT);
close(OUTCHROM);
print "\n\tDISCONNECTING FROM THE DATABASE : $dsn\n\n";
$dbh->disconnect();

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
print "\n\n*********DONE*********\n\n";
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - -T H E  E N D - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
exit;
