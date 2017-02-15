#!/usr/bin/perl
use DBI; 
use DBD::mysql; 
use Getopt::Long; 
use Pod::Usage; 
use lib '/home/modupe/SCRIPTS/SUB';
use passw;
use routine;

my ($dbh, $sth); my %HashDirectory;
$dbh = mysql();

open(IN, $ARGV[0]) or die "Input file is not specified\n";
#my $syntax = "insert into genes_fpkm (library_id, tracking_id, class_code, nearest_ref_id, gene_id, gene_short_name, tss_id, chrom_no, chrom_start, chrom_stop, length, coverage, fpkm, fpkm_conf_low, fpkm_conf_high, fpkm_status ) values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)";
my $syntax = "insert into genes_fpkm (library_id, gene_id, gene_short_name, chrom_no, chrom_start, chrom_stop, coverage, fpkm, fpkm_conf_low, fpkm_conf_high, fpkm_status ) values (?,?,?,?,?,?,?,?,?,?,?)";
my $sth = $dbh->prepare($syntax);
my @details = <IN>;
shift @details;
foreach (@details){
	#print $_;
        chomp;
	my ($lib, $track, $class, $ref_id, $gene, $gene_name, $tss, $chrom, $start, $stop, $length, $coverage, $fpkm, $fpkm_low, $fpkm_high, $fpkm_stat ) = split /\t/;
        $sth ->execute($lib, $gene, $gene_name, $chrom, $start, $stop, $coverage, $fpkm, $fpkm_low, $fpkm_high, $fpkm_stat) or die "\nERROR:\t Complication in genes_fpkm table, consult documentation\n";
}
