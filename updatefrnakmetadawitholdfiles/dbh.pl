#!/usr/bin/perl
#update frnak_metadata table with the results from files.txt as a std input
use DBI;
use lib '/home/modupe/SCRIPTS/SUB';
use routine;
use passw;
my ($dbh, $sth);


$dbh = mysql();
$sth = $dbh-> prepare("insert into frnak_metadata (library_id, ref_genome, ann_file, ann_file_ver, stranded, sequences,user ) values (?,?,?,?,?,?,?)"); 

open(IN, $ARGV[0]);
my @content = <IN>; close IN;

shift @content;

foreach (@content){
 my @c = split /\t/;
 print $c[0],"\n";
 $sth->execute($c[0], $c[1], $c[2], $c[3], $c[4], $c[5], $c[6]) or die "$DBI::errstr\n";
}
