#!/usr/bin/perl
#Extract standard input from the database as standard output
use strict;
use DBI;
use lib '/home/modupe/SCRIPTS/SUB';
use routine;
use passw;
my ($dbh, $sth);
$dbh = mysql();
@ARGV == 1 or die "input syntax\n";
my $syntax = $ARGV[0];

$sth = $dbh->prepare($syntax); $sth->execute or die "SQL Error: $DBI::errstr\n";
my $i = 0;
while ( my @row2 = $sth->fetchrow_array() ) { 
  $i++;
  foreach my $list(0..$#row2-1){
    print "$row2[$list]\t";
  }
  print "$row2[$#row2]\n";
}

print "\nTotal number of records $i\n";
