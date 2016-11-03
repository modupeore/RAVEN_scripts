#!/usr/bin/perl
use strict;
use DBI;
use lib '/home/modupe/SCRIPTS/SUB';
use routine;
my ($dbh, $sth);
$dbh = mysql();
my $syntax= $ARGV[0];

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
