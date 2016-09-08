#!/usr/bin/perl
use DBI; 
use DBD::mysql; 
use Getopt::Long; 
use Pod::Usage; 

# DATABASE ATTRIBUTES
my $dsn = 'dbi:mysql:transcriptatlas';
my $user = 'frnakenstein';
my $passwd = 'maryshelley';

$dbh = DBI->connect($dsn, $user, $passwd) or die "Connection Error: $DBI::errstr\n";
#open(IN, $ARGV[0]);
#$syntax = "select distinct library_id from variants_annotation order by library_id desc limit 5";
#$sth=$dbh->prepare($syntax);
#$sth->execute;
#while ( my @row = $sth->fetchrow_array() ) {
#        foreach my $jj (0..$#row-1){
#                print "$row[$jj]\t";
#        }
#        print "$row[$#row]\n";
#}
$syntax = "select count(*) from isoforms_fpkm  where library_id = 776";
print "\n$syntax\n";
$sth=$dbh->prepare($syntax);
$sth->execute;
while ( my @row = $sth->fetchrow_array() ) {
        foreach my $jj (0..$#row-1){
                print "$row[$jj]\t";
        }
        print "$row[$#row]\n";
}
