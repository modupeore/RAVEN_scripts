#!/usr/bin/perl
use strict;
use DBI;
use lib '/home/modupe/SCRIPTS/SUB';
use passw;
use routine;

#This is to extract variants and all the pertinent information into a txt file based on the Chicken line
#Output is stored in the folder specified

# DATABASE ATTRIBUTES
my $basepath = $ARGV[0]; 
#"/home/modupe/public_html/DBoutput";
my $statusfile = "$basepath/mystatus.txt";

open(STATUS,">>$statusfile") or die "$basepath does not exist\n"; close (STATUS);
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - G L O B A L  V A R I A B L E S- - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
my ($dbh, $sth); my %HashDirectory;
$dbh = mysql();
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - M A I N - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#GETTING ALL THE LIBRARIES FROM THE DATABASE.
my @thelines = undef;
my %LineHash;
my $line = "select distinct line from vw_libraryinfo where species = \"gallus\" and line is not null;";
$sth = $dbh->prepare($line); $sth->execute or die "SQL Error: $DBI::errstr\n";
while ( my @row = $sth->fetchrow_array() ) {
  push @thelines, @row;
}
foreach my $sline (@thelines) {
  my $libras = "select library_id from vw_libraryinfo where line = \"".$sline."\";";
  $sth = $dbh->prepare($libras); $sth->execute or die "SQL Error: $DBI::errstr\n";
  my $libraries = undef;
  while ( my $row = $sth->fetchrow_array() ) {
	$libraries .= $row.",";
  }
  $libraries = substr($libraries,0,-1);
  $LineHash{$sline} = $libraries;
}

#PARSING EACH FILE OUT OF THE DATABASE. 
foreach my $sublibraries (keys %LineHash){
  opendir (DIR, $basepath) or die "Folder \"$basepath\" doesn't exist\n"; 
  my @Directory = readdir(DIR); close(DIR);
  #pushing each subfolder
  foreach (@Directory) {
    if ($_ =~ /^(\w.*)\.txt$/) {
      $HashDirectory{$1}= $1;
    }
  }
  unless (exists $HashDirectory{$sublibraries}){if (length($sublibraries) > 1) {
    open(STATUS, ">>$statusfile");
    print STATUS "Processing: $sublibraries\tStarting: ",scalar (localtime),"\n";
    close(STATUS);
    open(OUT,">$basepath/$sublibraries\.txt") or die "cant open output\n";
    my $syntax= "select
        a.library_id, a.chrom, a.position, b.ref_allele, 
	b.alt_allele, b.quality, a.consequence, a.gene_name, 
	a.gene_id, a.feature, a.transcript, a.gene_type, 
	a.protein_position, a.aminoacid_change, a.codon_change, b.existing_variant, 
	b.variant_class, b.zygosity, c.tissue
        from variants_annotation a join variants_result b
                on a.library_id = b.library_id and 
		a.chrom = b.chrom and 
		a.position = b.position
        join vw_libraryinfo c
                on a.library_id = c.library_id
        where a.library_id in ($LineHash{$sublibraries})
	order by a.chrom, a.position";
    $sth = $dbh->prepare($syntax);
    $sth->execute or die "SQL Error: $DBI::errstr\n";
    
    open(STATUS, ">>$statusfile");
    print STATUS "Midway: $sublibraries\tSQL complete: ",scalar (localtime),"\n";
    close (STATUS);
    
    print OUT "library_id\tchrom\tposition\tref_allele\talt_allele\tquality\tconsequence\tgene_name\tgene_id\tfeature\ttranscript\tgene_type\tprotein_position\taminoacid_change\tcodon_change\texisting_variant\tvariant_class\tzygosity\ttissue\n";
    
    my $i = 0;
    while ( my @row2 = $sth->fetchrow_array() ) {
      $i++;
      foreach my $list (0..$#row2-1){
	  print OUT $row2[$list],"\t"
      }
      print OUT "$row2[$#row2]\n";
    }
    
    open(STATUS, ">>$statusfile");
    print STATUS "End: $sublibraries\t$i\tDone: ",scalar (localtime),"\n";
    close (OUT); close (STATUS);
  }}
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
print "\n\n*********DONE*********\n\n";
# - - - - - - - - - - - - - - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
exit;

