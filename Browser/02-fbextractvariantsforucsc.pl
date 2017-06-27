#!/usr/bin/perl
use strict;
use DBI;
use lib '/home/modupe/SCRIPTS/SUB';
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
our ($chickenpath, $mousepath, $alligatorpath) = FBPATHS();

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - M A I N - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#GETTING ALL THE LIBRARIES FROM THE DATABASE.
my @thelines = undef;
my %LineHash; my %RealHash;
my $ibis = "/home/modupe/.bin/bin/ibis -d $chickenpath -q '";
my $syntax = "select line,library where line != \"NULL\"' -o temp.txt";

`$ibis$syntax`; #execute ibis
print "\nDone with ibis " ,scalar (localtime), "\n";
open(IN, "<temp.txt");
while (<IN>){
	chomp;
	my ($linename, $libraryno) = split /,/;
	$linename =~ s/^"|"$//g;
	$LineHash{$libraryno} = $linename;
}
close (IN);
print "\nDone with parse 1st dictionary " ,scalar (localtime), "\n";
foreach my $libno (sort {$a <=> $b} keys %LineHash){
	if (exists $RealHash{$LineHash{$libno}}){
		$RealHash{$LineHash{$libno}} = $RealHash{$LineHash{$libno}}.",".$libno;
	} else {
		$RealHash{$LineHash{$libno}} = $libno;
	}
}

print "\nDONE with parse 2nd dictionary " ,scalar (localtime), "\n";
`rm -rf temp.txt`; #removing temp file

foreach (keys %RealHash) {
	print "$_\t$RealHash{$_}\n\n";
}

#PARSING EACH FILE OUT OF THE DATABASE. 
foreach my $sublibraries (keys %RealHash){
  opendir (DIR, $basepath) or die "Folder \"$basepath\" doesn't exist\n"; 
  my @Directory = readdir(DIR); close(DIR);
  #pushing each subfolder
  foreach (@Directory) {
    if ($_ =~ /^(\w.*)\.txt$/) {
      $HashDirectory{$1}= $1;
    }
  }
  unless (exists $HashDirectory{$sublibraries}){
		if (length($sublibraries) > 1) {
			open(STATUS, ">>$statusfile");
			print STATUS "Processing: $sublibraries\tStarting: ",scalar (localtime),"\n";
			close(STATUS);
			open(OUT,">$basepath/$sublibraries\.txt") or die "cant open output\n";
			$syntax= "select library,chrom,position,ref,alt,quality,consequence,genename,geneid,feature,transcript,genetype,proteinposition,aachange,codon,dbsnp,class,zygosity,tissue where library in ($RealHash{$sublibraries}) order by chrom,position' -o newtemp.txt";
	
			`$ibis$syntax`; #running ibis
			
			open(STATUS, ">>$statusfile");
			print STATUS "Midway: $sublibraries\tIBIS complete: ",scalar (localtime),"\n";
			close (STATUS);
    
			print OUT "library_id\tchrom\tposition\tref_allele\talt_allele\tquality\tconsequence\tgene_name\tgene_id\tfeature\ttranscript\tgene_type\tprotein_position\taminoacid_change\tcodon_change\texisting_variant\tvariant_class\tzygosity\ttissue\n";
    
			my $i = 0;
		
			open(IN2, "<newtemp.txt");
			while (<IN2>){
				chomp;
				$i++;
			
				my ($library, $chrom, $position,
					$ref, $alt, $quality,
					$consequence, $genename, $geneid,
					$feature, $transcript, $genetype,
					$pposition, $aachange, $codon,
					$dbsnp, $class, $zygosity,
					$tissue) = split(/\, /, $_, 19);
				#removing the quotation marks from the words
				$chrom = substr($chrom,1,-1); $ref = substr($ref,1,-1); $alt = substr($alt,1,-1);
				$consequence = substr($consequence,1,-1); $genename = substr($genename,1,-1); $geneid = substr($geneid,1,-1);
				$feature = substr($feature,1,-1); $transcript = substr($transcript,1,-1); $genetype = substr($genename,1,-1);
				$aachange = substr($aachange,1,-1); $codon = substr($codon,1,-1); $dbsnp = substr($dbsnp,1,-1);
				$class = substr($class,1,-1); $zygosity = substr($zygosity,1,-1); $tissue = substr($tissue,1,-1);
				if ($genename =~ /NULL/) { $genename = "-"; }
				if ($dbsnp =~ /NULL/) { $dbsnp = " "; }
			
				print OUT "$library\t$chrom\t$position\t$ref\t$alt\t$quality\t$consequence\t$genename\t$geneid\t$feature\t$transcript\t$genetype\t$pposition\t$aachange\t$codon\t$dbsnp\t$class\t$zygosity\t$tissue\n";
			}
		
			open(STATUS, ">>$statusfile");
			print STATUS "End: $sublibraries\t$i\tDone: ",scalar (localtime),"\n";
			close (OUT); close (STATUS);
		}
	}
}
`rm -rf newtemp.txt`; #remove newtemp file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
print "\n\n*********DONE*********\n\n";
# - - - - - - - - - - - - - - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
exit;

