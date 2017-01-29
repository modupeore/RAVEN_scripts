#!/usr/bin/perl
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# MANUAL for extracting stuff from the database
#10/26/2015

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - U S E R  V A R I A B L E S- - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
use strict;
use DBI;
use lib '/home/modupe/SCRIPTS/SUB';
use routine;
use passw;
my ($dbh, $sth);

# DATABASE ATTRIBUTES
my $basepath = "/home/modupe/TranscriptAtlas";
my $finalpath = "/home/modupe/public_html/FBTranscriptAtlas";
`mkdir -p $basepath`;
my $statusfile = "$basepath/mystatus.txt";
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - G L O B A L  V A R I A B L E S- - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
my %HashDirectory;
# CONNECT TO THE DATABASE
print "\n\n\tCONNECTING TO THE DATABASE \n\n";
$dbh = mysql();
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - M A I N - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CONNECT
#GETTING ALL THE LIBRARIES FROM THE DATABASE.
my $libraries;
print "Extracting library ids that don't have nosql\n";
my $libras = "select a.library_id from vw_libraryinfo a join transcripts_summary b on ((a.library_id = b.library_id)) join variants_summary c on ((a.library_id = c.library_id)) where c.nosql is null and b.status = 'done' and a.species = \"gallus\"";
$sth = $dbh->prepare($libras); $sth->execute or die "SQL Error: $DBI::errstr\n";

while ( my $row = $sth->fetchrow_array() ) {
	$libraries .= $row.",";
}
$libraries = substr($libraries,0,-1);
my @libraries = split (",", $libraries);
print "\nExtracting comtent for library\n";
#PARSING EACH FILE OUT OF THE DATABASE. 
foreach my $file (@libraries){
	opendir (DIR, $basepath) or die "Folder \"$basepath\" doesn't exist\n"; 
	my @Directory = readdir(DIR); close(DIR);
	#pushing each subfolder
	foreach (@Directory){
		if ($_ =~ /^\d*$/){
			$HashDirectory{$_}= $_;}
	}
	unless (exists $HashDirectory{$file}){
		mkdir "$basepath/$file"; print "Working on $file\n";
		$dbh = mysql();
		my $syntax= "select
			b.variant_class,b.zygosity,b.existing_variant,c.consequence,c.gene_id,
			c.gene_name,c.transcript,c.feature,c.gene_type,b.ref_allele,b.alt_allele,a.line,
			a.tissue,b.chrom,c.aminoacid_change,c.codon_change,a.species,
			a.notes,b.quality,a.library_id,a.total_VARIANTS,a.total_SNPS,a.total_INDELS,b.position,
			c.protein_position
			from vw_libraryinfo a join variants_result b 
				on a.library_id = b.library_id 
			join variants_annotation c 
				on b.library_id = c.library_id and b.chrom = c.chrom and b.position = c.position
			where a.species = \"gallus\" and a.library_id = $file";
	 
		$sth = $dbh->prepare($syntax);
		$sth->execute or die "SQL Error: $DBI::errstr\n";
		open(OUT,">$basepath/$file/$file\.txt");
			
		#TABLE FORMAT
		my $i = 0;
		while ( my @row2 = $sth->fetchrow_array() ) {
		$i++;
			foreach my $list (0..$#row2-1){
				if ((length($row2[$list]) < 1) || ($row2[$list] =~ /^\-$/) ){
					$row2[$list] = "NULL";
				}
				if ($list < 18) {
					print OUT "\'$row2[$list]\',";
				}
				else {
					print OUT "$row2[$list],";
				}
			}
			if ((length($row2[$#row2]) < 1) || ($row2[$#row2] =~ /^\-$/) ){
				$row2[$#row2] = "NULL";
			}
			print OUT "$row2[$#row2]\n";
		}
		open(STATUS, ">>$statusfile");
		print STATUS "$file\t$i\n";
		close (OUT);

		print "Importing to FastBit\n";
		my $execute = "ardea -d $finalpath/chickenvariants -m \"
						class:char,
						zygosity:char,
						dbsnp:char,
						consequence:char,
						geneid:char,
						genename:char,
						transcript:char,
						feature:char,
						genetype:char,
						ref:char,
						alt:char,
						line:char,
						tissue:char,
						chrom:char,
						aachange:char,
						codon:char,
						species:key,
						notes:text, 
						quality:double,
						library:int,
						variant:int,
						snp:int,
						indel:int,
						position:int,
						proteinposition:int\" -t $basepath/$file/$file\.txt";
		print $execute,"\n"; system($execute);
		
		#update nosql
		$dbh = mysql();
		$syntax = "update variants_summary set nosql = 'done' where library_id = $file";
		$sth->prepare($syntax); $sth->execute or die "$DBI::Errstr Failed to update variant summary\n";
		print "Done with $file\n";


	}
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
print "\n\n*********DONE*********\n\n";
# - - - - - - - - - - - - - - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
exit;

