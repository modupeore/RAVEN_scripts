#!/usr/bin/perl

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# MANUAL FOR databaseVariants.pl

=pod

=head1 NAME

$0 -- Insert results from variation data into the mysql database : PENGUIN and CHICKENSNPS

=head1 SYNOPSIS

databaseVariants.pl --in1 FOLDER -in2 [--help] [--manual]

=head1 DESCRIPTION

Insert variants for each library as well as annotation information into the database and 
	store as a table.
 
=head1 OPTIONS

=over 3

=item B<-i, -1, -a, --in, --folder>=FOLDER

Variant folder (e.g. the folder with the libraries to insert into the database).  (Required)

=item B<-h, --help>

Displays the usage message.  (Optional) 

=item B<-man, --manual>

Displays full manual.  (Optional) 

=back

=head1 DEPENDENCIES

Requires the following Perl libraries (all standard in most Perl installs).
   DBI
   File::Basename
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
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.  
This is free software: you are free to change and redistribute it.  
There is NO WARRANTY, to the extent permitted by law.  

Please acknowledge author and affiliation in published work arising from this script's 
usage <http://bioinformatics.udel.edu/Core/Acknowledge>.

=cut

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - U S E R  V A R I A B L E S- - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#code for databaseVariants.pl

use strict; use DBI; use DBD::mysql; use Getopt::Long; use Pod::Usage; use File::Basename;

# CREATING LOG FILES
my $std_out = '/home/modupe/.LOG/specificvt-'.`date +%m-%d-%y_%T`; chomp $std_out; $std_out = $std_out.'.log';
my $std_err = '/home/modupe/.LOG/specificvt-'.`date +%m-%d-%y_%T`; chomp $std_err; $std_err = $std_err.'.err';

open(STDOUT, '>', "$std_out") or die "Log file doesn't exist";
open(STDERR, '>', "$std_err") or die "Error file doesn't exist";


#CREATING EMAIL NOTIFICATION
my $notification = '/home/modupe/.LOG/note.txt';
my $email = 'amodupe@udel.edu';
open (NOTE, ">$notification"); 
print NOTE "Subject: Running specificDatabaseVariants.pl!!!\n\nName of files\n\t$std_out\n\t$std_err\n"; 
system "sendmail $email < $notification";
close NOTE;
system "rm -rf $notification";

#ARGUMENTS
my($folder, $library, $help, $manual);
GetOptions (
                                "1|a|in1|folder=s"	=>	\$folder,
                                "2|b|in2|lib=s"	=>	\$library,
                                "h|help"        =>      \$help,
                                "man|manual"	=>      \$manual );

# VALIDATE ARGS
pod2usage( -verbose => 2 )  if ($manual);
pod2usage( -verbose => 1 )  if ($help);
pod2usage( -msg  => "ERROR!  Required argument -1 (input folder and/or library no) not found. \n", 
-exitval => 2, -verbose => 1)  if (! $folder | ! $library );

# DATABASE ATTRIBUTES
my $user = 'modupe';
my $passwd = 'penguin123';
my $dsn = 'dbi:mysql:CHICKENSNPS';

# HASHtables for different works
my %VCFhash; my %VEPhash; my %galgal;

my $galgalgenes = "/home/modupe/DATA/genelist/Gallus_gallus.gene_info";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - G L O B A L  V A R I A B L E S- - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# DATABASE VARIABLES
my ($dbh, $sth, $syntax, $row, @row);

# RESULTS_HASH
my %Hashresults; my $number=0;


# CONNECT TO THE DATABASE

print "\n\n\tCONNECTING TO THE DATABASE : $dsn\n\n";
$dbh = DBI->connect($dsn, $user, $passwd) or die "Connection Error: $DBI::errstr\n";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - M A I N - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#working on the entrez genes list and information
&GALGAL($galgalgenes);

#processing each fastq file
print "\n. . .\n\tWORKING ON $library\n\n";
#processing number
$library =~ /^library-(\d*)$/;
my $code="DE";
my $libnumber = "$1$code";
my $information = "$folder/$library";
#VCF file
my $VCFinput = $information."_DP5.vcf";

#VEP file
my $VEPinput = $information."_VEP.txt";
my $VEPhtml = $information."_VEP.txt_summary.html";

#OUTPUT file
my $output = $information.".table";
open(OUT,">$output");

#other variant files
my $GTFinput = $information."_GTF.vcf";
my $GTFhtml = $information."_GTF.html";
#variant_files table
print ". . .\n\tINSERTING INTO THE DATABASE IN \"VARIANT_FILES\" TABLE\n\n";
$sth = $dbh->prepare("insert into VARIANT_FILES ( library_id, variant_file, variant_vep, vep_html, variant_gtf, 
gtf_html, variant_table ) values (?,?,?,?,?,?,?)");
$sth ->execute($libnumber, $VCFinput, $VEPinput, $VEPhtml, $GTFinput, $GTFhtml, $output);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - M A I N  W O R K F L O W - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
my $total = 0;

print OUT "#CHR\tPOS\tREF\tALT\tQUAL\t";
print OUT "ENSEMBL-GENE\tENTREZ-GENE\tSYMBOL\t";
print OUT "TRANSCRIPT\tFEATURE\tTYPE OF GENE\tCONQ\t";
print OUT "PROTEIN\tPROTEIN POSITION\tCODON CHANGE\tZYGOSITY\tdbSNP\n";
my $date = `date +%Y-%m-%d`;

#initializing the hash tables . . .
%VCFhash = ();
%VEPhash = ();
#running through subroutines . . . 
&VARIANT($VCFinput);
&VEP($VEPinput);
my $itsnp = 0; my $itindel = 0;
#printting to output table & variant_results table
print ". . .\n\tINSERTING INTO THE DATABASE IN \"VARIANT_RESULTS\" TABLE\n\n";
foreach my $abc (sort keys %VCFhash) {
	foreach my $def (sort {$a <=> $b} keys %{ $VCFhash{$abc} }) {
		my @vcf = split('\|', $VCFhash{$abc}{$def});
		my @vep = split('\|', $VEPhash{$abc}{$def});
		$total++;
		if ($vcf[3] =~ /,/){
			my $first = split(",",$vcf[3]);
			if (length $vcf[2] == length $first){
				$itsnp++;
			}
			else {
				$itindel++;
			}
		}
		elsif (length $vcf[2] == length $vcf[3]){
			$itsnp++;
		}
		else {
			$itindel++;
		}
		print OUT "$abc\t$def\t$vcf[2]\t$vcf[3]\t$vcf[4]\t";
		print OUT "$vep[4]\t$vep[11]\t$vep[12]\t";
		print OUT "$vep[5]\t$vep[3]\t$vep[13]\t$vep[6]\t";
		print OUT "$vep[8]\t$vep[7]\t$vep[9]\t$vcf[5]\t$vep[10]\n";

		$sth = $dbh->prepare("insert into VARIANT_RESULTS ( library_id, chrom, position, ref_allele, alt_allele, quality, 
		ensembl_gene, entrez_gene, gene_name, transcript, feature, gene_type,
		consequence, protein, protein_position, codon_change, zygosity, existing_variant ) values
		(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)");
		$sth ->execute($libnumber, $abc, $def, $vcf[2], $vcf[3], $vcf[4], $vep[4], $vep[11], $vep[12],
		$vep[5], $vep[3], $vep[13], $vep[6], $vep[8], $vep[7], $vep[9], $vcf[5], $vep[10]);
	}
}
close (OUT);
#variant_summary table
print ". . .\n\tINSERTING INTO THE DATABASE IN \"VARIANT_SUMMARY\" TABLE\n\n";
$sth = $dbh->prepare("insert into VARIANT_SUMMARY ( library_id, total_VARIANTS, total_SNPS, total_INDELS, date ) 
values (?,?,?,?,?)");
$sth ->execute($libnumber, $total, $itsnp, $itindel, $date);

#VARIANT_FILES table updating status column with 'done'
$sth = $dbh->prepare("update VARIANT_FILES set status='done' where library_id = ?");
$sth ->execute( $libnumber );

# DISCONNECT FROM THE DATABASE
print "\n\tDISCONNECTING FROM THE DATABASE : $dsn\n\n";
$dbh->disconnect();

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
print "\n\n*********DONE*********\n\n";
exit;
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - S U B R O U T I N E S - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#FORMATTING SNPEFF
sub VARIANT {
	unless(open(FILE,$_[0])){
		print "File \'$_[0]\' doesn't exist\n";
		exit;
	}
	my $verd;
	my @file = <FILE>;
	chomp @file;
	close (FILE);
	foreach my $chr (@file){
		unless ($chr =~ /^#/){
			my @chrdetails = split('\t', $chr);
			#removing the random chromosomes (scaffolds) - because no vital information can be found for them.
			if ($chrdetails[0] !~ /^chr.*random$/){
				my $chrIwant = $chrdetails[7];
				my @morechrsplit = split(';', $chrIwant);
				my $chrUwant = $chrdetails[9];
				my @morechrsplit2 = split(':', $chrUwant);
				if ($morechrsplit2[0] eq '0/1'){
					$verd = "heterozygous";
				}
				if ($morechrsplit2[0] eq '1/1'){
					$verd = "homozygous";
				}
				if ($morechrsplit2[0] eq '1/2'){
					$verd = "heterozygous alternate";
				}
				$VCFhash{$chrdetails[0]}{$chrdetails[1]} = "$chrdetails[0]|$chrdetails[1]|$chrdetails[3]|$chrdetails[4]|$chrdetails[5]|$verd";
			}
		}
	}
}

#parsing VEP file
sub VEP {
	unless(open(FILE,$_[0])){
		print "File \'$_[0]\' doesn't exist\n";
		exit;
	}
	%VEPhash = ();
	my @file = <FILE>;
	chomp @file;
	close (FILE);
	my %Geneinfo = ''; my %Transcriptinfo = '';
	my %Conqinfo = ''; my %Altinfo = '';
	my %Featureinfo = ''; my %Proinfo = '';
	my %Prochangeinfo = ''; my %Codoninfo = '';
	my %dbSNPinfo = ''; my %locinfo = '';
	my %Entrezinfo = ''; my %Galgal1info = '';
	my %Galgal2info = ''; my %Galgal3info = '';
	my $position; my $allprotein;

	foreach my $chr (@file){
		unless ($chr =~ /^#/){
			my @chrdetails = split('\t', $chr);
		
		#GENE - 1
			if (exists $Geneinfo{$chrdetails[0]}){
				if ($Geneinfo{$chrdetails[0]} eq $chrdetails[3]){
					$Geneinfo{$chrdetails[0]} = $chrdetails[3];
				}
				else {
					my $tempdetails = "$Geneinfo{$chrdetails[0]},$chrdetails[3]";
					$Geneinfo{$chrdetails[0]} = $tempdetails;
				}
			}
			else {	
				$Geneinfo{$chrdetails[0]} = $chrdetails[3];
			}
		
		#TRANSCRIPT - 2
			if (exists $Transcriptinfo{$chrdetails[0]}){
				if ($Transcriptinfo{$chrdetails[0]} eq $chrdetails[4]){
					$Transcriptinfo{$chrdetails[0]} = $chrdetails[4];
				}
				else {
					my $temp2details = "$Transcriptinfo{$chrdetails[0]},$chrdetails[4]";
					$Transcriptinfo{$chrdetails[0]} = $temp2details;
				}
			}
			else {	
				$Transcriptinfo{$chrdetails[0]} = $chrdetails[4];
			}
		
		#CONSEQUENCE - 3
			if (exists $Conqinfo{$chrdetails[0]}){
				if ($Conqinfo{$chrdetails[0]} eq $chrdetails[6]){
					$Conqinfo{$chrdetails[0]} = $chrdetails[6];
				}
				else {
					my $temp3details = "$Conqinfo{$chrdetails[0]},$chrdetails[6]";
					$Conqinfo{$chrdetails[0]} = $temp3details;
				}			
			}
			else {	
				$Conqinfo{$chrdetails[0]} = $chrdetails[6];
			}
		
		#ALTERNATIVE ALLELE - 4
			if (exists $Altinfo{$chrdetails[0]}){
				if ($Altinfo{$chrdetails[0]} eq $chrdetails[2]){
					$Altinfo{$chrdetails[0]} = $chrdetails[2];
				}
				else {
					my $temp4details = "$Altinfo{$chrdetails[0]},$chrdetails[2]";
					$Altinfo{$chrdetails[0]} = $temp4details;
				}			
			}
			else {	
				$Altinfo{$chrdetails[0]} = $chrdetails[2];
			}		
				
		#FEATURE TYPE - 5
			if (exists $Featureinfo{$chrdetails[0]}){
				if ($Featureinfo{$chrdetails[0]} eq $chrdetails[5]){
					$Featureinfo{$chrdetails[0]} = $chrdetails[5];
				}
				else {
					my $temp5details = "$Featureinfo{$chrdetails[0]},$chrdetails[5]";
					$Featureinfo{$chrdetails[0]} = $temp5details;
				}			
			}
			else {	
				$Featureinfo{$chrdetails[0]} = $chrdetails[5];
			}
				
		#PROTEIN POSITION - 6
			if (exists $Proinfo{$chrdetails[0]}){
				if ($Proinfo{$chrdetails[0]} eq $chrdetails[9]){
					$Proinfo{$chrdetails[0]} = $chrdetails[9];
				}
				else {
					my $temp6details = "$Proinfo{$chrdetails[0]},$chrdetails[9]";
					$Proinfo{$chrdetails[0]} = $temp6details;
				}			
			}
			else {	
				$Proinfo{$chrdetails[0]} = $chrdetails[9];
			}
		
		#PROTEIN CHANGE - 7
			if (exists $Prochangeinfo{$chrdetails[0]}){
				if ($Prochangeinfo{$chrdetails[0]} eq $chrdetails[10]){
					$Prochangeinfo{$chrdetails[0]} = $chrdetails[10];
				}
				else {
					my $temp7details = "$Prochangeinfo{$chrdetails[0]},$chrdetails[10]";
					$Prochangeinfo{$chrdetails[0]} = $temp7details;
				}			
			}
			else {	
				$Prochangeinfo{$chrdetails[0]} = $chrdetails[10];
			}
		
		#CODON - 8
			if (exists $Codoninfo{$chrdetails[0]}){
				if ($Codoninfo{$chrdetails[0]} eq $chrdetails[11]){
					$Codoninfo{$chrdetails[0]} = $chrdetails[11];
				}
				else {
					my $temp8details = "$Codoninfo{$chrdetails[0]},$chrdetails[11];";
					$Codoninfo{$chrdetails[0]} = $temp8details;
				}			
			}
			else {	
				$Codoninfo{$chrdetails[0]} = $chrdetails[11];
			}
		
		#dbSNP - 9
			if (exists $dbSNPinfo{$chrdetails[0]}){
				if ($dbSNPinfo{$chrdetails[0]} eq $chrdetails[12]){
					$dbSNPinfo{$chrdetails[0]} = $chrdetails[12];
				}
				else {
					my $temp9details = "$dbSNPinfo{$chrdetails[0]},$chrdetails[12];";
					$dbSNPinfo{$chrdetails[0]} = $temp9details;
				}			
			}
			else {	
				$dbSNPinfo{$chrdetails[0]} = $chrdetails[12];
			}
		
		#location - 10
			if (exists $locinfo{$chrdetails[0]}){
				if ($locinfo{$chrdetails[0]} eq $chrdetails[1]){
					$locinfo{$chrdetails[0]} = $chrdetails[1];
				}
				else {
					my $temp10details = "$locinfo{$chrdetails[0]},$chrdetails[1];";
					$locinfo{$chrdetails[0]} = $temp10details;
				}			
			}
			else {	
				$locinfo{$chrdetails[0]} = $chrdetails[1];
			}			
		#GALGAL - 11
			my @identity = split("\=", $galgal{$chrdetails[3]});
			#1
			if (exists $Galgal1info{$chrdetails[0]}){
				if ($Galgal1info{$chrdetails[0]} eq $identity[0]){
					$Galgal1info{$chrdetails[0]} = $identity[0];
				}
				else {
					my $temp12Adetails = "$Galgal1info{$chrdetails[0]},$identity[0]";
					$Galgal1info{$chrdetails[0]} = $temp12Adetails;
				}
			}
			else {	
				$Galgal1info{$chrdetails[0]} = $identity[0];
			}
			#2
			if (exists $Galgal2info{$chrdetails[0]}){
				if ($Galgal2info{$chrdetails[0]} eq $identity[1]){
					$Galgal2info{$chrdetails[0]} = $identity[1];
				}
				else {
					my $temp12Bdetails = "$Galgal2info{$chrdetails[0]},$identity[1]";
					$Galgal2info{$chrdetails[0]} = $temp12Bdetails;
				}
			}
			else {	
				$Galgal2info{$chrdetails[0]} = $identity[1];
			}
			#3
			if (exists $Galgal3info{$chrdetails[0]}){
				if ($Galgal3info{$chrdetails[0]} eq $identity[2]){
					$Galgal3info{$chrdetails[0]} = $identity[2];
				}
				else {
					my $temp12Cdetails = "$Galgal3info{$chrdetails[0]},$identity[2]";
					$Galgal3info{$chrdetails[0]} = $temp12Cdetails;
				}
			}
			else {	
				$Galgal3info{$chrdetails[0]} = $identity[2];
			}	
		}
	}			
	foreach my $i (keys %locinfo){
		if ($locinfo{$i} =~ /^(.*)\:(\d*)\-*\d*$/){
			my $tag = 0;
			my @splits = split("\,",$Altinfo{$i});
			foreach (@splits){
				if ($_ eq "-"){$tag = 1;}
			}
			if($tag == 1){
				$position = $2 - 1; $tag = 0; 
			}
			else {$position = $2;}
		my $chromosome = "chr$1";
		$VEPhash{$chromosome}{$position} = "$chromosome|$position|$Altinfo{$i}|$Featureinfo{$i}|$Geneinfo{$i}|$Transcriptinfo{$i}|$Conqinfo{$i}|$Proinfo{$i}|$Prochangeinfo{$i}|$Codoninfo{$i}|$dbSNPinfo{$i}|$Galgal1info{$i}|$Galgal2info{$i}|$Galgal3info{$i}";
		}
	}
}

#parsing Gallus_gallus.gene_info file 
sub GALGAL {
	unless(open(FILE,$_[0])){
		print "File \'$_[0]\' doesn't exist\n";
		exit;
	}
	my @file = <FILE>;
	chomp @file;
	close (FILE);
	my %GalgalA = ''; my %GalgalB = '';
	my %GalgalC = ''; my %GalgalD = '';
	foreach my $line (@file){
		unless ($line =~ /^#/){
			my @details = split('\t', $line);
			my @ensemblid = split("embl\:",$details[5]);
			my $formerensembl = $ensemblid[$#ensemblid];
			$formerensembl =~ /^(ENSGALG\d*)/;
			my $idensembl = $1;
			if ($idensembl =~ /^ENSGAL.*/){
				#GeneID
				if(exists $GalgalA{$idensembl}){
					unless ($GalgalA{$idensembl} eq $details[1]){
						$GalgalA{$idensembl} = "$GalgalA{$idensembl},$details[1]";
					}
				}
				else {
					$GalgalA{$idensembl} = $details[1];
				}
				#Symbol
				if(exists $GalgalB{$idensembl}){
					unless ($GalgalB{$idensembl} eq $details[2]){
						$GalgalB{$idensembl} = "$GalgalB{$idensembl},$details[2]";
					}
				}
				else {
					$GalgalB{$idensembl} = $details[2];
				}
				#Type of gene
				if(exists $GalgalC{$idensembl}){
					unless ($GalgalC{$idensembl} eq $details[9]){
						$GalgalC{$idensembl} = "$GalgalC{$idensembl},$details[9]";
					}
				}
				else {
					$GalgalC{$idensembl} = $details[9];
				}
			}
		}
	}	
	foreach (keys %GalgalA){
		$galgal{$_} = "$GalgalA{$_}=$GalgalB{$_}=$GalgalC{$_}";
	}
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - -T H E  E N D - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


