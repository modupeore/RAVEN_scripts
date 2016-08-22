#!/usr/bin/perl
use strict;
use File::Basename;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#OBJECTIVE


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - U S E R  V A R I A B L E S- - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
my %GTFhash;
my %GFFhash;
my %VEPhash;
		
my $folder = "/home/modupe/NEWTTT/VARIANTS/TEST";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - G L O B A L  V A R I A B L E S- - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

opendir(DIR,$folder) or die "Folder \"$folder\" doesn't exist\n";
my @Directory = readdir(DIR);
close(DIR);

#processing each fastq file
foreach my $library (@Directory){
	if ($library =~ /^library-(\d*)$/) { # parsing library name 
		
		my $information = "$folder/$library/$library";
		#GTF file
		my $GTFinput = $information."_GTF.vcf";

		#GFF file
		my $GFFinput = $information."_GFF.vcf";

		#VEP file
		my $VEPinput = $information."_VEP.txt";

		#OUTPUT file
		my $output = $information.".table";
		open(OUT,">$output");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - M A I N  W O R K F L O W - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		my $total = 0;
		
		print OUT "#CHR\tPOS\tREF\tALT\tVEP-ALT\tQUAL\t";
		print OUT "GTF-GENE\tGFF-GENE\tVEP-GENE\tGTF-TRANSCRIPT\t";
		print OUT "GFF-TRANSCRIPT\tVEP-TRANSCRIPT\tGTF-FEATURE\tGFF-FEATURE\t";
		print OUT "VEP-FEATURE\tGTF-CONQ\tGFF-CONQ\tVEP-CONQ\t";
		print OUT "GTF-IMPACT\tGFF-IMPACT\tGTF-MUTATION\tGFF-MUTATION\tGTF-PROTEIN|PROTEIN POSITION\t";
		print OUT "GFF-PROTEIN|PROTEIN POSITION\tVEP-PROTEIN|PROTEIN POSITION\tGTF-CODON CHANGE\t";
		print OUT "GFF-CODON CHANGE\tVEP-CODON CHANGE\tGTF-ZYGOSITY\tGFF-ZYGOSITY\tVEP-dbSNP\n";
		
		#initializing the hash tables . . .
		%GTFhash = ();
		%GFFhash = ();
		%VEPhash = ();
		&SNPEFF($GTFinput,"gtf");
		&SNPEFF($GFFinput,"gff");
		&VEP($VEPinput);
	
		foreach my $abc (sort keys %GTFhash) {
			foreach my $def (sort {$a <=> $b} keys %{ $GTFhash{$abc} }) {
				my @gtf = split('\|', $GTFhash{$abc}{$def});
				my @gff = split('\|', $GFFhash{$abc}{$def}); #print @gff;die;
				my @vep = split('\|', $VEPhash{$abc}{$def});
				$total++;
				print OUT "$abc\t$def\t$gtf[2]\t$gtf[3]\t$vep[2]\t$gtf[4]\t";
				print OUT "$gtf[10]\t$gff[10]\t$vep[4]\t$gtf[13]\t$gff[13]\t$vep[5]\t";
				print OUT "$gtf[11]|$gtf[12]\t$gff[11]|$gff[12]\t$vep[3]\t$gtf[5]\t$gff[5]\t";
				print OUT "$vep[6]\t$gtf[6]\t$gff[6]\t$gtf[7]\t$gff[7]\t$gtf[9]\t$gff[9]\t";
				print OUT "$vep[8]|$vep[7]\t$gtf[8]\t$gff[8]\t$vep[9]\t$gtf[14]\t$gff[14]\t$vep[10]\n";
			}
		} 
		close (OUT);
	}
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
print "\n\n*********DONE*********\n\n";
exit;
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - S U B R O U T I N E S - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#FORMATTING SNPEFF
sub SNPEFF {
	unless(open(FILE,$_[0])){
		print "File \'$_[0]\' doesn't exist\n";
		exit;
	}
	my ($verd,@finalchrsplit,@secondchrsplit);
	my @file = <FILE>;
	chomp @file;
	close (FILE);
	foreach my $chr (@file){
		unless ($chr =~ /^#/){
			my @chrdetails = split('\t', $chr);
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
			foreach my $Imptchr (@morechrsplit){
				if ($Imptchr =~ m/^EFF/) {
					my @addchrsplit = split('=', $Imptchr);
					@finalchrsplit = split('\(',$addchrsplit[1]);
					@secondchrsplit = split('\|', $finalchrsplit[1]);
				}
			}
			if ($_[1] eq "gtf") {
				$GTFhash{$chrdetails[0]}{$chrdetails[1]} = "$chrdetails[0]|$chrdetails[1]|$chrdetails[3]|$chrdetails[4]|$chrdetails[5]|$finalchrsplit[0]|$secondchrsplit[0]|$secondchrsplit[1]|$secondchrsplit[2]|$secondchrsplit[3]|$secondchrsplit[4]|$secondchrsplit[5]|$secondchrsplit[6]|$secondchrsplit[7]|$verd";
			}
			elsif ($_[1] eq "gff") {
				$GFFhash{$chrdetails[0]}{$chrdetails[1]} = "$chrdetails[0]|$chrdetails[1]|$chrdetails[3]|$chrdetails[4]|$chrdetails[5]|$finalchrsplit[0]|$secondchrsplit[0]|$secondchrsplit[1]|$secondchrsplit[2]|$secondchrsplit[3]|$secondchrsplit[4]|$secondchrsplit[5]|$secondchrsplit[6]|$secondchrsplit[7]|$verd";
			}
			else {
				print "$chrdetails[0]\t$chrdetails[1]\t$chrdetails[3]\t$chrdetails[4]\t$finalchrsplit[0]\t$secondchrsplit[0]\t$secondchrsplit[1]\t$secondchrsplit[2]\t$secondchrsplit[3]\t$secondchrsplit[4]\t$secondchrsplit[5]\t$secondchrsplit[6]\t$secondchrsplit[7]\t$verd\n";
				die;
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
	my $position;

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
		$VEPhash{$chromosome}{$position} = "$chromosome|$position|$Altinfo{$i}|$Featureinfo{$i}|$Geneinfo{$i}|$Transcriptinfo{$i}|$Conqinfo{$i}|$Proinfo{$i}|$Prochangeinfo{$i}|$Codoninfo{$i}|$dbSNPinfo{$i}";
		}
	}
}





# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - -T H E  E N D - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

