#! usr/bin/perl
use File::Basename;

#OBJECTIVE
#tabulating SNPeff output all of it
print "\n\tTabulating VEP text input\n";

my $input = $ARGV[0];
unless(open(FILE,$input)){
	print "File \'$input\' doesn't exist\n";
	exit;
}

my $out = fileparse($input, qr/\.[^.]*(\.txt)?$/);

my $output = "$out".".VEPtable";
open(OUT,">$output");

my @file = <FILE>;
chomp @file;
close (FILE);
my $count = 0; my $total=0;
my %Geneinfo = ''; my %Transcriptinfo = '';
my %Conqinfo = ''; my %Altinfo = '';
my %Featureinfo = ''; my %Proinfo = '';
my %Prochangeinfo = ''; my %Codoninfo = '';
my %dbSNPinfo = ''; my %locinfo = '';
my $position;
print OUT "#CHR\tPOS\tALT\tFEATURE\tGENE\tTRANSCRIPT\tCONQ\tPROTEIN\tPROTEIN CHANGE\tCODON\tdbSNP\n";
foreach my $chr (@file){
	unless ($chr =~ /^#/){
		$total++;
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
		$count++; my $tag = 0;
		my @splits = split("\,",$Altinfo{$i});
		foreach (@splits){
			if ($_ eq "-"){$tag = 1;}
		}
		if($tag == 1){
			$position = $2 - 1; $tag = 0; 
		}
		else {$position = $2;}

		print OUT "chr$1\t$position\t$Altinfo{$i}\t$Featureinfo{$i}\t$Geneinfo{$i}\t$Transcriptinfo{$i}\t";
		print OUT "$Conqinfo{$i}\t$Proinfo{$i}\t$Prochangeinfo{$i}\t$Codoninfo{$i}\t$dbSNPinfo{$i}\n";
	}
}

close (OUT);
print "\nTotal VEP Annotations\t$total\nTotal SNPs\t$count\n";
exit;
