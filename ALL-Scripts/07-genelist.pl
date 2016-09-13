#! usr/bin/perl
use File::Basename;

#OBJECTIVE
#tabulating from SNPeff "table" initially created, to a better output.
print "\n\t Creating a Gene list from the Table file for Downstream analysis\n";

my $input = $ARGV[0];
my $i= 0;
unless(open(FILE,$input)){
	print "File \'$input\' doesn't exist\n";
	exit;
}

my $out = fileparse($input, qr/\.[^.]*(\.vcf)?$/);

my $output = "$out".".annotable";
my $output2 = "$out".".genes";
open(OUT,">$output");
open(OUT2,">$output2");
my @file = <FILE>;
chomp @file;
close (FILE);
my $count = 0; my $newcount=0;
print OUT "COUNT\tCHR\tPOS\tR-A\tEFFECT\tEFFECT_IMPACT\tFUNC_CLASS\tCODON_CHANGE\tAA_CHANGE\tGENE_NAME\n";
foreach my $chr (@file){
	if ($chr =~ /^\d/){
		my @chrdetails = split('\t', $chr);
		my $chrIwant = $chrdetails[5];
		my @morechrsplit = split('\(', $chrIwant);
		@finalchrsplit = split('\|',$morechrsplit[1]); # print "$morechrsplit[1]\n@finalchrsplit\n"; die;
		$count++;;
		print OUT "$chrdetails[0]\t$chrdetails[1]\t$chrdetails[2]\t$chrdetails[3]-$chrdetails[4]\t$morechrsplit[0]";
                foreach my $newstart (0..3){print OUT "\t$finalchrsplit[$newstart]";}
                print OUT "\t$finalchrsplit[4]\n";
		$Hashdetails{$finalchrsplit[4]} = 0;
	}
}
close (OUT);
foreach my $newkey (keys %Hashdetails){
	if (length $newkey >=1){
		print OUT2 "$newkey\n"; $newcount++;
	}
}
close (OUT2);
print "\nTotal lines = $count\nTotal number of genes = $newcount\n";
exit;
