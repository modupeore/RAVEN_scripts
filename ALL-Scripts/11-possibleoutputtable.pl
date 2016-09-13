#! usr/bin/perl
use File::Basename;

#OBJECTIVE
#tabulating SNPeff output all of it
print "\n\tTabulating SNPeff annotations from the VCF file input\n";

my $input = $ARGV[0];
my $i= 0;
unless(open(FILE,$input)){
	print "File \'$input\' doesn't exist\n";
	exit;
}

my $out = fileparse($input, qr/\.[^.]*(\.vcf)?$/);

my $output = "$out".".submittable";
open(OUT,">$output");

my @file = <FILE>;
chomp @file;
close (FILE);
my $count = 0; my $total=0;
print OUT "CHROM\tPOS\tREF\tALT\tQUAL\tDP\tGT\n";
foreach my $chr (@file){
	if ($chr =~ /^chr/){
		my @chrdetails = split('\t', $chr);
		my $chrIwant = $chrdetails[7];
		my @morechrsplit = split(';', $chrIwant);
		my $secchrIwant = $chrdetails[9];
                my @secmorechrsplit = split(':', $secchrIwant);

		foreach my $Imptchr (@morechrsplit){
			if ($Imptchr =~ m/^DP/) {
				my @addchrsplit = split('=', $Imptchr);
			
				$count++;
				print OUT "$chrdetails[0]\t$chrdetails[1]\t$chrdetails[3]\t$chrdetails[4]\t$chrdetails[5]\t$addchrsplit[1]\t$secmorechrsplit[0]\n";
				
			}
		}
	}
}
close (OUT);
print "\nTotal SNPs\t$count\n";
exit;
