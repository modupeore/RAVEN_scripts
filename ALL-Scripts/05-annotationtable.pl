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

my $output = "$out".".table";
open(OUT,">$output");

my @file = <FILE>;
chomp @file;
close (FILE);
my $count = 0; my $total=0;

foreach my $chr (@file){
	if ($chr =~ /^chr/){
		my @chrdetails = split('\t', $chr);
		my $chrIwant = $chrdetails[7];
		my @morechrsplit = split(';', $chrIwant);
		
		foreach my $Imptchr (@morechrsplit){
			if ($Imptchr =~ m/^EFF/) {
				my @addchrsplit = split('=', $Imptchr);
				@finalchrsplit = split("\,",$addchrsplit[1]);
				my $finale = $#finalchrsplit+1; my $abc = $count++;
				foreach my $newstart (0..$#finalchrsplit-1){print OUT "$abc\t$chrdetails[0]\t$chrdetails[1]\t$chrdetails[3]\t$chrdetails[4]\t$finalchrsplit[$newstart]\n";}
				print OUT "$abc\t$chrdetails[0]\t$chrdetails[1]\t$chrdetails[3]\t$chrdetails[4]\t$finalchrsplit[$#finalchrsplit]\n";
				$total = $total+$finale;
			}
		}
	}
}
close (OUT);
print "\nTotal SNPs\t$count\nTotal Annotations\t$total\n";
exit;
