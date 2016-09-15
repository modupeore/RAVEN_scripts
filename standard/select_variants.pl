#!/usr/bin/perl
use File::Basename;

#OBJECTIVE
#Selecting variants from GATK: VariantFiltration. 
# & saves the header information in another file

my $usage="$0
Only needs the Variant Filtered vcf file
\t<Filtered file>
\t<outputdir>
";

@ARGV==2 or die $usage;

my $input = $ARGV[0];
my $wkdir = $ARGV[1];

unless(open(FILE,$input)){
	print "File \'$input\' doesn't exist\n";
	exit;
}

my $out = fileparse($input, qr/(\.vcf)?$/);

my $output = "$out"."_PASS.vcf";
open(OUT,">$wkdir/$output");
my $output2 = "$out"."_header.vcf";
open(OUT2,">$wkdir/$output2");

my @file = <FILE>;
chomp @file;
close (FILE);
my ($count,$fail) = (0,0);
foreach my $chr (@file){
	if ($chr =~ /^chr/){
		my @chrdetails = split('\t', $chr);
		my $chrIwant = $chrdetails[6];
		if ($chrIwant eq "PASS" ) {
			print OUT "$chr\n"; $count++;
		} else {$fail++;}
	}
	else {
		print OUT "$chr\n";
		print OUT2 "$chr\n";
	}
}
close (OUT); close (OUT2);
print "\tPass filter : $count\n";
print "\tFailed filter : $fail\n";
exit;
