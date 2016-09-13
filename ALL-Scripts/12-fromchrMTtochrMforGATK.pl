#! usr/bin/perl
use File::Basename;

#OBJECTIVE
#working on the vcf file to change a chromosome name
print "\n\tChanging chromosome MT \(chrMT\) to \"chrM\"\n";

my $input = $ARGV[0];
my $i= 0;
unless(open(FILE,$input)){
	print "File \'$input\' doesn't exist\n";
	exit;
}

my $out = fileparse($input, qr/\.[^.]*(\.vcf)?$/);

my $output = "$out".".out-vcf";
open(OUT,">$output");

my @file = <FILE>;
chomp @file;
close (FILE);
my $count = 0;

foreach my $chr (@file){
	if ($chr =~ /^chrMT.*/){
		my @thisisit = split('\t', $chr);
		my $ending = $#thisisit;
		print OUT "chrM\t"; $count++;
		foreach my $abc (1..$ending-1){
			print OUT "$thisisit[$abc]\t";
		}
			print OUT "$thisisit[$ending]\n";
	}
	else {
		print OUT "$chr\n"; $count++;
	}
}

print "\n$count\n";
close (OUT);
exit;
