#! usr/bin/perl
use File::Basename;

#OBJECTIVE
#working on the genes output
print "\n\tChanging genes with \"id\" in front of them and counting them\n";

my $input = $ARGV[0];
my $i= 0;
unless(open(FILE,$input)){
	print "File \'$input\' doesn't exist\n";
	exit;
}

my $out = fileparse($input, qr/\.[^.]*(\.vcf)?$/);

my $output = "$out".".editedgenes";
open(OUT,">$output");

my @file = <FILE>;
chomp @file;
close (FILE);
my $count = 0;

foreach my $chr (@file){
	if ($chr =~ /^id(.*)/){
		print OUT "$1\n"; $count++;
	}
	else {
		print OUT "$chr\n";
	}
}

print "\n$count\n";
close (OUT);
exit;
