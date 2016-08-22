#! usr/bin/perl
use File::Basename;

#OBJECTIVE
#getting header information only


my $input = $ARGV[0];
my $i= 0;
unless(open(FILE,$input)){
	print "File \'$input\' doesn't exist\n";
	exit;
}

my $out = fileparse($input, qr/(\.vcf)?$/);

my $output = "$out"."_header.vcf";
open(OUT,">$output");

my @file = <FILE>;
chomp @file;
close (FILE);
my $count = 0;
foreach my $chr (@file){
	if ($chr !~ /^chr/){
		print OUT "$chr\n";
	}
}
close (OUT);
exit;
