#! usr/bin/perl
use File::Basename;

#OBJECTIVE
#filtering DP > 5.
#saving header information in another file

my $input = $ARGV[1];
my $wkdir = $ARGV[0];
my $i= 0;
unless(open(FILE,$input)){
	print "File \'$input\' doesn't exist\n";
	exit;
}

my $out = fileparse($input, qr/(\.vcf)?$/);

my $output = "$out"."_DP5.vcf";
open(OUT,">$wkdir/$output");
my $output2 = "$out"."_header.vcf";
open(OUT2,">$wkdir/$output2");

my @file = <FILE>;
chomp @file;
close (FILE);
my $count = 0;
foreach my $chr (@file){
	if ($chr =~ /^chr/){
		my @chrdetails = split('\t', $chr);
		my $chrIwant = $chrdetails[7];
		my @morechrsplit = split(';', $chrIwant);
		foreach my $Imptchr (@morechrsplit){
			if ($Imptchr =~ m/^DP/) {
				my @addchrsplit = split('=', $Imptchr);
				if ($addchrsplit[1] > 4){print OUT "$chr\n"; $count++;}
			}
		}
	}
	else {
		print OUT "$chr\n";
		print OUT2 "$chr\n";
	}
}
close (OUT); close (OUT2);
print "\n$count\n";
exit;
