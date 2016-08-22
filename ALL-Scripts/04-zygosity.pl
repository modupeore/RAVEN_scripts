#! usr/bin/perl
use File::Basename;

#OBJECTIVE
#zygosity, extracting the homozygouse and heterzygous SNPs into different VCF files.


my $input = $ARGV[0];
my $i= 0;
unless(open(FILE,$input)){
	print "File \'$input\' doesn't exist\n";
	exit;
}

my $out = fileparse($input, qr/(\.vcf)?$/);

my $output1 = "$out"."_hetero.vcf";
open(OUT1,">$output1");

my $output2 = "$out"."_homo.vcf";
open(OUT2,">$output2");

my @file = <FILE>;
chomp @file;
close (FILE);
my $count = 0; my $homo = 0; $hetero=0;;
foreach my $chr (@file){
	if ($chr =~ /^chr/){
		my @chrdetails = split('\t', $chr);
		my $chrIwant = $chrdetails[9];
		my @morechrsplit = split('\:', $chrIwant);

		if ($morechrsplit[0] eq "1/1"){print OUT2 "$chr\n"; $homo++;}
		elsif ($morechrsplit[0] eq "0/1"){print OUT1 "$chr\n"; $hetero++;}
		elsif ($morechrsplit[0] eq "1/2"){print OUT1 "$chr\n"; $hetero++;}
		$count++;
	}
	else {
		print OUT1 "$chr\n";
		print OUT2 "$chr\n";
	}
}
close (OUT1);
close (OUT2);
print "Homozygous = $homo\nHeterozygous = $hetero\nTotal SNPs = $count\n";
exit;
