#!/usr/bin/perl
#removing duplicates from annovar results.
use File::Basename;

my $input = $ARGV[0];

unless(open(FILE,$input)){
        print "File \'$input\' doesn't exist\n";
        exit;
}

my $out = fileparse($input, qr/(\.txt)?$/);

my $output = "$out"."_cleaned.txt";
open(OUT,">$output");


my @file = <FILE>;
chomp @file;
close (FILE);
my $count = 0;
my %Details;
my %Altdetails;
my %Refdetails;
foreach my $chr (@file){
        unless ($chr =~ /^\s/){
                my @chrdetails = split('\t', $chr, 6);
                my @newchr = split('\t', $chr);
                $Details{$chrdetails[0]}{$newchr[19]} = $chrdetails[5];
                if (exists $Altdetails{$chrdetails[0]}{$newchr[19]}) {
                	$count++;
                	$Altdetails{$chrdetails[0]}{$newchr[19]} = "$Altdetails{$chrdetails[0]}{$newchr[19]},$chrdetails[4]";
                } else {
                	$Altdetails{$chrdetails[0]}{$newchr[19]} = $chrdetails[4];
                }
                if (exists $Refdetails{$chrdetails[0]}{$newchr[19]}) {
                        $Refdetails{$chrdetails[0]}{$newchr[19]} = "$Refdetails{$chrdetails[0]}{$newchr[19]},$chrdetails[3]";
                } else {
                        $Refdetails{$chrdetails[0]}{$newchr[19]} = $chrdetails[3];
                }
                if (exists $Odadetails{$chrdetails[0]}{$newchr[19]}) {
                        $Odadetails{$chrdetails[0]}{$newchr[19]} = "$Odadetails{$chrdetails[0]}{$newchr[19]},$chrdetails[2]";
                } else {
                        $Odadetails{$chrdetails[0]}{$newchr[19]} = $chrdetails[2];
                }
        }
}

foreach my $a (sort keys %Details){
	foreach my $b (sort keys %{ $Details{$a} }){
		print OUT "$a\t$b\t$Odadetails{$a}{$b}\t$Refdetails{$a}{$b}\t$Altdetails{$a}{$b}\t$Details{$a}{$b}\n";
	}
}

print "Total number of duplicates $count\n";
exit;

