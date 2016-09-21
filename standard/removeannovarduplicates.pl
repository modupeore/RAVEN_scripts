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
foreach my $chr (@file){
        unless ($chr =~ /^\s/){
                my @chrdetails = split('\t', $chr, 6);
                my @newchr = split('\t', $chr);
                $Details{$chrdetails[0]}{$chrdetails[1]}{$chrdetails[2]}{$chrdetails[3]} = $chrdetails[5];
                if (exists $Altdetails{$chrdetails[0]}{$chrdetails[1]}{$chrdetails[2]}{$chrdetails[3]}) {
                	$count++;
                	$Altdetails{$chrdetails[0]}{$chrdetails[1]}{$chrdetails[2]}{$chrdetails[3]} = "$Altdetails{$chrdetails[0]}{$chrdetails[1]}{$chrdetails[2]}{$chrdetails[3]},$chrdetails[4]";
                } else {
                	$Altdetails{$chrdetails[0]}{$chrdetails[1]}{$chrdetails[2]}{$chrdetails[3]} = $chrdetails[4];
                }
        }
}

foreach my $a (sort keys %Details){
	foreach my $b (sort keys %{ $Details{$a} }){
		foreach my $c (sort keys %{ $Details{$a}{$b} } ){
			foreach my $d (sort keys %{ $Details{$a}{$b}{$c} } ){
				print OUT "$a\t$b\t$c\t$d\t$Altdetails{$a}{$b}{$c}{$d}\t$Details{$a}{$b}{$c}{$d}\n";
			}
		}
    }
}

print "Total number of duplicates $count\n";
exit;

