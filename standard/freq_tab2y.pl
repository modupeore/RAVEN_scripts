#!/usr/bin/perl
use strict;
use File::Basename;

#OBJECTIVE
my $usage="$0
\t<Tab-delimited file>
\t<Number of column to be sorted>
";

@ARGV==2 or die $usage;
my $filen = $ARGV[0];
my $column = $ARGV[1]-1;
print "\n\tFrequency of annotations called from Tabulated file\n";
my $input = $filen;
my $i= 0;
unless(open(FILE,$input)){
	print "File \'$input\' doesn't exist\n";
	exit;
}

my @file = <FILE>;
chomp @file;
close (FILE);
my $count = 0;
my %Hashdetails;
foreach my $chr (@file){
	unless ($chr =~ /^\s/){
		my @chrdetails = split('\t', $chr);
		my $chrIwant = "$chrdetails[0]-$chrdetails[1]";
		my @morechrsplit = split('\(', $chrIwant);
		$Hashdetails{$morechrsplit[0]} = 0;
	}
}

foreach my $newcount (@file){
        unless ($newcount =~ /^\s/){
                my @details = split('\t', $newcount);
                my $chragain = "$details[0]-$details[1]";
                my @morechr = split('\(', $chragain);
                $Hashdetails{$morechr[0]} ++; $count++;
		#print $count++."\t$details[0]\t$details[1]\t\t$morechr[0]\t$Hashdetails{$morechr[0]}\n";
        }
}

print "ABC\n";
foreach my $newkey (sort keys %Hashdetails){
	print "$newkey\t=$Hashdetails{$newkey}\n";
}

print "\n$count\n";
exit;
