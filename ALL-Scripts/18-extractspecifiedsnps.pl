#!/usr/bin/perl
use strict;
# to extract the not listed SNPs
# 


# - - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
$a = "\t**SELECTING THE SNPS THAT FIT THE SPECIFIED REGION**\n\n";
my $b = "Type in the name\n\t\t1. The \"SNP file\".\n\t\t2. The \"chromosome\".\n\t\t3. The \"starting location\".\n\t\t4. The \"ending location\".\n";

my $snpfile = $ARGV[0];
my $chr = $ARGV[1];
my $first = $ARGV[2];
my $second = $ARGV[3];

# - - - - - U S E R    V A R I A B L E S  - - - - - - - - - -
# open SNP file
open (FILE, "<$snpfile") or die "$a$b Cannot find file $snpfile\n";
my @snps = <FILE>; close FILE;

# the output file
my $output = $chr."_".$first."-".$second.".vcf";
open (OUTPUTFILE, ">$output");

# - - - - - G L O B A L  V A R I A B L E S  - - - - - - - - -
my %Hashtable;
my $i = 0;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#1. getting the list of all the snps
foreach my $sepsnps (@snps) {
	if ($sepsnps =~ m/^$chr\t/ || $sepsnps != m/^\#.*/){
		my @newsepsnps = split('\t',$sepsnps);
			$Hashtable{$newsepsnps[1]} = $sepsnps;
	} elsif ($sepsnps =~ m/^\#.*/) {
		print OUTPUTFILE $sepsnps;
	}
}

# 2. extracting snps of my specified region
foreach my $tocheck (keys %Hashtable) {
	if ($tocheck >= $first && $tocheck <= $second){
		print OUTPUTFILE $Hashtable{$tocheck};
            	$i++;
	}
}
print "Total number of read matching $chr:$first-$second is \"$i\".\n";
print "Successfully saved in \"$output\"\n\n";
close OUTPUTFILE;
print "\n\n*****************DONE*****************\n\n";
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - S U B R O U T I N E S - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -



