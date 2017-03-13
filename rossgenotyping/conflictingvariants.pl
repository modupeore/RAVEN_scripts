#!/usr/bin/perl
use strict;
# identify conflicting variants

#- - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
print "\t**Identify conflicting variants**\n";
my $usage = "
To use : '$0'
\ttype in the merged-vcf file to find conflicting variants.
\ttype in the output file name.\n
";
@ARGV == 2 or die $usage;

# - - - - - G L O B A L V A R I A B L E S - - - - - - - - -

my $i = 0;
my ($line, $header);
my $out = $ARGV[1];
open (OUT,">", $out);
# - - - - - M A I N - - - - - - - - - - - - - - - - - - -

open (COMMON, $ARGV[0]) or die "Can't open '$ARGV[0]'\n$usage";
while (<COMMON>){
	if (/^chr.*/){
		$i++;
		$line = $_;
    my @commonline = split (/\t/, $line);
		if ($commonline[4] =~ /\,/){
			if ($commonline[6] =~ /PASS/){
				$commonline[6] = "CONFLICT";
			} else {
				$commonline[6] .= ";CONFLICT";
			}
		}
		print OUT $commonline[0]; foreach (1..$#commonline-1) { print OUT "\t", $commonline[$_]; } print OUT "\t", $commonline[$#commonline];
	} else {
		print OUT $_;
	}
}
close COMMON;
