#!/usr/bin/perl
use strict;
use File::Basename;

#OBJECTIVE
my $usage="$0
\t<Tab-delimited file>
";

@ARGV==1 or die $usage;
my $filen = $ARGV[0];
print "\n\tGetting the genes\n";
my $input = $filen;
my $i= 0;
unless(open(FILE,$input)){
	print "File \'$input\' doesn't exist\n";
	exit;
}
my $out = fileparse($input, qr/(\.txt)?$/);

my @file = <FILE>;
chomp @file;
close (FILE);
my $count = 0;
my (@exonicarray, @unknownarray, @nonunknownarray) ;
my (@unexonicarray,@unnonsyn, @unsyn, @unstop, @unrest);
my (@nonsyn, @syn, @stop, @rest);
my (%NONSYN, %SYN, %STOP);

#get exonic
foreach my $chr (@file){
	unless ($chr =~ /^\s/){
		my @chrdetails = split('\t', $chr);
		my $chrIwant = $chrdetails[5];
		if ($chrIwant =~ /^exonic/){
			push (@exonicarray, $chr);
		}
	}
}
#get unknown from exonic
foreach my $newchr (@exonicarray){
	unless ($newchr =~ /^\s/){
		my @chrdetails = split('\t', $newchr);
		my $chrIwant = $chrdetails[8];
		if ($chrIwant =~ /unknown/){
			push (@unknownarray, $newchr);
			my $odachrIwant = $chrdetails[10];
			if ($odachrIwant =~ /exonic/){
				push (@unexonicarray, $newchr);
				my $newchrIwant = $chrdetails[13];
				if($newchrIwant =~ /^non/){
					push (@unnonsyn, $newchr);
				}
				elsif ($newchrIwant =~ /^syn/) {
					push (@unsyn, $newchr);
				}
				elsif ($newchrIwant =~ /^stop/) {
					push (@unstop, $newchr);
				}
				else {
					push (@unrest, $newchr);
				}	
			}
			
		} else {
			push (@nonunknownarray, $newchr);
			my @chrdetails = split('\t', $newchr);
			my $newchrIwant = $chrdetails[8];
			if($newchrIwant =~ /^non/){
				push (@nonsyn, $newchr);
			}
			elsif ($newchrIwant =~ /^syn/) {
				push (@syn, $newchr);
			}
			elsif ($newchrIwant =~ /^stop/) {
				push (@stop, $newchr);
			}
			else {
				push (@rest, $newchr);
			}
		}
	}
}
#nonsyn
foreach my $spec (@nonsyn){
	unless ($spec =~ /^\s/){
		my @chrdetails = split('\t', $spec);
		#print $chrdetails[5];
		$NONSYN{$chrdetails[6]} ++;
	}
}
foreach my $spec (@unnonsyn){
	unless ($spec =~ /^\s/){
		my @chrdetails = split('\t', $spec);
		#print $chrdetails[5];
		$NONSYN{$chrdetails[11]} ++;
	}
}
#syn
foreach my $spec (@syn){
	unless ($spec =~ /^\s/){
		my @chrdetails = split('\t', $spec);
		#print $chrdetails[5];
		$SYN{$chrdetails[6]} ++;
	}
}
foreach my $spec (@unsyn){
	unless ($spec =~ /^\s/){
		my @chrdetails = split('\t', $spec);
		#print $chrdetails[5];
		$SYN{$chrdetails[11]} ++;
	}
}
#stop
foreach my $spec (@stop){
	unless ($spec =~ /^\s/){
		my @chrdetails = split('\t', $spec);
		#print $chrdetails[5];
		$STOP{$chrdetails[6]} ++;
	}
}
foreach my $spec (@unstop){
	unless ($spec =~ /^\s/){
		my @chrdetails = split('\t', $spec);
		#print $chrdetails[5];
		$STOP{$chrdetails[11]} ++;
	}
}

open(NON,">", $out."_nonsyn.genes");
open(SYN,">", $out."_syn.genes");
open(STOP,">", $out."_stop.genes");

foreach my $newkey (sort keys %NONSYN){
	print NON "$newkey\t$NONSYN{$newkey}\n";
}

foreach my $newkey (sort keys %SYN){
	print SYN "$newkey\t$SYN{$newkey}\n";
}

foreach my $newkey (sort keys %STOP){
	print STOP "$newkey\t$STOP{$newkey}\n";
}
#
#print "\n$count\n";
#exit;
