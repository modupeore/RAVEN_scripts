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
my (@unexonicarray,@unnonsyn, @unsyn, @unstop, @unrest, @unframeshift, @unnonframeshift);
my (@nonsyn, @syn, @stop, @rest, @frameshift, @nonframeshift);
my (%NONSYN, %SYN, %STOP, %FRAME, %NONFRAME);

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
			if ($odachrIwant =~ /^exonic/){
				push (@unexonicarray, $newchr);
				my $newchrIwant = $chrdetails[13];
				if($newchrIwant =~ /^nonsy/){
					push (@unnonsyn, $newchr);
				}
				elsif ($newchrIwant =~ /^syn/) {
					push (@unsyn, $newchr);
				}
				elsif ($newchrIwant =~ /^stop/) {
					push (@unstop, $newchr);
				}
				elsif ($newchrIwant =~ /^frame/) {
					push (@unframeshift, $newchr);
				}
				elsif ($newchrIwant =~ /^nonframe/) {
					push (@unnonframeshift, $newchr);
				}
				else {
					push (@unrest, $newchr);
					print "Unrest == > $newchr\n";
				}	
			}
			
		} else {
			push (@nonunknownarray, $newchr);
			my @chrdetails = split('\t', $newchr);
			my $newchrIwant = $chrdetails[8]; 
			if($newchrIwant =~ /^nonsy/){
				push (@nonsyn, $newchr);
			}
			elsif ($newchrIwant =~ /^syn/) {
				push (@syn, $newchr);
			}
			elsif ($newchrIwant =~ /^stop/) {
				push (@stop, $newchr);
			}
			elsif ($newchrIwant =~ /^frame/) {
				push (@frameshift, $newchr);
			}
			elsif ($newchrIwant =~ /^nonfra/) {
				push (@nonframeshift, $newchr); 
			}
			else {
				push (@rest, $newchr);
				print "Rest == > $newchr\n";
			}
		}
	}
}
#nonsyn
foreach my $spec (@nonsyn){
	unless ($spec =~ /^\s/){
		my @chrdetails = split('\t', $spec);
		my $thechr = $chrdetails[6];
		if($thechr =~ /\,/) {
			my @newthechr = split("\,",$thechr);
			foreach (@newthechr) {
				$NONSYN{$_} ++;
			}
		} else {
			$NONSYN{$thechr} ++;
		}
	}
}
foreach my $spec (@unnonsyn){
	unless ($spec =~ /^\s/){
		my @chrdetails = split('\t', $spec);
		my $thechr = $chrdetails[11]; 
		if($thechr =~ /\,/) {
			my @newthechr = split("\,",$thechr);
			foreach (@newthechr) {
				my @notagain = split("\\.", $_);
				$NONSYN{$notagain[0]} ++;
			}
		} else {
			my @notagain = split("\\.", $thechr); 
			$NONSYN{$notagain[0]} ++;
		}
	}
}
#syn
foreach my $spec (@syn){
	unless ($spec =~ /^\s/){
		my @chrdetails = split('\t', $spec);
		my $thechr = $chrdetails[6];
		if($thechr =~ /\,/) {
			my @newthechr = split("\,",$thechr);
			foreach (@newthechr) {
				$SYN{$_} ++;
			}
		} else {
			$SYN{$thechr} ++;
		}
	}
}
foreach my $spec (@unsyn){
	unless ($spec =~ /^\s/){
		my @chrdetails = split('\t', $spec);
		my $thechr = $chrdetails[11];
		if($thechr =~ /\,/) {
			my @newthechr = split("\,",$thechr);
			foreach (@newthechr) {
				my @notagain = split("\\.", $_);
				$SYN{$notagain[0]} ++;
			}
		} else {
			my @notagain = split("\\.", $thechr);
			$SYN{$notagain[0]} ++;
		}
	}
}
#stop
foreach my $spec (@stop){
	unless ($spec =~ /^\s/){
		my @chrdetails = split('\t', $spec);
		my $thechr = $chrdetails[6];
		if($thechr =~ /\,/) {
			my @newthechr = split("\,",$thechr);
			foreach (@newthechr) {
				$STOP{$_} ++;
			}
		} else {
			$STOP{$thechr} ++;
		}
	}
}
foreach my $spec (@unstop){
	unless ($spec =~ /^\s/){
		my @chrdetails = split('\t', $spec);
		my $thechr = $chrdetails[11];
		if($thechr =~ /\,/) {
			my @newthechr = split("\,",$thechr);
			foreach (@newthechr) {
				my @notagain = split("\\.", $_);
				$STOP{$notagain[0]} ++;
			}
		} else {
			my @notagain = split("\\.", $thechr);
			$STOP{$notagain[0]} ++;
		}
	}
}
#frameshift
foreach my $spec (@frameshift){
	unless ($spec =~ /^\s/){
		my @chrdetails = split('\t', $spec);
		my $thechr = $chrdetails[6];
		if($thechr =~ /\,/) {
			my @newthechr = split("\,",$thechr);
			foreach (@newthechr) {
				$FRAME{$_} ++;
			}
		} else {
			$FRAME{$thechr} ++;
		}
	}
}
foreach my $spec (@unframeshift){
	unless ($spec =~ /^\s/){
		my @chrdetails = split('\t', $spec);
		my $thechr = $chrdetails[11];
		if($thechr =~ /\,/) {
			my @newthechr = split("\,",$thechr);
			foreach (@newthechr) {
				my @notagain = split("\\.", $_);
				$FRAME{$notagain[0]} ++;
			}
		} else {
			my @notagain = split("\\.", $thechr);
			$FRAME{$notagain[0]} ++;
		}
	}
}
#nonframeshift
foreach my $spec (@nonframeshift){
	unless ($spec =~ /^\s/){
		my @chrdetails = split('\t', $spec);
		my $thechr = $chrdetails[6];
		if($thechr =~ /\,/) {
			my @newthechr = split("\,",$thechr);
			foreach (@newthechr) {
				$NONFRAME{$_} ++;
			}
		} else {
			$NONFRAME{$thechr} ++;
		}
	}
}
foreach my $spec (@unnonframeshift){
	unless ($spec =~ /^\s/){
		my @chrdetails = split('\t', $spec);
		my $thechr = $chrdetails[11];
		if($thechr =~ /\,/) {
			my @newthechr = split("\,",$thechr);
			foreach (@newthechr) {
				my @notagain = split("\\.", $_);
				$NONFRAME{$notagain[0]} ++;
			}
		} else {
			my @notagain = split("\\.", $thechr);
			$NONFRAME{$notagain[0]} ++;
		}
	}
}

open(NON,">", $out."_nonsyn.genes.txt");
open(SYN,">", $out."_syn.genes.txt");
open(STOP,">", $out."_stop.genes.txt");
open(FRAME,">", $out."_frameshift.genes.txt");
open(NONFRAME,">", $out."_nonframeshift.genes.txt");

foreach my $newkey (sort keys %NONSYN){
	print NON "$newkey\n";
}

foreach my $newkey (sort keys %SYN){
	print SYN "$newkey\n";
}

foreach my $newkey (sort keys %STOP){
	print STOP "$newkey\n";
}
foreach my $newkey (sort keys %FRAME){
	print FRAME "$newkey\n";
}
foreach my $newkey (sort keys %NONFRAME){
	print NONFRAME "$newkey\n";
}
#
#print "\n$count\n";
#exit;
