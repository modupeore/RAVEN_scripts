#!/usr/bin/perl
use strict;
use File::Basename;
use POSIX qw(ceil);
use Sort::Key::Natural qw(natsort);
# merge similar variants based on number of vcfs provided

#- - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
print "\t**MERGE SIMILAR VARIANTS**\n";
my $usage = "
To use : '$0'
\ttype in the AF tables (more than 1) separated by space.\n
";
@ARGV >= 2 or die $usage;

# - - - - - G L O B A L V A R I A B L E S - - - - - - - - -
my (%COMPARE, %header, @headerline, %ALTCOMPARE, %REFCOMPARE, %AFCOMPARE, %AFSCOMPARE, %ANN, %FUNC); 
my ($line, $header);
my @prefix;
my $out = "merge.table.txt";
open (OUT,">", $out);
open (OUT2, ">final-merge.table.txt");
my $libs = 1; #select 1 [1 or more] or 6 [2 or more] or 12  [all libs]
# - - - - - M A I N - - - - - - - - - - - - - - - - - - -

foreach my $no (0..$#ARGV){
	print "Working on $ARGV[$no]\n";
	open (COMMON, $ARGV[$no]) or die "Can't open '$ARGV[$no]'\n$usage";
	my $prefix = fileparse($ARGV[$no], qr/\.[^.]*(\.table.txt)?$/);
	while (<COMMON>){ #reading the input file
		if (/^chr.*/){
			chomp;
			$line = $_;
			my @commonline = split (/\t/, $line);
			#indexing the SNPs in the table provided
			if ($commonline[$header{'AN'}] > $libs) { #to make sure more than 2 libraries for each tissue are represented.
				if (exists $COMPARE{$commonline[$header{"CHROM"}]}{$commonline[$header{"POS"}]}{$prefix}) {
					die "This should not happen\n";
				} else {
					$COMPARE{$commonline[$header{"CHROM"}]}{$commonline[$header{"POS"}]}{$prefix} = "$commonline[$header{'ALT'}]\t$commonline[$header{'FILTER'}]\t$commonline[$header{'AF'}]\t$commonline[$header{'NEW-AF'}]";
				}
				
				#checking the REF allele
				if (exists $REFCOMPARE{$commonline[$header{"CHROM"}]}{$commonline[$header{"POS"}]}) {
					unless ($REFCOMPARE{$commonline[$header{"CHROM"}]}{$commonline[$header{"POS"}]} eq $commonline[$header{'REF'}]) { #making sure alternate alleles called aren't repeated.
						die "Ref allele is different => '$REFCOMPARE{$commonline[$header{'CHROM'}]}{$commonline[$header{'POS'}]}' is not '$commonline[$header{'REF'}]'\n";
					}
				} else {
					$REFCOMPARE{$commonline[$header{'CHROM'}]}{$commonline[$header{'POS'}]} = $commonline[$header{'REF'}];
				}
				
				#getting annotation information
				if (exists $header{'GENE.REFGENE'}) {
					$ANN{$commonline[$header{'CHROM'}]}{$commonline[$header{'POS'}]} = "$commonline[$header{'GENE.REFGENE'}]";
					if ($commonline[$header{'FUNC.REFGENE'}] eq "exonic"){
						$FUNC{$commonline[$header{'CHROM'}]}{$commonline[$header{'POS'}]} = "$commonline[$header{'EXONICFUNC.REFGENE'}]";
					}else {
						$FUNC{$commonline[$header{'CHROM'}]}{$commonline[$header{'POS'}]} = "$commonline[$header{'FUNC.REFGENE'}]";
					}
				}
				
				#summing the AF(s) from GATK & custom only if the filter has PASS, so as not to skew the AF.
				if ($commonline[$header{'FILTER'}] eq "PASS") {
					
					#checking the ALT allele
					my $verdict = "no"; #check if the alt allele is already in the dictionary
					if (exists $ALTCOMPARE{$commonline[$header{"CHROM"}]}{$commonline[$header{"POS"}]}) {
						my @alts = split(",", $ALTCOMPARE{$commonline[$header{"CHROM"}]}{$commonline[$header{"POS"}]});
						foreach (@alts) {
							if ($_ eq $commonline[$header{'ALT'}]) {
								$verdict .= "yes";
							}
						}
						unless ($verdict =~ /yes/) { #making sure alternate alleles called aren't repeated.
							$ALTCOMPARE{$commonline[$header{"CHROM"}]}{$commonline[$header{"POS"}]} = "$ALTCOMPARE{$commonline[$header{'CHROM'}]}{$commonline[$header{'POS'}]},$commonline[$header{'ALT'}]";
						}
					} else {
						$ALTCOMPARE{$commonline[$header{'CHROM'}]}{$commonline[$header{'POS'}]} = $commonline[$header{'ALT'}];
					}
					#summing AF
					unless ($commonline[$header{'AF'}] =~ /\,/) { #this shouldn't happen but just incase
						if (exists $AFCOMPARE{$commonline[$header{"CHROM"}]}{$commonline[$header{"POS"}]}) {
							$AFCOMPARE{$commonline[$header{"CHROM"}]}{$commonline[$header{"POS"}]} = "$AFCOMPARE{$commonline[$header{'CHROM'}]}{$commonline[$header{'POS'}]},$commonline[$header{'AF'}]";
							$AFSCOMPARE{$commonline[$header{"CHROM"}]}{$commonline[$header{"POS"}]} = "$AFSCOMPARE{$commonline[$header{'CHROM'}]}{$commonline[$header{'POS'}]},$commonline[$header{'NEW-AF'}]";
						} else {
							$AFCOMPARE{$commonline[$header{"CHROM"}]}{$commonline[$header{"POS"}]} = $commonline[$header{'AF'}];
							$AFSCOMPARE{$commonline[$header{"CHROM"}]}{$commonline[$header{"POS"}]} = $commonline[$header{'NEW-AF'}];
						}
					} else {
						die "AF has two values.\n";
					}
				} #end if PASS
			} #end if all libraries for each tissue is represented
		} elsif (/^CHROM/) {
			chomp;
			@headerline = split /\t/;
			foreach (0..$#headerline){ $header{uc($headerline[$_])} = $_; }
		} else {
			die "$_ is invalid\n";
		} #end if ^chr
	} #end while
	push @prefix, $prefix;
	close COMMON;
}

#Exporting the table;
#print OUT headers
print OUT "\t\t\t"; #chrom, pos, ref
foreach my $prefix (sort {$a cmp $b} @prefix){
	print OUT "\t\t$prefix\t\t"; #tissues having alt, filter, af, newAF
}
print OUT "\t\t\n"; #mergeALT, avgAD, avgnewAF;
print OUT "CHROM\tPOS\tREF\t";
	
	
foreach my $prefix (sort {$a cmp $b} @prefix){
	print OUT "ALT\tFILTER\t";
	print OUT "AF\tNEW-AF\t";
}
print OUT "mergeALT\t";
if (exists $header{'GENE.REFGENE'}) { print OUT "GENE\tFUNC\t"; }
print OUT "avgAF\tmedAF\tavgnew-AF\tmednew-AF\n";
print OUT2 "CHROM\tPOS\tREF\tmergeALT\t";
if (exists $header{'GENE.REFGENE'}) { print OUT2 "GENE\tFUNC\t"; }
print OUT2 "avgnew-AF\tmednew-AF\n";

#print OUT content
my ($homoalt, $homoref, $hetero) = (0,0,0);
my ($medhomoalt, $medhomoref, $medhetero) = (0,0,0);

open (OUT3, ">homo-alt-average-freq.txt");
open (OUT4, ">homo-ref-average-freq.txt");
open (OUT5, ">hetero-average-freq.txt");
open (OUT6, ">homo-alt-median-freq.txt");
open (OUT7, ">homo-ref-median-freq.txt");
open (OUT8, ">hetero-median-freq.txt");
print OUT3 "CHROM\tPOS\tREF\tmergeALT\t"; if (exists $header{'GENE.REFGENE'}) { print OUT3 "GENE\tFUNC\t"; } print OUT3 "avgnew-AF\n";
print OUT4 "CHROM\tPOS\tREF\tmergeALT\t"; if (exists $header{'GENE.REFGENE'}) { print OUT4 "GENE\tFUNC\t"; } print OUT4 "avgnew-AF\n";
print OUT5 "CHROM\tPOS\tREF\tmergeALT\t"; if (exists $header{'GENE.REFGENE'}) { print OUT5 "GENE\tFUNC\t"; } print OUT5 "avgnew-AF\n";
print OUT6 "CHROM\tPOS\tREF\tmergeALT\t"; if (exists $header{'GENE.REFGENE'}) { print OUT6 "GENE\tFUNC\t"; } print OUT6 "avgnew-AF\n";
print OUT7 "CHROM\tPOS\tREF\tmergeALT\t"; if (exists $header{'GENE.REFGENE'}) { print OUT7 "GENE\tFUNC\t"; } print OUT7 "avgnew-AF\n";
print OUT8 "CHROM\tPOS\tREF\tmergeALT\t"; if (exists $header{'GENE.REFGENE'}) { print OUT8 "GENE\tFUNC\t"; } print OUT8 "avgnew-AF\n";

foreach my $chrom (natsort keys %COMPARE) {
	foreach my $position (sort {$a <=> $b} keys %{ $COMPARE{$chrom} } ) {
		if (exists $ALTCOMPARE{$chrom}{$position}) { #if filter(s) has at least PASS
			print OUT "$chrom\t$position\t$REFCOMPARE{$chrom}{$position}\t";
			print OUT2 "$chrom\t$position\t$REFCOMPARE{$chrom}{$position}\t";
			foreach my $prefix (sort {$a cmp $b} @prefix ){
				if (exists $COMPARE{$chrom}{$position}{$prefix}){
					print OUT "$COMPARE{$chrom}{$position}{$prefix}\t";
				} else {
					print OUT "-\t-\t-\t-\t";
				}
			}
			print OUT "$ALTCOMPARE{$chrom}{$position}\t";
			print OUT2 "$ALTCOMPARE{$chrom}{$position}\t";
			
			if (exists $header{'GENE.REFGENE'}){
				print OUT "$ANN{$chrom}{$position}\t$FUNC{$chrom}{$position}\t";
				print OUT2 "$ANN{$chrom}{$position}\t$FUNC{$chrom}{$position}\t";
			}

			my @total = split(",",$AFCOMPARE{$chrom}{$position}); #print "$AFCOMPARE{$chrom}{$position} , this is total@total\n";
			my $average = sprintf("%.3f", &average(\@total));
			my $median = sprintf("%.3f", &median(\@total)); #print "here $median\n";
			@total = split(",",$AFSCOMPARE{$chrom}{$position}); #print "this is total@total\n";
			my $afaverage = sprintf("%.3f", &average(\@total));
			my $afmedian = sprintf("%.3f", &median(\@total)); #print "there $median\n"; 
			#my $average = sprintf("%.3f", ($AFCOMPARE{$chrom}{$position} / $AFcount{$chrom}{$position}));
			#my $afaverage = sprintf("%.3f", ($AFSCOMPARE{$chrom}{$position} / $AFcount{$chrom}{$position}));
			print OUT "$average\t$median\t$afaverage\t$afmedian\n";
			print OUT2 "$afaverage\t$afmedian\n";
			
			#Frequencies of AFs
			if ($afaverage > 0.99) {	 #homo alt
				print OUT3 "$chrom\t$position\t$REFCOMPARE{$chrom}{$position}\t$ALTCOMPARE{$chrom}{$position}\t";
				if (exists $header{'GENE.REFGENE'}){ print OUT3 "$ANN{$chrom}{$position}\t$FUNC{$chrom}{$position}\t"; }
				print OUT3 "$afaverage\n";
				$homoalt++;
			} elsif ($afaverage < 0.01) { #homo ref
				print OUT4 "$chrom\t$position\t$REFCOMPARE{$chrom}{$position}\t$ALTCOMPARE{$chrom}{$position}\t";
				if (exists $header{'GENE.REFGENE'}){ print OUT4 "$ANN{$chrom}{$position}\t$FUNC{$chrom}{$position}\t"; }
				print OUT4 "$afaverage\n";
				$homoref++;
			} else {
				print OUT5 "$chrom\t$position\t$REFCOMPARE{$chrom}{$position}\t$ALTCOMPARE{$chrom}{$position}\t";
				if (exists $header{'GENE.REFGENE'}){ print OUT5 "$ANN{$chrom}{$position}\t$FUNC{$chrom}{$position}\t"; }
				print OUT5 "$afaverage\n";
				$hetero++;
			}
			
			if ($afmedian > 0.99) { #homo alt
				print OUT6 "$chrom\t$position\t$REFCOMPARE{$chrom}{$position}\t$ALTCOMPARE{$chrom}{$position}\t";
				if (exists $header{'GENE.REFGENE'}){ print OUT6 "$ANN{$chrom}{$position}\t$FUNC{$chrom}{$position}\t"; }
				print OUT6 "$afaverage\n";
				$medhomoalt++;
			} elsif ($afmedian < 0.01) { #homo ref
				print OUT7 "$chrom\t$position\t$REFCOMPARE{$chrom}{$position}\t$ALTCOMPARE{$chrom}{$position}\t";
				if (exists $header{'GENE.REFGENE'}){ print OUT7 "$ANN{$chrom}{$position}\t$FUNC{$chrom}{$position}\t"; }
				print OUT7 "$afaverage\n";
				$medhomoref++;
			} else {
				print OUT8 "$chrom\t$position\t$REFCOMPARE{$chrom}{$position}\t$ALTCOMPARE{$chrom}{$position}\t";
				if (exists $header{'GENE.REFGENE'}){ print OUT8 "$ANN{$chrom}{$position}\t$FUNC{$chrom}{$position}\t"; }
				print OUT8 "$afaverage\n";
				$medhetero++;
			}
		} 
	}
}

close OUT;
close OUT2;
close OUT3;
close OUT4;
close OUT5;
close OUT6;
close OUT7;
close OUT8;
print "Finished\n";
print "
criteria 0.01;
Average
\t0/0 = $homoref
\t0/1 = $hetero
\t1/1 = $homoalt

Median
\t0/0 = $medhomoref
\t0/1 = $medhetero
\t1/1 = $medhomoalt
";

sub range {
			my ($data) = @_;
			my $max = (sort {$b <=> $a} @$data)[0];
			my $min = (sort {$a <=> $b} @$data)[0];
			return $max, $min;
}
sub median {
				my $median; my ($data) = @_;
				my $length = @$data;
				if ($length > 1) {
				my @median = (sort {$a <=> $b} @$data)[ int($length/2), ceil($length /2) ];
				$median = sprintf("%.3f", (&sum(\@median)/2));
				} else {
					$median = sprintf("%.3f",@$data[0]);
				}
				return $median;
}
sub sum {
				my ($data) = @_;
				if (not @$data) {
                die("Empty arrayn");
        }
				my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
				return $total;
}
sub average{
        my ($data) = @_;
        if (not @$data) {
                die("Empty arrayn");
        }
        my $total = &sum($data);
        my $average = $total / @$data;
        return $average;
}
sub stdev{
        my($data) = @_;
        if(@$data == 1){
                return 0;
        }
        my $average = &average($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}
