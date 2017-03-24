#!/usr/bin/perl
use strict;
use Pod::Usage;
use Getopt::Long;
use List::MoreUtils 'pairwise';
# convert to table and estimate AF

#- - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
print "\t**Estimate AF and convert to table based on columns specified excluding the FORMAT fields**\n";
my $usage = "
To use : '$0'
\t-V|v VCF file.
\t-F|f columns (multiple columns are separated by comma).
\t-o|O output file name.
";
my ($variant, $columns, $output);
GetOptions('V|v=s'=>\$variant, 'F|f=s'=>\$columns , 'O|o=s'=>\$output) or die $usage;
die "Incomplete Arguments!\n$usage" if(!$variant || !$columns || !$output);

# - - - - - G L O B A L V A R I A B L E S - - - - - - - - -

my (%header, %info, %format, %filter, @headerline, $line, $verdict, @sumAC);
open (OUT,">", $output);

# - - - - - M A I N - - - - - - - - - - - - - - - - - - -

#print out header
print OUT "CHROM\tPOS\t";

#Parse columns specified
$columns =~ s/\s+|\s+//g;
my @filter = split (",", $columns);
my $i = 0; foreach (@filter){ $filter{$i++} = uc($_); unless ( ($_ =~ /CHROM/i) || ($_ =~ /POS/i)){ print OUT uc($_),"\t"; } }
print OUT "new-AC\tnew-AF\n";
#variant file
open (COMMON, $variant) or die "Can't open '$variant'\n$usage";
while (<COMMON>){
	chomp;
	if (/^chr.*/) {
		undef $verdict;
		my ($newAC, $newAF) = (0,0);
		my (@newAC);
		$line = $_;
    my @commonline = split (/\t/, $line);
		
		#extracting all the columns into a hash key based on the info column.
		my $selectfilter = $commonline[$header{"INFO"}];
		my @info = split (/\;/, $selectfilter);
		undef %info; #clean the info hash
		foreach (@info) {
			my @theinfo = split("\=", $_, 2);
			$info{uc($theinfo[0])} = $theinfo[1];
		}
		
		#estimating AF from FORMAT column
		my $selectformat = $commonline[$header{"FORMAT"}];
		my @format = split (/\:/, $selectformat);
		undef %format; #clean format hash
		foreach (0..$#format){ $format{$format[$_]} = $_; }
		if (exists $format{"AD"}) {
			my $cols = $header{"FORMAT"} + 1;
			foreach my $colloc ($cols..$#headerline) {
				if ($commonline[$colloc] =~ /^\d/) {
					my $allAD = (split (/\:/,$commonline[$colloc]))[$format{"AD"}];
					if ($allAD =~ /^\d/) {
						my ($ADref, $ADalt) = split(",", $allAD,2); #prefect diploid cases
						if ($ADalt =~ /,/){
							#print $line,"\n";
							my @tmpAC = @newAC;
							my $sumAD = $ADref;
							$verdict = "yes";
							my @alts = split(/,/, $ADalt);
							foreach (@alts) { $sumAD +=$_; }
							
							#print "$ADalt\t@alts\t$sumAD\t@tmpAC\n";
							my @sumAC = map { sprintf("%.3f", ($_ / $sumAD)) } @alts;
							@newAC = pairwise { $a + $b } @tmpAC, @sumAC;
						} else {
							$verdict = "no";
							#print "$commonline[$colloc]\t$ADalt\t$ADref\n";
							if ($ADalt == 0 && $ADref == 0){ $newAC += 0 ; }
							else { $newAC += (sprintf ("%.3f", ($ADalt/($ADref+$ADalt)))); }
						}
					#} else {
					#	if($commonline[$header{"ALT"}] =~ /,/) {
					#		#die;
					#		my @tmpAC = @newAC;
					#		$verdict = "yes";
					#		my @alts = split(/,/, $commonline[$header{"ALT"}]);
					#		my $sumAC; foreach (@alts){push @sumAC, 1;}
					#		@newAC = pairwise { $a + $b } @tmpAC, @sumAC;
					#	} else {
					#		$verdict = "no";
					#		$newAC += 1;
					#	}
					}
				}
			} 
			if ($verdict eq "no") {
				#print "ADalt $sumADalt\tAD ",$sumAD,"\t","AC",$newAC,"\tAN ", $info{"AN"},"\n";
				$newAF = (sprintf ("%.3f", ($newAC/($info{"AN"}/2)))); #convert it to a haploid case, to estimate AF
			}elsif ($verdict eq "yes") {
				my @newAF = map { sprintf("%.3f", ($_ / ($info{"AN"}/2))) } @newAC;
				$newAC = join(",", @newAC);
				$newAF = join(",", @newAF);
			}
		} else {
			$newAC = "null";
			$newAF = "null";
		}
		#print out columns specified and the new file
		
		print OUT $commonline[$header{"#CHROM"}],"\t", $commonline[$header{"POS"}];
		foreach (sort {$a <=> $b} keys %filter){
			unless ( ($filter{$_} =~ /CHROM/) || ($filter{$_} =~ /POS/)){
				if (exists $header{$filter{$_}}){
					print OUT "\t$commonline[$header{$filter{$_}}]";
				} elsif (exists $info{$filter{$_}}){
					print OUT "\t$info{$filter{$_}}";
				} else { die "'$filter{$_}' does not exist in the variant file\n"; }
			}
		}
		print OUT "\t$newAC\t$newAF\n";
	} elsif (/^#CH/) { #get the header column
		@headerline = split /\t/;
		foreach (0..$#headerline){ $header{$headerline[$_]} = $_; }	
	}
}
close COMMON;
