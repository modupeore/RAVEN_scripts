#!/usr/bin/perl
use strict;
use POSIX qw(ceil);
use File::Basename;

#trying to identify the mean of DPs
my $usage = "$0 <vcf file>\n";
@ARGV==1 || die "Forgot to include vcf file\n$usage";
my ($ind, $i) = (7, 0); my $aaa;
my @total;
my (%VALUES);
open(FILE,$ARGV[0]) or die "File '$ARGV[0]' doesn't exist\n"; #input file
undef @total;
while (<FILE>) {
	chomp;
	my $file = $_;
	unless ($file =~ /^#/) {
		$i++; #if ($i == 10){ last; }
		my @chrdetails = split(/\t/,$file);
		$ind = 7; if ($chrdetails[$ind] =~ /DP=([0-9]*)\;/) { $aaa = $1;} 
		$ind = 7; if ($chrdetails[$ind] =~ /QD=([0-9]*\.[0-9]*)/) { push @total, $1; print $aaa,"\t", $1,"\n";$VALUES{$1}{$i} = $file;} #Mapping Quality
		#$ind = 7; if ($chrdetails[$ind] =~ /AN=([0-9]*)\;/) { push @total, $1;  $VALUES{$1}{$i} = $file;} #Read depth
		#$ind = 9; my @d = split('\:', $chrdetails[$ind]); push @total, $d[2]; #Read depth
	}
}
close (FILE);
my $ave = &average(\@total);
my $std = &stdev(\@total);
my $sum = &sum(\@total);
my $median = &median(\@total);
my ($max, $min) = &range(\@total);
#print join "\n" , @total;
my $above = $ave + (3 * $std);
my $below = $ave - (3 * $std);
print "\nMedian : $median,
Sum : $sum,
Average : $ave,
Stdev : $std,
+3SD : $above,
-3SD : $below,
Max : $max,
Min : $min\n"; die;
open (OUT, ">output.vcf");
foreach my $key (keys %VALUES){
	if ($key >= $below && $key <= $above) {
		foreach my $odakey (keys %{$VALUES{$key}}){
			print OUT $VALUES{$key}{$odakey},"\n";
		}
	}
}



sub range {
			my ($data) = @_;
			my $max = (sort {$b <=> $a} @$data)[0];
			my $min = (sort {$a <=> $b} @$data)[0];
			return $max, $min;
}
sub median {
				my ($data) = @_;
				my $length = @$data;
				my @median = (sort {$a <=> $b} @$data)[ int($length/2), ceil($length /2) ];
				my $median = &sum(\@median)/2;
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
