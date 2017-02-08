#!/usr/bin/perl
use strict;
use IO::File;
use Pod::Usage;
use Getopt::Long;
use Sort::Key::Natural qw(natsort);
use Statistics::R;

my $usage="Arguments needed for '$0'
\t-variant <Variant file (multiple separated by comma>
\t-fasta <Reference fasta file>
\t-bin <Sliding window size>
\t[-chr <chromosome name (multiple separated by comma)>]
\t[-output <Output file name>]
";

my (%CHR, %GENOME);
my($variantfile, $fastafile,$outputfile,$bin, $chrom, %chrom, @variantfile, $new);
GetOptions('variant|v=s'=>\$variantfile, 'fasta|f=s'=>\$fastafile,
		'output|o=s'=>\$outputfile, 'bin|b=s'=>\$bin, 'chr|chrom|c=s'=>\$chrom ) or die $usage;

die $usage unless ($variantfile && $fastafile && $bin);

$outputfile ||="output.txt"; $outputfile = @{ open_unique($outputfile) }[1]; $outputfile =~ /^(.*)\..+$/;
my $picture = $1.".png"; $picture = @{ open_unique("$picture") }[1]; $picture =~ /^(.*)\.png$/; my $index = $1;
@ARGV == 0 || die $usage;

GENOME(); #read genome file

if ($chrom){
	print "chromosomes selected $chrom\n";
	$chrom =~ s/\s+|\s+//g;
	%chrom = map {$_ => 1} (split (",", $chrom));
} else {
	print "Processing All chromosomes\n";
}
if ($variantfile){
	print "Variant file selected $variantfile\n";
	@variantfile = split (",", $variantfile);
}
foreach my $files (@variantfile) {
	$files =~ s/^\s+|\s+$//g;
	open(IN,"<", $files);
	while (<IN>){
		unless (/^#/){
			my @all = split /\t/;
			unless ($all[0] =~ /random/){
				if ($chrom){
					if (exists $chrom{$all[0]}) {
						my $verdict = int($all[1]/$bin);
						unless ($verdict == 0) { $verdict = $verdict * $bin; } else {$verdict = 1;}
						$CHR{$all[0]}{$verdict}{$files} = $CHR{$all[0]}{$verdict}{$files}+1;
					}
				} else {
					my $verdict = int($all[1]/$bin);
					unless ($verdict == 0) { $verdict = $verdict * $bin; } else {$verdict = 1;}
					$CHR{$all[0]}{$verdict}{$files} = $CHR{$all[0]}{$verdict}{$files}+1;
				}
			}
		}
	}
}
open (OUT, ">$outputfile"); #save to outputfile
print OUT "chrom\tposition";
unless ($#variantfile >= 1){
	print OUT "\tcount\n";
} else { 
	foreach my $files (@variantfile){ $files =~ s/^\s+|\s$//g; print OUT "\t$files"; }
	print OUT "\n";
}
foreach my $a (natsort keys %CHR) {
	for (my $x = 0; $x <= $GENOME{$a}; $x=$x+$bin){
		my $v;if ($x ==0) { $v = 1; } else { $v = $x; }
		print OUT "$a\t$v";
		foreach my $files (@variantfile){
			$files =~ s/^\s+|\s$//g;
			if (exists $CHR{$a}{$v}{$files}){ print OUT "\t",$CHR{$a}{$v}{$files}; } else { print OUT "\t0"; }
		}
		print OUT "\n";
	}
}
close (OUT);
Rscript($outputfile, $picture);

my $tempout = @{ open_unique(".tempout.txt")}[1]; #ploting 10 chrs in one pdf
my ($count, $counter) = (1,1);
my $pictures= $index."-".$count.".png";
open (TEMPOUT, ">$tempout");
print TEMPOUT "chrom\tposition";
unless ($#variantfile >= 1){
	print TEMPOUT "\tcount\n";
} else { 
	foreach my $files (@variantfile){ $files =~ s/^\s+|\s$//g; print TEMPOUT "\t$files"; }
	print TEMPOUT "\n";
}
foreach my $a (natsort keys %CHR) {
	unless ($a eq $new) {
		if ($counter > 10) {
			$pictures= $index."-".$count.".png";
			$count++;
			Rscript($tempout, $pictures);
			open (TEMPOUT, ">$tempout");
			print TEMPOUT "chrom\tposition";
			unless ($#variantfile >= 1){
				print TEMPOUT "\tcount\n";
			} else { 
			foreach my $files (@variantfile){ $files =~ s/^\s+|\s$//g; print TEMPOUT "\t$files"; }
			print TEMPOUT "\n";
			}
			$counter = 1;
		}
		else {
			$counter++;
		}
		$new = $a;
	}
	for (my $x = 0; $x <= $GENOME{$a}; $x=$x+$bin){
		my $v;if ($x ==0) { $v = 1; } else { $v = $x; }
		print TEMPOUT "$a\t$v";
		foreach my $files (@variantfile){
			$files =~ s/^\s+|\s$//g;
			if (exists $CHR{$a}{$v}{$files}){ print TEMPOUT "\t",$CHR{$a}{$v}{$files}; } else { print TEMPOUT "\t0"; }
		}
		print TEMPOUT "\n";
	}
}
close (TEMPOUT);
$pictures= $index."-".$count.".png";
Rscript($tempout,$pictures);
`rm -rf $tempout`;


##SUBROUTINES
sub GENOME{
	$/ = ">";
	open (REF, "<", $fastafile);
	my @genomefile = <REF>; close REF;
	shift @genomefile;
	foreach (@genomefile){
		my @pieces = split /\n/;
		my $total = 0;
		foreach my $num (1..$#pieces){
			$total = $total + length($num);
		}
		my $newtotal = int($total/$bin)*$bin; if ($newtotal < $total){ $newtotal= $newtotal+$bin;}
		$GENOME{$pieces[0]} = $newtotal;
	}
	$/ = "\n";
}

sub open_unique {
	  my $file = shift || '';
    unless ($file =~ /^(.*?)(\.[^\.]+)$/) {
        print "Bad file name: '$file'\n";
        return;
    }
    my $io;
    my $seq = '';
    my $base = $1;
    my $ext = $2;
    until (defined ($io = IO::File->new($base.$seq.$ext
                                  ,O_WRONLY|O_CREAT|O_EXCL))) {
        last unless $!{EEXIST};
        $seq = '_0' if $seq eq '';
        $seq =~ s/(\d+)/$1 + 1/e;
    }
    return [$io,$base.$seq.$ext] if defined $io;
		#http://stackoverflow.com/questions/5188806/perl-open-file-but-not-overwrite-existing-one-but-append-number #SOURCE;
}

sub Rscript {
	my $path = $ENV{'PWD'};
	my $R_code = <<"RCODE";
setwd("$path");
gg = read.table('$_[0]', header=T)
library(ggplot2)
png(filename="$_[1]", width = 1000, height = 500);
RCODE
	unless ($#variantfile >= 1){
		$R_code .= <<"RCODE";
myplot = ggplot(gg) +
geom_line(aes(x=position,y=count,group=chrom, color=chrom))+
facet_grid(.~chrom,scales="free")
RCODE
	} else {
		$R_code .= <<"RCODE";
library("reshape2")
hh <- melt(gg, id=c("chrom","position")) #convert to long format
myplot = ggplot(hh) +
geom_line(aes(x=position,y=value,group=chrom, color=variable)) +
facet_grid(.~chrom,scales="free")
RCODE
	}
$R_code .= <<"RCODE";
myplot
myplot + theme_bw()
dev.off()
RCODE

my $R = Statistics::R->new();
$R->startR;
$R->send($R_code);
$R->stopR();
}