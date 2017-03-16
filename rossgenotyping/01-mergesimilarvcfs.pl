#!/usr/bin/perl
use strict;
use File::Basename;
use Sort::Key::Natural qw(natsort);
# merge similar variants based on number of vcfs provided

#- - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
print "\t**MERGE SIMILAR VARIANTS**\n";
my $usage = "
To use : '$0'
\ttype in the vcfs files (more than 1) separated by space.\n
";
@ARGV >= 2 or die $usage;

# - - - - - G L O B A L V A R I A B L E S - - - - - - - - -

my $GATK = "/home/modupe/.software/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar";
my $REF = "/home/modupe/.big_ten/chicken/chicken/chicken.fa";

# - - - - - U S E R V A R I A B L E S - - - - - - - - - - -

my (%comparison, %commonfile, %inputfile, %GATKs); #, @new);
my $i = 0;
my ($line, $header);
my $out = "merge-".fileparse($ARGV[0], qr/\.[^.]*(\.vcf)?$/).".vcf";
open (OUT,">", $out);
# - - - - - M A I N - - - - - - - - - - - - - - - - - - -

open (COMMON, $ARGV[0]) or die "Can't open '$ARGV[0]'\n$usage";
while (<COMMON>){
	if (/^chr.*/){
		$i++;
		$line = $_;
    my @commonline = split (/\t/, $line);
		#push @new, $commonline[0], $commonline[1];
		#my $new = "$commonline[0]\.$commonline[1]";
		$GATKs{$i}{$commonline[0]}{$commonline[1]} = $line;
		#"$commonline[0]\t$commonline[1]";
    $commonfile{$commonline[0]}{$commonline[1]} = $line;
  } elsif (/^#CH/) {
		$header = $_;
	} else {
		print OUT $_;
	}
}
close COMMON;
print OUT "##SUB=<ID=Merged,Parent=\"$ARGV[0]\"";
print "Total number of variants in the parent file '$ARGV[0]' : $i\n"; $i = 0;
my $GATKstd = "java -jar $GATK -T CombineVariants -R $REF --variant $ARGV[0] ";

#Working on the other vcf files
foreach my $no (1..$#ARGV){
	open (FILE, $ARGV[$no]) or die "Can't open '$ARGV[$no]'\n$usage";
	while (<FILE>){
		if (/^chr.*/){
			$line = $_;
			my @commonline = split (/\t/, $line);
			my $new = "$commonline[0]\.$commonline[1]";
			if (exists $commonfile{$commonline[0]}{$commonline[1]}) {
				$i++;
				$inputfile{$commonline[0]}{$commonline[1]}= $commonfile{$commonline[0]}{$commonline[1]};
			}
		}
	}
	close FILE;
	print OUT ",Branch=\"$ARGV[$no]\"";
	print "Total number of variants that match from '$ARGV[$no]' : $i\n"; $i = 0;
	undef %commonfile; %commonfile = %inputfile; undef %inputfile; #move common to the input file
	$GATKstd .= "--variant $ARGV[$no] ";
}
$GATKstd .= "-o gatk_$out -genotypeMergeOptions UNIQUIFY";
print OUT ">\n$header";
#Output
foreach my $newkey (sort {$a <=> $b} keys %GATKs) {
	foreach my $key1 (keys %{$GATKs{$newkey}}) {
		foreach my $key2 (keys %{$GATKs{$newkey}{$key1}}) {
			if (exists $commonfile{$key1}{$key2}) { print OUT $commonfile{$key1}{$key2}; }
		}
	}
}
#foreach my $key (natsort keys %commonfile) {
#	foreach my $key2 (sort {$a <=> $b} keys %{$commonfile{$key}}) {
#		print OUT $commonfile{$key}{$key2};
#	}
#}
close OUT;
`$GATKstd 2> .error`; #run GATK combine variants

#Select SNPs
#`java -jar $GATK -T CombineVariants -R $REF --variant $out --variant $out -o sort_$out 2>> .error`;
`java -jar $GATK -T SelectVariants -selectType SNP -R $REF -V $out -o snps-$out 2>> .error`;
`java -jar $GATK -T SelectVariants -selectType SNP -R $REF -V gatk_$out -o snps-gatk_$out 2>> .error`;

#Select INDELs
`java -jar $GATK -T SelectVariants -selectType INDEL -R $REF -V $out -o indel-$out 2>> .error`;
`java -jar $GATK -T SelectVariants -selectType INDEL -R $REF -V gatk_$out -o indel-gatk_$out 2>> .error`;

print "Merged variant file output : '$out'\n";
print "Merged variant file output 'GATK' : 'gatk_$out'\n";
print "SNP-Merged variant file output : 'snp-$out'\n";
print "SNP-Merged variant file output 'GATK' : 'gatk_snp-$out'\n";
print "InDEL-Merged variant file output : 'indel-$out'\n";
print "InDEL-Merged variant file output 'GATK' : 'gatk_indel-$out'\n";
