#!/usr/bin/perl
use strict;
use IO::File;
use Pod::Usage;
use Getopt::Long;
use POSIX qw(ceil);
use Sort::Key::Natural qw(natsort);
use Statistics::R;

##### Usage Documentation
my $usage="
=> Measuring heterozygosity ratios following Rubin et al. Nature 2010 & others equation.
		and plotting the ZH scores using R.
Arguments needed for '$0'
\t-v|variant\t< Variant file (AF flag will be used) > 
\t\t\t[ Multiple variant files should be separated by comma ]
\t\tor
\t-t|table\t< Tab delimited file (with 1st col = CHROM, 2nd col = POSITION, and nth col = AF) >
\t\t\t[ Multiple table files should be separated by comma ]

Other required arguments:
\t-f|fasta\t< Reference fasta file >
\t-w|window\t< Sliding window size >
\t-s|step\t< Sliding window step >

Optional:
\t-c|chr <chromosome name (multiple separated by comma)>
\t-r|-region <chromosome region (specific location or region e.g 1-100000, only works if -chr is specified)>
\t-o|output <Output file name>

Copyright: 2017 MOA
Author: Modupeore Adetunji.
Institution: University of Delaware

";

##### Options of usage
my($variantfile, $tablefile, $fastafile,$outputfile,$window,$step, $chrom,$region);
my ($tables, $variants) = (0,0);
GetOptions('variant|v=s'=>\$variantfile, 't|table=s'=>\$tablefile, 'fasta|f=s'=>\$fastafile,
		'output|o=s'=>\$outputfile, 'w|window=s'=>\$window, 's|step=s'=>\$step, 'chr|chrom|c=s'=>\$chrom, 'region|r=s'=>\$region ) or die $usage;

die "Syntax Error$usage" unless ($fastafile && $window && $step);
die "Input file not provided$usage" unless ($variantfile || $tablefile);
die "Specify one type of input$usage" if ($variantfile && $tablefile);

####Logs
my $date = `date +%m-%d-%y_%T`; chomp $date;
my $std_out = "htrzy-$date.log";
my $std_err = "htrzy-$date.err";
open(STDOUT, '>', "$std_out") or die "Log file doesn't exist";
open(STDERR, '>', "$std_err") or die "Error file doesn't exist";

##### Global options
my ( %chrom, @variantfile, @heterozygosity, @regions);
my (%CHR, %GENOME, %REALGENOME, %header, %info, %GENENAME);
my (%AFmax, %AFmin, %SNPcount, %RUNTHRU, %NEWSTEP);

my %removehash  = map {$_ => 1} qw|W LGE64 16|; #chromosomes to remove

##### Code
print "... Processing genome ";
GENOME(); #read genome file
print "... Done ",`date`;

#if a specific chromosome is selected
if ($chrom){
	print "chromosomes selected $chrom\n";
	$chrom =~ s/\s+|\s+//g;
	%chrom = map {$_ => 1} (split (",", $chrom));
} else {
	print "Processing All chromosomes\n";
}

#if a specific region of the chromosome is selected
if ($region) {
	print "region specified $region\n";
	$region =~ s/\s+|\s+//g;
	if ($region =~ /\-/){
		@regions = split("\-", $region, 2);
	} else {
		push @regions, $region-int($window/2), $region+int($window/2);
	}
	if ($regions[0] < 1){$regions[0] = 1;}
	print "... Finished region ",`date`,"\n";
}

#reading the variant files or table files
if ($variantfile){
	print "Variant file selected $variantfile\n";
	@variantfile = split (",", $variantfile);
	$variants = 1;
}
if ($tablefile){
	print "Table file selected $tablefile\n";
	@variantfile = split (",", $tablefile);
	$tables = 1;
}

#sorting the chromosomes into steps.
print "... Processing steps for the chromosome with Window = $window and Step = $step ";
foreach my $gchr (natsort keys %REALGENOME){
	my $checklimit = 1;
	if ($chrom){
		if (exists $chrom{$gchr}) {
			my $steps = ceil($REALGENOME{$gchr}/$step);
			for (my $i = 0; $i <= $steps; $i++){
				my $newstep = $i * $step; if ($newstep == 0){$newstep = 1;}
				my $windows = $newstep + $window;
				if ($windows > $REALGENOME{$gchr}){ $windows = $REALGENOME{$gchr}; }
				if (($windows - $newstep) < $window){$checklimit = 2;}
				unless ($checklimit == 2){
					$RUNTHRU{$gchr}{$i} = "$newstep-$windows";
					$NEWSTEP{$gchr}{$newstep} = $i;
				}
			}
		}
	} else {
		my $steps = ceil($REALGENOME{$gchr}/$step);
		for (my $i = 0; $i <= $steps; $i++){
			my $newstep = $i * $step; if ($newstep == 0){$newstep = 1;}
			my $windows = $newstep + $window;
			if ($windows > $REALGENOME{$gchr}){ $windows = $REALGENOME{$gchr}; }
			if (($windows - $newstep) < $window){$checklimit = 2;}
			unless ($checklimit == 2){
				$RUNTHRU{$gchr}{$i} = "$newstep-$windows";
				$NEWSTEP{$gchr}{$newstep} = $i;
			}
		}
	}
} #end foreach real genome
print "... Done ",`date`,"\n";
#processing each variant file
foreach my $files (@variantfile) {
	$files =~ s/^\s+|\s+$//g;
	open(IN,"<",$files) or die "The file \'$files\' does not exist\n";
	#if ($tablefile) { print "yes"; die; } else { print "no"; die;}
	#die "main";
	if ($tablefile) {
		print "... Processing table file => $files";
		print "\n... Processing table file => $files";
		MAIN: while (<IN>){ chomp; 
			my $line = $_; 
			my @all = split (/\t/, $line); 
			if (exists $RUNTHRU{$all[0]}) { #checking the chromosome exists in the reference genome file
				my $ident = (int($all[1]/$step) * $step);
				my $newident = $ident - $window;
				my @array = map {$_; } sort {$a <=> $b} keys %{ $RUNTHRU{$all[0]} };
				if ($newident < 0){$newident = 0};
				foreach my $stepz ($NEWSTEP{$all[0]}{$newident}..@array){ #getting the index to search through
					my ($first, $last) = split(/\-/,$RUNTHRU{$all[0]}{$stepz},2);
					if ($all[1] >= $first && $all[1] <= $last) {
						my $AF = $all[$#all]; #computing AFs the last column
						unless ($AF =~ /,/){ #removing multiple AF for 1 position.
							if ($AF < 0.01 || $AF > 0.99) {print "homozygous allele @ $all[0]:$all[1] won't be computed\n";}
							else {
								#getting the gene information
								$all[4] =~ s/^\s+|\s+$//g;
								if (exists $GENENAME{$all[0]}{$stepz}){ 
									$GENENAME{$all[0]}{$stepz} = "$GENENAME{$all[0]}{$stepz},$all[4]";
								} else {
									$GENENAME{$all[0]}{$stepz} = $all[4];
								} # end else getting the genename
								#end of getting the gene information
								
								push my @newAF, $AF, 1-$AF;
								my ($maxAF, $minAF)= &range(\@newAF); #getting the maximum and minimum AF
								$AFmin{$all[0]}{$stepz}= $AFmin{$all[0]}{$stepz} + $minAF;
								$AFmax{$all[0]}{$stepz}= $AFmax{$all[0]}{$stepz} + $maxAF;
								$SNPcount{$all[0]}{$stepz}++;
							} #end if AF
						} #end unless multiple AF
					} #end if window region
					#next MAIN if $first > $all[1];
					#elsif ($first > $all[1]){next MAIN;}
				} #end foreach step
			} #print $line,"\n";
		} #end while
		print "... Done ",`date`;
	} #end if table file
	elsif ($variantfile){ #if variant file
		print "... Processing variant file => $files";
		MAINV: while (<IN>){ chomp;
			unless (/^#/){ #removing previous files
				my @all = split /\t/;
				if (exists $RUNTHRU{$all[0]}) {
					my $ident = (int($all[1]/$step) * $step);
					my $newident = $ident - $window;
					my @array = map {$_; } sort {$a <=> $b} keys %{ $RUNTHRU{$all[0]} };
					
					my $selectfilter = $all[$header{"INFO"}]; #getting the allele frequency from the INFO column
					my @info = split (/\;/, $selectfilter);
					undef %info; #clean the info hash
					foreach (@info) {
						my @theinfo = split("\=", $_, 2);
						$info{$theinfo[0]} = $theinfo[1];
					}
					
					if ($newident < 0){$newident = 0};
					foreach my $stepz ($NEWSTEP{$all[0]}{$newident}..@array){ #getting the index to search through
						my ($first, $last) = split(/\-/,$RUNTHRU{$all[0]}{$stepz},2);
						if ($all[1] >= $first && $all[1] <= $last) {
							my $AF = $info{'AF'}; #getting the AF.
							unless ($AF =~ /,/){ #removing multiple AF for 1 position
								if ($AF < 0.01 || $AF > 0.99) {print "homozygous allele @ $all[0]:$all[1] won't be computed\n";}
								else {
									push my @newAF, $AF, 1-$AF;
									my ($maxAF, $minAF)= &range(\@newAF);
									$AFmin{$all[0]}{$stepz}= $AFmin{$all[0]}{$stepz} + $minAF;
									$AFmax{$all[0]}{$stepz}= $AFmax{$all[0]}{$stepz} + $maxAF;
									$SNPcount{$all[0]}{$stepz}++;
								} #end if AF
							} #end unless multiple AF
						} #end if window region
						elsif ($first > $all[1]){next MAINV;}
					} #end foreach step
				}
			} #end unless
			if (/^#CH/) { #get the header column
				my @headerline = split /\t/;
				foreach (0..$#headerline){ $header{$headerline[$_]} = $_; }	
			}
		} #end while
		print "... Done ",`date`;
	} #end if variant file
	else { die "Not matching variables\n"; }
	my ($j,$i) = (0,0);
	#computing windows & printing OUT
	print "... Computing Heterozygosity ";
	$outputfile ||="output.txt"; my $fileout = @{ open_unique($outputfile) }[1];#outputfile name
	open (OUT, ">", $fileout) or die;
	print OUT "CHROM\tSTART\tEND\tGENE\tnumber\tSNPcount\tHeterozygosity\n";

	foreach my $schr (natsort %SNPcount){
		foreach my $stepy (sort {$a <=> $b} keys %{$SNPcount{$schr}}){
			my $finalAF = $AFmin{$schr}{$stepy} + $AFmax{$schr}{$stepy};
			unless ($SNPcount{$schr}{$stepy} < 10 ){
				my $heterozygosity = sprintf("%.3f",((2*$AFmin{$schr}{$stepy}*$AFmax{$schr}{$stepy})/($finalAF*$finalAF)));
				push @heterozygosity, $heterozygosity;
				my ($start, $end) = split(/\-/,$RUNTHRU{$schr}{$stepy});
				my %alls = map {$_ => 1} split(",", $GENENAME{$schr}{$stepy});
				my @array = map {$_; } sort {$a <=> $b} keys %alls ;
				my $genes = join (",", @array);
				print OUT "$schr\t$start\t$end\t$genes\t$finalAF\t$SNPcount{$schr}{$stepy}\t$heterozygosity\n";
				$i++;
				
				#unless ($SNPcount{$schr}{$stepy} == $finalAF){ print "Dont match s$SNPcount{$schr}{$stepy}","c\ts$finalAF","c\n"; }
			} else {
				#print OUT "$schr\t$RUNTHRU{$schr}{$stepy}\t$finalAF\t$SNPcount{$schr}{$stepy}\n";
				$j++;
			}
		}
	} #end foreach compute window
	print "... Done ",`date`;
	
	print "Number that passed SNPcount $i\nElse $j\n";
	close OUT;
	
	#computing ZH & printing OUT
	print "... Computing Z transformed heterozygosity ";
	my $mean = &average(\@heterozygosity);
	my $stddev = &stdev(\@heterozygosity);
	my $fileout2 = "ZHresult-".$fileout; #outputfile name
	open (OUT2, ">", $fileout2) or die;
	print OUT2 "CHROM\tSTART\tEND\tGENE\tnumber\tSNPcount\tHeterozygosity\tZHeterozygosity\n";
	
	foreach my $schr (natsort %SNPcount){
		foreach my $stepy (sort {$a <=> $b} keys %{$SNPcount{$schr}}){
			my $finalAF = $AFmin{$schr}{$stepy} + $AFmax{$schr}{$stepy};
			unless ($SNPcount{$schr}{$stepy} < 10 ){
				my $heterozygosity = sprintf("%.3f",((2*$AFmin{$schr}{$stepy}*$AFmax{$schr}{$stepy})/($finalAF*$finalAF)));
				my $zh = sprintf("%.3f",(($heterozygosity - $mean)/$stddev));
				my ($start, $end) = split(/\-/,$RUNTHRU{$schr}{$stepy});
				my %alls = map {$_ => 1} split(",", $GENENAME{$schr}{$stepy});
				my @array = map {$_; } sort {$a <=> $b} keys %alls ;
				my $genes = join (",", @array);
				print OUT2 "$schr\t$start\t$end\t$genes\t$finalAF\t$SNPcount{$schr}{$stepy}\t$heterozygosity\t$zh\n";
			} else {
				#print OUT2 "$schr\t$RUNTHRU{$schr}{$stepy}\t$finalAF\t$SNPcount{$schr}{$stepy}\n";
			}
		}
	} #end foreach compute window
	close OUT2;
	print "... Done ",`date`;
	
	#R scripts
	open(IN, $fileout2) or die "Can't open $fileout2\n"; my @inputall = <IN>; close(IN);
	$fileout2 = "forR-".$fileout; #outputfile name
	open (OUT, ">", $fileout2) or die;
	my $originalfileout = "original-$fileout"; open (OUT3, ">", $originalfileout) or die;
	my $imptfileout = "impt-$fileout"; open (OUT2, ">", $imptfileout) or die;
	my $header = shift @inputall;
	$i = 0;
	print OUT "No\t$header"; print OUT3 "No\t$header";
	foreach my $for (@inputall){
	  chomp $for;
	  my @id = split(/\t/, $for);
	  if ($id[0] =~ /chr/) { $id[0] =~ s/chr//ig };
	  $id[0] = uc($id[0]);
	  unless (exists $removehash{$id[0]}){
			$i++; 
	    print OUT $i,"\t",$id[0],"\t"; print OUT3 $i,"\t",$id[0],"\t";
	    foreach my $j (1..$#id-1){ print OUT $id[$j],"\t"; print OUT3 $id[$j],"\t";}
	    print OUT3 $id[$#id],"\n";
			if ($id[$#id] > 0){ $id[$#id] = $id[$#id]*-1; }
	    if ($id[$#id] <= -4) {print OUT2 $for,"\n"; }
	    print OUT $id[$#id],"\n";
	  }
	}
	close (OUT); close (OUT2); close(OUT3);
	
	#picture
	my $picture = "ZH-$fileout"; my $picture2 = "plot-$fileout";
	$picture =~ s/\..+$/\.png/g; $picture2 =~ s/\..+$/\.png/g;
	my $path = $ENV{'PWD'};
	my $R_code = <<"RCODE";
setwd("$path");
library(ggplot2);
png(filename="$picture", width=960, height=480, res = 200);
info = read.table("$fileout2", header=T)
info\$CHROM <- factor(info\$CHROM, levels=unique(as.character(info\$CHROM)) )
p <- ggplot(data = info, aes(x = No, y = ZHeterozygosity, color=CHROM)) +
  geom_point(stat = "identity") + labs(x = "Chromosome", y=expression(ZH[ W]))
if (min(info\$ZHeterozygosity) <= -6) { p = p + geom_hline(yintercept=c(-4,-6), color="black", linetype="dashed") 
} else { p = p + geom_hline(yintercept=-4, color="black", linetype="dashed") }
p + facet_grid(~CHROM, scales = 'free_x', space = 'free_x', switch = 'x') + 
  theme_classic() + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0, "lines")) +
  theme(legend.position="none") +
  scale_x_continuous(expand = c(0, 0))
dev.off()

png(filename="$picture2", width=480, height=480, res = 200);
info = read.table("$originalfileout", header=T)
temp <- paste("mu == ", 0)
other <- paste("sigma == ", 1)
h <- ggplot(info,aes(x = ZHeterozygosity)) + 
  theme_classic() + 
  theme(legend.position="none") +
  geom_histogram(binwidth = .5,aes(fill =..count..)) +
  scale_fill_gradient("", low = "darkgrey", high = "navy") +
  labs(x=expression(ZH[ W]), y = "frequency")
h
dev.off()
RCODE

	my $R = Statistics::R->new();
	$R->startR;
	$R->send($R_code);
	$R->stopR();	
	
} #end foreach input file

print "... finished\n";
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
			$total = $total + length($pieces[$num]);
		}
		#my $newtotal = int($total/$bin)*$bin; if ($newtotal < $total){ $newtotal= $newtotal+$bin;}
		#$GENOME{$pieces[0]} = $newtotal;
		$REALGENOME{$pieces[0]} = $total;
		#print "$pieces[0]\t$total\t$newtotal\n";
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
