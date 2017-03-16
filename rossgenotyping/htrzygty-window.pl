#!/usr/bin/perl
use strict;
use IO::File;
use Pod::Usage;
use Getopt::Long;
use POSIX qw(ceil);
use Sort::Key::Natural qw(natsort);

##### Usage Documentation
my $usage="
=> Measuring heterozygosity ratios following Rubin et al. Nature 2010 equation.
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
my($variantfile, $tablefile, $fastafile,$outputfile,$window, $step, $chrom,$region);
GetOptions('variant|v=s'=>\$variantfile, 't|table=s'=>\$tablefile, 'fasta|f=s'=>\$fastafile,
		'output|o=s'=>\$outputfile, 'w|window=s'=>\$window, 's|step=s'=>\$step, 'chr|chrom|c=s'=>\$chrom, 'region|r=s'=>\$region ) or die $usage;

die $usage unless ($fastafile && $window && $step);
die $usage unless ($variantfile || $tablefile);
die $usage if ($variantfile && $tablefile);

####Logs
my $date = `date +%m-%d-%y_%T`; chomp $date;
my $std_out = "htrzy-$date.log";
my $std_err = "htrzy-$date.err";
open(STDOUT, '>', "$std_out") or die "Log file doesn't exist";
open(STDERR, '>', "$std_err") or die "Error file doesn't exist";

##### Global options
my ( %chrom, @variantfile, @heterozygosity, @regions);
my (%CHR, %GENOME, %REALGENOME, %header, %info);
my (%AFmax, %AFmin, %SNPcount, %RUNTHRU, %NEWSTEP);

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
}
if ($tablefile){
	print "Table file selected $tablefile\n";
	@variantfile = split (",", $tablefile);
}

#sorting the chromosomes into steps.
print "... Processing steps for the chromosome ";
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
	open(IN,"<", $files) or die "The file '$files' does not exist\n";
	if ($tablefile) { #if table file
		print "... Processing table file => $files";
		MAIN: while (<IN>){ chomp;
			my @all = split /\t/; 
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
								push my @newAF, $AF, 1-$AF;
								my ($maxAF, $minAF)= &range(\@newAF); #getting the maximum and minimum AF
								$AFmin{$all[0]}{$stepz}= $AFmin{$all[0]}{$stepz} + $minAF;
								$AFmax{$all[0]}{$stepz}= $AFmax{$all[0]}{$stepz} + $maxAF;
								$SNPcount{$all[0]}{$stepz}++;
							} #end if AF
						} #end unless multiple AF
					} #end if window region
					elsif ($first > $all[1]){next MAIN;}
				} #end foreach step
			}
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
	my ($j,$i) = (0,0);
	#computing windows & printing OUT
	print "... Computing Heterozygosity ";
	$outputfile ||="output.txt"; $outputfile = @{ open_unique($outputfile) }[1]; $outputfile =~ /^(.*)\..+$/; #outputfile name
	open (OUT, ">", $outputfile) or die;
	print OUT "CHROM\tSTART\tEND\tnumber\tSNPcount\tHeterozygosity\n";

	foreach my $schr (natsort %SNPcount){
		foreach my $stepy (sort {$a <=> $b} keys %{$SNPcount{$schr}}){
			my $finalAF = $AFmin{$schr}{$stepy} + $AFmax{$schr}{$stepy};
			unless ($SNPcount{$schr}{$stepy} < 10 ){
				my $heterozygosity = sprintf("%.3f",((2*$AFmin{$schr}{$stepy}*$AFmax{$schr}{$stepy})/($finalAF*$finalAF)));
				push @heterozygosity, $heterozygosity;
				my ($start, $end) = split(/\-/,$RUNTHRU{$schr}{$stepy});
				print OUT "$schr\t$start\t$end\t$finalAF\t$SNPcount{$schr}{$stepy}\t$heterozygosity\n";
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
	$outputfile ||="output.txt"; $outputfile = @{ open_unique($outputfile) }[1]; $outputfile =~ /^(.*)\..+$/; #outputfile name
	open (OUT2, ">", $outputfile) or die;
	print OUT2 "CHROM\tSTART\tEND\tnumber\tSNPcount\tHeterozygosity\tZHeterozygosity\n";
	
	foreach my $schr (natsort %SNPcount){
		foreach my $stepy (sort {$a <=> $b} keys %{$SNPcount{$schr}}){
			my $finalAF = $AFmin{$schr}{$stepy} + $AFmax{$schr}{$stepy};
			unless ($SNPcount{$schr}{$stepy} < 10 ){
				my $heterozygosity = sprintf("%.3f",((2*$AFmin{$schr}{$stepy}*$AFmax{$schr}{$stepy})/($finalAF*$finalAF)));
				my $zh = sprintf("%.3f",(($heterozygosity - $mean)/$stddev));
				my ($start, $end) = split(/\-/,$RUNTHRU{$schr}{$stepy});
				print OUT2 "$schr\t$start\t$end\t$finalAF\t$SNPcount{$schr}{$stepy}\t$heterozygosity\t$zh\n";
			} else {
				#print OUT2 "$schr\t$RUNTHRU{$schr}{$stepy}\t$finalAF\t$SNPcount{$schr}{$stepy}\n";
			}
		}
	} #end foreach compute window
	close OUT2;
	print "... Done ",`date`;
} #end foreach input file

#		
#	}
#	}
#}
#open (OUT, ">$outputfile"); #save to outputfile
#print OUT "chrom\tposition";
#unless ($#variantfile >= 1){
#	print OUT "\tcount\n";
#} else { 
#	foreach my $files (@variantfile){ $files =~ s/^\s+|\s$//g; print OUT "\t$files"; }
#	print OUT "\n";
#}
#foreach my $a (natsort keys %CHR) {
#	for (my $x = 0; $x <= $GENOME{$a}; $x=$x+$window){
#		my $v;if ($x ==0) { $v = 1; } else { $v = $x; }
#		print OUT "$a";
#		print OUT "($REALGENOME{$a})";
#		print OUT "\t$v";
#		foreach my $files (@variantfile){
#			$files =~ s/^\s+|\s$//g;
#			if (exists $CHR{$a}{$v}{$files}){ print OUT "\t",$CHR{$a}{$v}{$files}; } else { print OUT "\t0"; }
#		}
#		print OUT "\n";
#	}
#}
#close (OUT);


#Rscript($outputfile, $picture);

#my $tempout = @{ open_unique(".tempout.txt")}[1]; #ploting 10 chrs in one pdf
#my ($count, $counter) = (1,1);
#my $pictures= $index."-".$count.".png";
#open (TEMPOUT, ">$tempout");
#print TEMPOUT "chrom\tposition";
#unless ($#variantfile >= 1){
#	print TEMPOUT "\tcount\n";
#} else { 
#	foreach my $files (@variantfile){ $files =~ s/^\s+|\s$//g; print TEMPOUT "\t$files"; }
#	print TEMPOUT "\n";
#}
#foreach my $a (natsort keys %CHR) {
#	unless ($a eq $new) {
#		if ($counter > 10) {
#			$pictures= $index."-".$count.".png";
#			$count++;
#			Rscript($tempout, $pictures);
#			open (TEMPOUT, ">$tempout");
#			print TEMPOUT "chrom\tposition";
#			unless ($#variantfile >= 1){
#				print TEMPOUT "\tcount\n";
#			} else { 
#			foreach my $files (@variantfile){ $files =~ s/^\s+|\s$//g; print TEMPOUT "\t$files"; }
#			print TEMPOUT "\n";
#			}
#			$counter = 1;
#		}
#		else {
#			$counter++;
#		}
#		$new = $a;
#	}
#	for (my $x = 0; $x <= $GENOME{$a}; $x=$x+$bin){
#		my $v;if ($x ==0) { $v = 1; } else { $v = $x; }
#		print TEMPOUT "$a";
#		print TEMPOUT "($REALGENOME{$a})";
#		print TEMPOUT "\t$v";
#		foreach my $files (@variantfile){
#			$files =~ s/^\s+|\s$//g;
#			if (exists $CHR{$a}{$v}{$files}){ print TEMPOUT "\t",$CHR{$a}{$v}{$files}; } else { print TEMPOUT "\t0"; }
#		}
#		print TEMPOUT "\n";
#	}
#}
#close (TEMPOUT);
#$pictures= $index."-".$count.".png";
#Rscript($tempout,$pictures);
#`rm -rf $tempout`;


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
