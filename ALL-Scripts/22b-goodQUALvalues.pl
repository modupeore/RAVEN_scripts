#! usr/bin/perl
use File::Basename;

#OBJECTIVE
#filtering DP > 5.
#saving header information in another file

opendir(DIR,$ARGV[0]) or die "Folder \"$fastqfolder\" doesn't exist\n";
my @Directory = readdir(DIR);
close(DIR);
my $final = "all-there.data";
open(FINAL,">$final");
print FINAL "Lib. No\tMin.\t1st Qu.\tMedian\tMean\t3rd Qu.\tMax.\n";
print "Total SNPs\tLib. No\tMin.\t1st Qu.\tMedian\tMean\t3rd Qu.\tMax.\n";
		
foreach my $FILE (@Directory){ 
	if ($FILE =~ /^l.*/){
		my $input = $ARGV[0]."/".$FILE."/".$FILE."\_DP5\.vcf";
		unless(open(FILE,$input)){
			print "File \'$input\' doesn't exist\n";
			next;
		}
		$number++;
		#parsing NAME of file
		my $out = fileparse($input, qr/(\.vcf)?$/);
		
		my $output = "$out".".txt"; 
		open(OUT,">$output");
		
		my $Routput = "$out".".R";
		open(ROUT,">$Routput");
	
		my %DPvalue = '';

		my @file = <FILE>;
		chomp @file;
		close (FILE);
		my $count = 0;
		foreach my $chr (@file){
			if ($chr =~ /^chr/){
				my @chrdetails = split('\t', $chr);
				my $chrIwant = $chrdetails[5];
				
						$count++;
						$DPvalue{$count} = $chrIwant;
			}
		}
		print OUT "Number\tDPvalue\n";
		foreach my $i ( keys %DPvalue){
			print OUT "$i\t$DPvalue{$i}\n";
		}
		close (out);
		#running R
		print ROUT "d <- read.table(\"$output\",header=T,sep=\"\\t\")\n"; 
		print ROUT "min(d\$DPvalue, na.rm = T)\n";
		print ROUT "quantile(d\$DPvalue,0.25, names = F, na.rm = T)\n";
		print ROUT "median(d\$DPvalue, na.rm = T)\n";
		print ROUT "mean(d\$DPvalue, na.rm = T)\n";
		print ROUT "quantile(d\$DPvalue,0.75, names = F, na.rm = T)\n";
		print ROUT "max(d\$DPvalue, na.rm = T)\n";
		close (ROUT);

		#processing R
		my $Rresult = `/usr/bin/Rscript $Routput`;
		my @allR = split('\n', $Rresult);
		my $realR1 = substr $allR[0], 4;
		my $realR2 = substr $allR[1], 4;
		my $realR3 = substr $allR[2], 4;
		my $realR4 = substr $allR[3], 4;
		my $realR5 = substr $allR[4], 4;
		my $realR6 = substr $allR[5], 4;

		print FINAL "$FILE\t$realR1\t$realR2\t$realR3\t$realR4\t$realR5\t$realR6\n";
		print "$count\t$FILE\t$realR1\t$realR2\t$realR3\t$realR4\t$realR5\t$realR6\n";
		`rm -rf $Routput $output`;
		
	}
}
print $number;

exit;