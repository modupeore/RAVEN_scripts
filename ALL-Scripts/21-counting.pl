#! usr/bin/perl
use File::Basename;

#OBJECTIVE

print "\n\tFrequency of variants called per chromosome \n";

my $input = $ARGV[0];
my $id= $ARGV[1];
unless(open(FILE,$input)){
	print "File \'$input\' doesn't exist\n";
	exit;
}

my @file = <FILE>;
chomp @file;
close (FILE);
my $count = 0;
my %Hashdetails;
foreach my $chr (@file){
	if ($chr !~ /^\#/){
		my @chrdetails = split('\t', $chr);
		my $chrIwant = $chrdetails[$id];
		$Hashdetails{$chrIwant} = 0;
	}
}

foreach my $newcount (@file){
        if ($newcount !~ /^\#/){
                my @details = split('\t', $newcount);
                my $chragain = $details[$id];
                $Hashdetails{$chragain}++; $count++;
        }
}

print "The distribution of variants per chromosome\n";
foreach my $newkey (sort keys %Hashdetails){
	print "$newkey\t$Hashdetails{$newkey}\n";
}

print "\nTotal number of variants : $count\n";
exit;
