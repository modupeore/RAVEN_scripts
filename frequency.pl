#! usr/bin/perl
use File::Basename;

#OBJECTIVE
#continuation of tabulationSNPS-2.pl
print "\n\tFrequency of annotations called from Tabulated file\n";

my $input = $ARGV[0];
my $i= 0;
unless(open(FILE,$input)){
	print "File \'$input\' doesn't exist\n";
	exit;
}

my @file = <FILE>;
chomp @file;
close (FILE);
my $count = 0;
my %Hashdetails;
$ind = $ARGV[1];
foreach my $chr (@file){
	if ($chr =~ /^#/){
		next;
	}
	my @chrdetails = split('\t', $chr);
	my $chrIwant = $chrdetails[$ind];
	#my @morechrsplit = split('\(', $chrIwant);
	$Hashdetails{$chrIwant} = 0;
	
}

foreach my $newcount (@file){
        if ($newcount =~ /^#/){
		next;
	}
        my @details = split('\t', $newcount);
        my $chragain = $details[$ind];
        #my @morechr = split('\(', $chragain);
        $Hashdetails{$chragain} ++; $count++;
	#print $count++."\t$details[0]\t$details[1]\t\t$morechr[0]\t$Hashdetails{$morechr[0]}\n";
        
}

print "ABC\n";
foreach my $newkey (sort keys %Hashdetails){
	print "$newkey\t$Hashdetails{$newkey}\n";
}

print "\n$count\n";
exit;
