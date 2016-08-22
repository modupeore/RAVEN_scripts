#! usr/bin/perl
use File::Basename;

#OBJECTIVE
#working on the genes output
print "\n\tGetting the correct gene_ids from the genes.gff file\n";

my $input = $ARGV[0];
my $out = fileparse($input, qr/\.[^.]*(\.genes)?$/);
my $output = "$out"."-newgeneids.txt";
my $output2 = "$out"."-converted.old";
open(OUT,">$output");
open(OUT2,">$output2");

my $i= 0;
unless(open(FILE,$input)){
	print "File \'$input\' doesn't exist\n";
	exit;
}

my @file = <FILE>;
chomp @file;
close (FILE);
my $count = 0;
my $check = 0;
foreach my $genes (@file){
	`egrep \"=$genes\" /home/amodupe/CARL/PERMANENT/genes.gff > temp.txt`;
	open(TEMP,"<temp.txt");
	my @temporary = <TEMP>;
	chomp @temporary;
	close (TEMP);
	$check = 0;
	foreach my $tempoutput (@temporary){
		if ($tempoutput =~ /Dbxref.*GeneID\:(\d*)/ && $check == 0){
			print OUT "$1\n"; $count++; $check = 1;
			print OUT2 "$genes\n";
		}
	}
}

print "\n$count\n";
`rm -rf temp.txt`;
close (OUT); close (OUT2);
exit;
