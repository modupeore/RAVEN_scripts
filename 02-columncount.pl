#!/usr/bin/perl
use strict;
# to extract similar or unique based on 1st column

#- - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
print "\t**SELECTING THE SIMILAR AND NOT SIMILAR column 1**\n\n";
my $usage = "
For $0
Type in the name of the :
\t\t1.\"common file\".
\t\t2. The \"testfile\".
\t\t3. The unique_outputfilename.
\t\t4. The common_outputfilename
";
@ARGV >= 2 or die $usage;

# - - - - - U S E R V A R I A B L E S - - - - - - - - - -
my (%comparison, %commonfile, %inputfile);
my ($inputcommon, $input, $output_only, $output_common);
# the comparison file
$inputcommon = $ARGV[0]; open (FILE, "<$inputcommon") or die "Cannot find file $ARGV[0]\n$usage\n";
while (<FILE>){	chomp;	my @array = (split(" ",$_,2))[0];    $commonfile{$array[0]} = $array[1];	}	close FILE;
#the test file
$input = $ARGV[1]; open (INPUTFILE, "<$input") or die "Cannot find file $input\n";
while (<INPUTFILE>){	chomp;	my @array = (split(" ",$_,2))[0];   $inputfile{$array[0]} = $array[1];  } close INPUTFILE;

# the output files
if ($ARGV[2] || $ARGV[3]) {
$output_only = $ARGV[2]; open (OUTPUTFILE, ">$output_only");
$output_common = $ARGV[3]; open(OUTPUTFILE2, ">$output_common");
}
# - - - - - G L O B A L V A R I A B L E S - - - - - - - - -
# counting the reads
my ($unique_count, $common_count, $thesame) = (0,0,0);
# common variables
my $commonlist;
#test file variables
my $list;
my %sequencename;
# all variables
my @concat_list_commonlist;
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - M A I N - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

#1. getting the list for the common regions
foreach my $commonline (sort keys %commonfile){
    if (exists $inputfile{$commonline}) {
        if ($commonfile{$commonline} == $inputfile{$commonline}){
            $thesame++;
        }
        if ($output_common) {print OUTPUTFILE2 $commonline;}
        $common_count++;
    } else {
        if ($output_only) {print OUTPUTFILE $commonline;}
        $unique_count++;
    }
}

#            OUTPUTFILE $sequencename{$commonline} = $i;
#}
## 2. getting the list for the test file
#foreach my $inputline (@inputfile){
#    if (exists $sequencename{$inputline}){
#            print OUTPUTFILE2 $inputline;
#            $common_count++;
#    } else {
#            print OUTPUTFILE $inputline;
#            $unique_count++;
#    }
#}
print "Total number of unique reads in \"$ARGV[1]\" that is not in \"$ARGV[0]\" are \"$unique_count\".\n";
print "Total number of common reads in both \"$ARGV[0]\" and \"$ARGV[1]\" are \"$common_count\".\n";
print "The same count\t'$thesame'.\n";
if ($ARGV[2] || $ARGV[3]) {
    print "Successfully saved in the outputfile $output_only & $output_common";
    close OUTPUTFILE;
    close OUTPUTFILE2;
}
print "\n\n*****************DONE*****************\n\n";
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - S U B R O U T I N E S - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -

