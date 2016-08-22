#!/usr/bin/perl

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# MANUAL FOR insert_results.pl

=pod

=head1 NAME

$0 -- Insert results from tophat and cufflinks analysis into the mysql database : PENGUIN

=head1 SYNOPSIS

insert_results.pl --in FOLDER [--help] [--manual]

=head1 DESCRIPTION

Accepts the folders containing the output from "top_cuff_analysis.pl"
 
=head1 OPTIONS

=over 3

=item B<-i, -1, -a, --in, --folder>=FOLDER

Results folder (e.g. the folder with results to insert into the database).  (Required)

=item B<-h, --help>

Displays the usage message.  (Optional) 

=item B<-man, --manual>

Displays full manual.  (Optional) 

=back

=head1 DEPENDENCIES

Requires the following Perl libraries (all standard in most Perl installs).
   DBI
   DBD::mysql
   Getopt::Long
   Pod::Usage

=head1 AUTHOR

Written by Modupe Adetunji, 
Center for Bioinformatics and Computational Biology Core Facility, University of Delaware.

=head1 REPORTING BUGS

Report bugs to amodupe@udel.edu

=head1 COPYRIGHT

Copyright 2013 MOA.  
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.  
This is free software: you are free to change and redistribute it.  
There is NO WARRANTY, to the extent permitted by law.  

Please acknowledge author and affiliation in published work arising from this script's 
usage <http://bioinformatics.udel.edu/Core/Acknowledge>.

=cut

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - U S E R  V A R I A B L E S- - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# CODE FOR insert_results.pl

use strict;
use DBI;
use DBD::mysql;
use Getopt::Long;
use Pod::Usage;

#ARGUMENTS
my($in1,$help,$manual);
GetOptions (	
# 				"i|1|a|in|folder=s"	=>	\$in1,
                                "h|help"        =>      \$help,
                                "man|manual"	=>      \$manual );

# VALIDATE ARGS
pod2usage( -verbose => 2 )  if ($manual);
pod2usage( -verbose => 1 )  if ($help);
# pod2usage( -msg  => "ERROR!  Required argument -1 (input folder) not found. \n", -exitval => 2, -verbose => 1)  if (! $in1 );

#making sure the input file is parsable
my @temp = split('',$in1); $in1 = undef; my $checking = pop(@temp); push (@temp, $checking); unless($checking eq "/"){ push (@temp,"/")}; foreach(@temp){$in1 .= $_};

# DATABASE ATTRIBUTES
my $dsn = 'dbi:mysql:PENGUIN';
my $user = 'modupe';
my $passwd = 'penguin123';

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - G L O B A L  V A R I A B L E S- - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# FOLDER VARIABLES
my ($top_folder, $cuff_folder);

# RESULTS_HASH
my %Hashresults; my $number=0;

# DATABASE VARIABLES
my ($dbh, $sth, $syntax, $row, @row);
my $code="DE"; 
#PARSING VARIABLES
my @parse; my $len;

# TABLE VARIABLES
my ($lib_id, $total, $mapped, $unmapped, $deletions, $insertions, $junctions, $genes, $isoforms, $prep, $date); #RESULTS_SUMMARY TABLE
my ($raw_reads, $fastqc, $accepted, $unmapped_bam, $deletions_bed, $insertions_bed, $junctions_bed, $skipped_gtf, $transcripts_gtf, $isoforms_fpkm, $genes_fpkm, $run_log); #ACTUAL_FILES TABLE
my ($track, $class, $ref_id, $gene, $gene_name, $tss, $locus, $chrom_no, $chrom_start, $chrom_stop, $length, $coverage, $fpkm, $fpkm_low, $fpkm_high, $fpkm_stat); # GENES_FPKM & ISOFORMS_FPKM TABLE


my @newarray = ( 
	123
	# 123,124,125,126,127,128,133,134,135,136,
	# 137,138,139,140,141,142,143,144,145,146,
	# 147,148,152,153,154,155,156,157,159,15,
	# 160,161,162,163,164,16,172,173,174,175,
	# 176,177,179,17,180,181,182,183,184,186,
	# 187,188,189,18,190,191,19,200,201,203,
	# 204,205,206,207,208,20,210,211,212,213,
	# 214,215,217,218,219,220,221,222,224,225,
	# 226,227,228,229,231,232,234,235,236,238,
	# 239,240,241,242,243,245,246,247,248,249,
	# 250,252,256,257,258,259,260,261,262,263,
	# 264,265,266,267,26,270,271,272,273,274,
	# 275,276,277,278,279,280,281,282,283,284,
	# 285,286,287,288,289,290,291,292,293,299,
	# 320,321,322,323,325,326,327,328,329,32,
	# 330,332,333,334,335,336,337,339,340,341,
	# 342,343,344,347,348,349,350,351,352,354,
	# 355,356,357,358,359,361,362,363,364,365,
	# 366,372,373,375,376,377,378,379,380,381,
	# 382,383,384,385,
	
	#38109.pts-10.raven
	# 386,387,388,389,390,391,
	# 392,398,399,400,401,402,403,404,405,406,
	
	#22150.pts-22.raven
	# 407,408,409,429,430,431,437,438,439,440,
	
	
	# 443,445,446,447,448,449,450,451,453,454,
	# 455,456,457,458,459,460,468,469,470,471,
	# 472,476,477,480,481,484,485,486,487,488,
	# 489,490,491,492,493,494,496,497,498,
	
	#12744.pts-10.raven
	# 499,
	# 500,501,503,504,505,506,507,508,510,511,
	# 512,513,514,515,522,523,524,525,526,527,
	
	#15184.pts-10.raven
	# 529,530,531,533,537,538,539,540,541,542,
	# 543,544,546,547,549,550,551,552,553,554,
	
	
	# 556,557,558,559,560,562,563,564,565,566,
	# 568,569,570,572,573,576,577,578,579,580,
	# 581,583,584,585,586,587,588,590,591,592,
	# 593,594,595,611,612,613,614,615,616,617,
	# 618,619,620,621,622,623,626,627,628,629,
	# 630,631,632,634,635,636,637,638,647,648,
	# 649,650,651,652,653,654,655,656,657,658,
	# 661,662,663,664,665,666,667,668,669,66,
	# 670,671,672,675,676,677,678,679,67,680,
	# 681,682,684,685,686,687,688,689,701,702,
	# 703,704,705,706,708,709,70,710,711,712,
	# 713,715,716,717,718,719,71,720,721,728,
	# 72,739,73,774,79,807,808,809,80,810,
	# 811,814,819,81,825,826,827,828,82,834,
	# 836,839,83,841,842,843,844,845,846,84,
	# 863,864,865,867,868,869,879,880,881,883,
	# 884,886,894,895,897,899,900,901,903,912,
	# 913,914,915,917,918,919,921,925,926,927,
	# 928,931,932,933,935,936,955,956,957,958,
	# 95,960,973,
	
	#24608.pts-25.raven
	# 974,976,977,979,981,983,984,
	# 986,987,98,990,992,993,994,997,99,725,
	# 743,922
);



# CONNECT TO THE DATABASE
print "\n\n\tCONNECTING TO THE DATABASE : $dsn\n\n";
$dbh = DBI->connect($dsn, $user, $passwd) or die "Connection Error: $DBI::errstr\n";


foreach my $spec (@newarray){
	print "working on $spec\n\n";
	my $outfile = "/home/modupe/genesfpkm/$spec"."genes.fpkm_tracking";
	open(GENES, "<$outfile") or die "Can't open file $outfile\n";
	$lib_id = $spec."DE";

# # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# # - - - - - - - - - - - - - - - - M A I N - - - - - - - - - - - - - - - - - - - - - -
# # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 
	print ". . .\n\tINSERTING INTO THE DATABASE IN \"GENES_FPKM\" TABLE\n\n";
	$sth = $dbh->prepare("insert into GENES_FPKM (library_id, tracking_id, class_code, nearest_ref_id, gene_id, 
		gene_short_name, tss_id, chrom_no, chrom_start, chrom_stop, length, coverage, fpkm, fpkm_conf_low, 
		fpkm_conf_high, fpkm_status ) values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)");
	while (<GENES>){
	    chomp;
		my ($track, $class, $ref_id, $gene, $gene_name, $tss, $locus, $length, $coverage, $fpkm, $fpkm_low, $fpkm_high, $fpkm_stat ) = split /\t/;
		unless ($track eq "tracking_id"){ #check & specifying undefined variables to null
			if($class =~ /-/){$class = undef;} if ($ref_id =~ /-/){$ref_id = undef;}
			if ($length =~ /-/){$length = undef;} if($coverage =~ /-/){$coverage = undef;}
			$locus =~ /^(.+)\:(.+)\-(.+)$/;
			$chrom_no = $1; $chrom_start = $2; $chrom_stop = $3;
			$sth ->execute($lib_id, $track, $class, $ref_id, $gene, $gene_name, $tss, $chrom_no, $chrom_start, $chrom_stop, $length, $coverage, $fpkm, $fpkm_low, $fpkm_high, $fpkm_stat );
		}
	}
	close GENES;
}

# SUMMARY OF RESULTS
# print "\n\tEXECUTING SELECT STATEMENT ON THE DATABASE TABLES \n";


# $syntax = "select  fpkm, tracking_id, gene_id, chrom_no from GENES_FPKM where fpkm >= 99999.99999 and library_id like \'$lib_id\'";
# $sth = $dbh->prepare($syntax);
# $sth->execute or die "SQL Error: $DBI::errstr\n";
# while (@row = $sth->fetchrow_array() ) {print "@row\n";}

# DISCONNECT FROM THE DATABASE
print "\n\tDISCONNECTING FROM THE DATABASE : $dsn\n\n";
$dbh->disconnect();

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
print "\n\n*********DONE*********\n\n";
# - - - - - - - - - - - - - - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
exit;

