#!/usr/bin/perl

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# MANUAL FOR DBinsert4mapcount.pl

=pod

=head1 NAME

$0 -- Insert results from tophat and cufflinks analysis into the mysql database : PENGUIN

=head1 SYNOPSIS

DBinsert4mapcount.pl [--help] [--manual]

=head1 DESCRIPTION

Accepts the folders containing the output from Allen's frnakenstein tool
 
=head1 OPTIONS

=over 3

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

Copyright 2015 MOA.  
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

use strict; use DBI; use DBD::mysql; use Getopt::Long; use Pod::Usage; 

# CREATING LOG FILES
my $std_out = '/home/modupe/.LOG/output-'.`date +%m-%d-%y_%T`; chomp $std_out; $std_out = $std_out.'.log';
my $std_err = '/home/modupe/.LOG/output-'.`date +%m-%d-%y_%T`; chomp $std_err; $std_err = $std_err.'.err';

open(STDOUT, '>', "$std_out") or die "Log file doesn't exist";
open(STDERR, '>', "$std_err") or die "Error file doesn't exist";


#CREATING EMAIL NOTIFICATION
my $notification = '/home/modupe/.LOG/note.txt';
my $email = 'amodupe@udel.edu';
open (NOTE, ">$notification"); print NOTE "Subject: RAVEN Script Running!!!\n\nName of files\n\t$std_out\n\t$std_err\n"; 
system "sendmail $email < $notification";
close NOTE;
system "rm -rf $notification";

#ARGUMENTS
my($help,$manual);
GetOptions (
                                "h|help"        =>      \$help,
                                "man|manual"	=>      \$manual );

# VALIDATE ARGS
pod2usage( -verbose => 2 )  if ($manual);
pod2usage( -verbose => 1 )  if ($help);

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

#EXCLUDING MICRORNAs
my $microRNA = "778,640,641,642,643,644";
my @excludedmicroRNA = split(",",$microRNA);
foreach my $v (0..$#excludedmicroRNA){
	$Hashresults{$excludedmicroRNA[$v]} = $number; $number++;
}

# DATABASE VARIABLES
my ($dbh, $sth, $syntax, $row, @row);
my $code="DE";
my $inputdirectory='/var/www/subdirectories_for_interface/';


#PARSING VARIABLES
my @parse; my $len;

# TABLE VARIABLES
my ($lib_id, $total, $mapped, $unmapped, $deletions, $insertions, $junctions, $genes, $isoforms, $prep, $date); #RESULTS_SUMMARY TABLE
my ($raw_reads, $fastqc, $accepted, $unmapped_bam, $deletions_bed, $insertions_bed, $junctions_bed, $skipped_gtf, $transcripts_gtf, $isoforms_fpkm, $genes_fpkm, $run_log); #ACTUAL_FILES TABLE
my ($track, $class, $ref_id, $gene, $gene_name, $tss, $locus, $chrom_no, $chrom_start, $chrom_stop, $length, $coverage, $fpkm, $fpkm_low, $fpkm_high, $fpkm_stat); # GENES_FPKM & ISOFORMS_FPKM TABLE

#OPENING mapcount FOLDER '/var/www/subdirectories_for_interface/mapcount_output'
my $mapcount = $inputdirectory.'mapcount_output';
opendir(DIR,$mapcount) or die "Folder \"$mapcount\" doesn't exist\n";
my @Directory = readdir(DIR);
close(DIR);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - M A I N - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# CONNECT TO THE DATABASE
print "\n\n\tCONNECTING TO THE DATABASE : $dsn\n\n";
$dbh = DBI->connect($dsn, $user, $passwd) or die "Connection Error: $DBI::errstr\n";

#CHECKING THE LIBRARIES ALREADY IN THE DATABASE
$syntax = "select library_id from RESULTS_SUMMARY";
$sth = $dbh->prepare($syntax);
$sth->execute or die "SQL Error: $DBI::errstr\n";
while ($row = $sth->fetchrow_array() ) {
        my @hashtablearray = split('', $row);
	pop (@hashtablearray); pop (@hashtablearray);
	my $newrow = join ('',@hashtablearray);
	$Hashresults{$newrow} = $number; $number++;
}
# EACH FOLDER
foreach my $FOLDER (@Directory ){ 
	if ($FOLDER =~ /^\w.*_(\d.*)$/){
        	unless (exists $Hashresults{$1}) {
			$lib_id = "$1$code";
        		print "\n. . .\n\tWORKING ON LIBRARY-$1\n\n";
        		$top_folder = "$mapcount/$FOLDER/tophat_out"; $cuff_folder = "$mapcount/$FOLDER/cufflinks_out";
        		open(ALIGN, "<$top_folder/align_summary.txt") or die "Can't open file $top_folder/align_summary.txt\n";
        		open(GENES, "<$cuff_folder/genes.fpkm_tracking") or die "Can't open file $cuff_folder/genes.fpkm_tracking\n";
        		open(ISOFORMS, "<$cuff_folder/isoforms.fpkm_tracking") or die "Can't open file $cuff_folder/isoforms.fpkm_tracking\n";

# PARSER FOR RESULTS_SUMMARY TABLE
			while (<ALIGN>){
				chomp;
        			if (/Input/){my $line = $_; $line =~ /Input.*:\s+(\d+)$/;$total = $1;}
        			if (/Mapped/){my $line = $_; $line =~ /Mapped.*:\s+(\d+).*$/;$mapped = $1;}
        		}
 			$unmapped = $total-$mapped;
	 		$deletions = `cat $top_folder/deletions.bed | wc -l`; $deletions--;
	 		$insertions = `cat $top_folder/insertions.bed | wc -l`; $insertions--;
 			$junctions = `cat $top_folder/junctions.bed | wc -l`; $junctions--;
 			$genes = `cat $cuff_folder/genes.fpkm_tracking | wc -l`; $genes--;
 			$isoforms = `cat $cuff_folder/isoforms.fpkm_tracking | wc -l`; $isoforms--;
 			$prep = `cat $top_folder/prep_reads.info`;
 			$date = `date +%Y-%m-%d`;

#PARSER FOR RAW READS & FASTQC
                        $raw_reads = `find /home/schmidt_fastq/ -name '*$1*gz'`;
                        my @allraw_reads = split("\n", $raw_reads);
                        if (length $raw_reads < 1){print "FASTQ file isn't present\t==>\tSkipping library_$1!!\n\n"; next;} # 4 ACTUAL_FILES table

                        print "\tRunning FASTQC on library_$1\n\n"; # 4 ACTUAL_FILES table
                 #RUNNING FASTQC 
                        `mkdir $mapcount/$FOLDER/fastqc_out/`;
                        my $fastqc_folder = "$mapcount/$FOLDER/fastqc_out";
			if ($#allraw_reads > 0){
                        	foreach my $i (0..$#allraw_reads){
                                	`/home/tjcarr/FastQC/fastqc -t 4 -o $fastqc_folder $allraw_reads[$i]`
                          	}
                        }
                        else {
               	                `/home/tjcarr/FastQC/fastqc -t 4 -o $fastqc_folder $raw_reads`;
                        }
                      
                        #geting FASTQC output zip file
                        $fastqc = `ls $mapcount/$FOLDER/fastqc_out/*zip`;
                        if (length $fastqc < 1){print "\tFASTQC file doesn't exist\t==>\tSkipping library_$1!!\n\n"; next;}# 4 ACTUAL_FILES table

# PARSER FOR ACTUAL_FILES TABLE	
			$accepted = "$top_folder/accepted_hits.bam"; @parse = split('\/\/',$accepted); $accepted = undef; $len = $#parse+1; foreach(@parse){$accepted .= $_; if($len>1){$accepted .="\/"; $len--;}};
			$unmapped_bam = "$top_folder/unmapped.bam"; @parse = split('\/\/',$unmapped_bam); $unmapped_bam = undef; $len = $#parse+1; foreach(@parse){$unmapped_bam .= $_; if($len>1){$unmapped_bam .="\/"; $len--;}};
			$deletions_bed = "$top_folder/deletions.bed"; @parse = split('\/\/',$deletions_bed); $deletions_bed = undef; $len = $#parse+1; foreach(@parse){$deletions_bed .= $_; if($len>1){$deletions_bed .="\/"; $len--;}};
			$insertions_bed = "$top_folder/insertions.bed"; @parse = split('\/\/',$insertions_bed); $insertions_bed = undef; $len = $#parse+1; foreach(@parse){$insertions_bed .= $_; if($len>1){$insertions_bed .="\/"; $len--;}};
			$junctions_bed = "$top_folder/junctions.bed"; @parse = split('\/\/',$junctions_bed); $junctions_bed = undef; $len = $#parse+1; foreach(@parse){$junctions_bed .= $_; if($len>1){$junctions_bed .="\/"; $len--;}};
			$skipped_gtf = "$cuff_folder/skipped.gtf"; @parse = split('\/\/',$skipped_gtf); $skipped_gtf = undef; $len = $#parse+1; foreach(@parse){$skipped_gtf .= $_; if($len>1){$skipped_gtf .="\/"; $len--;}};
			$transcripts_gtf = "$cuff_folder/transcripts_gtf"; @parse = split('\/\/',$transcripts_gtf); $transcripts_gtf = undef; $len = $#parse+1; foreach(@parse){$transcripts_gtf .= $_; if($len>1){$transcripts_gtf .="\/"; $len--;}};
			$isoforms_fpkm = "$cuff_folder/isoforms.fpkm_tracking"; @parse = split('\/\/',$isoforms_fpkm); $isoforms_fpkm = undef; $len = $#parse+1; foreach(@parse){$isoforms_fpkm .= $_; if($len>1){$isoforms_fpkm .="\/"; $len--;}};
			$genes_fpkm = "$cuff_folder/genes.fpkm_tracking"; @parse = split('\/\/',$genes_fpkm); $genes_fpkm = undef; $len = $#parse+1; foreach(@parse){$genes_fpkm .= $_; if($len>1){$genes_fpkm .="\/"; $len--;}};
			$run_log = "$top_folder/logs/run.log"; my @parse = split('\/\/',$run_log); $run_log = undef; $len = $#parse+1; foreach(@parse){$run_log .= $_; if($len>1){$run_log .="\/"; $len--;}};

# INSERT INTO DATABASE : PENGUIN
	        	#RESULTS_SUMMARY table
    			print ". . .\n\tINSERTING INTO THE DATABASE IN \"RESULTS_SUMMARY\" TABLE\n\n";
        		$sth = $dbh->prepare("insert into RESULTS_SUMMARY (library_id, total_reads, mapped_reads, unmapped_reads, deletions, insertions, junctions, isoforms, genes, info_prep_reads, date ) values (?,?,?,?,?,?,?,?,?,?,?)");
        		$sth ->execute($lib_id, $total, $mapped, $unmapped, $deletions, $insertions, $junctions, $isoforms, $genes, $prep, $date);

        		#ACTUAL_FILES table
			print ". . .\n\tINSERTING INTO THE DATABASE IN \"ACTUAL_FILES\" TABLE\n\n";
        	       	$sth = $dbh->prepare("insert into ACTUAL_FILES (library_id, RAW_READS, FASTQC_html, ACCEPTED_HITS_bam, UNMAPPED_bam, DELETIONS_bed, INSERTIONS_bed, JUNCTIONS_bed, SKIPPED_gtf, TRANSCRIPTS_gtf, ISOFORMS_FPKM, GENES_FPKM, RUN_LOG ) values (?,?,?,?,?,?,?,?,?,?,?,?,?)");
             		$sth ->execute($lib_id, $raw_reads, $fastqc, $accepted, $unmapped_bam, $deletions_bed, $insertions_bed, $junctions_bed, $skipped_gtf, $transcripts_gtf, $isoforms_fpkm, $genes_fpkm, $run_log);

			#GENES_FPKM table
	        	print ". . .\n\tINSERTING INTO THE DATABASE IN \"GENES_FPKM\" TABLE\n\n";
        		$sth = $dbh->prepare("insert into GENES_FPKM (library_id, tracking_id, class_code, nearest_ref_id, gene_id, gene_short_name, tss_id, chrom_no, chrom_start, chrom_stop, length, coverage, fpkm, fpkm_conf_low, fpkm_conf_high, fpkm_status ) values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)");
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
        		} close GENES;

        		#ISOFORMS_FPKM table
        		print "...\n\tINSERTING INTO THE DATABASE IN \"ISOFORMS_FPKM\" TABLE\n\n";
        		$sth = $dbh->prepare("insert into ISOFORMS_FPKM (library_id, tracking_id, class_code, nearest_ref_id, gene_id, gene_short_name, tss_id, chrom_no, chrom_start, chrom_stop, length, coverage, fpkm, fpkm_conf_low, fpkm_conf_high, fpkm_status ) values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)");
        		while (<ISOFORMS>){
	                	chomp;
        		       	my ($track, $class, $ref_id, $gene, $gene_name, $tss, $locus, $length, $coverage, $fpkm, $fpkm_low, $fpkm_high, $fpkm_stat ) = split /\t/;
                		unless ($track eq "tracking_id"){
                        	        if ($class =~ /-/){$class = undef;} if ($ref_id =~ /-/){$ref_id = undef;}
                               		if ($length =~ /-/){$length = undef;} if($coverage =~ /-/){$coverage = undef;}
                      	       	  	$locus =~ /^(.+)\:(.+)\-(.+)$/;
                        	        $chrom_no = $1; $chrom_start = $2; $chrom_stop = $3;
                               		$sth ->execute($lib_id, $track, $class, $ref_id, $gene, $gene_name, $tss, $chrom_no, $chrom_start, $chrom_stop, $length, $coverage, $fpkm, $fpkm_low, $fpkm_high, $fpkm_stat );
                     		}
        		} close ISOFORMS;
	                #RESULTS_SUMMARY table updating status column with 'done'
        	        $sth = $dbh->prepare("update RESULTS_SUMMARY set status='done' where library_id = ?");
                	$sth ->execute( $lib_id );
        	}
	}	 
}

# DISCONNECT FROM THE DATABASE
print "\n\tDISCONNECTING FROM THE DATABASE : $dsn\n\n";
$dbh->disconnect();

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
print "\n\n*********DONE*********\n\n";
# - - - - - - - - - - - - - - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
exit;
