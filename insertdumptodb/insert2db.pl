#!/usr/bin/perl
use DBI; 
use DBD::mysql; 
use Getopt::Long; 
use Pod::Usage;
use threads;
use Thread::Queue;
use lib '/home/modupe/SCRIPTS/SUB';
use passw;
use routine;

my ($dbh, $sth, $syntax, @VAR); my %HashDirectory;
$dbh = mysql();

open(IN, $ARGV[0]) or die "Input file is not specified\n";
#my $syntax = "insert into genes_fpkm (library_id, tracking_id, class_code, nearest_ref_id, gene_id, gene_short_name, tss_id, chrom_no, chrom_start, chrom_stop, length, coverage, fpkm, fpkm_conf_low, fpkm_conf_high, fpkm_status ) values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)";
my $syntax = "insert into genes_fpkm (library_id, gene_id, gene_short_name, chrom_no, chrom_start, chrom_stop, coverage, fpkm, fpkm_conf_low, fpkm_conf_high, fpkm_status ) values (?,?,?,?,?,?,?,?,?,?,?)";
$sth = $dbh->prepare($syntax);
my @details = <IN>;
shift @details;
push @VAR, [ splice @details, 0, 200 ] while @details; #sub array the genes into a list of 2000
	
#foreach (0..$#VAR){ $newfile .= "tadtmp/tmp_".$tmpname."-".$_.".zzz "; } #foreach sub array create a temporary file
my $queue = new Thread::Queue();
my $builder=threads->create(\&main); #create thread for each subarray into a thread
push @threads, threads->create(\&processor) for 1..5; #execute 5 threads
$builder->join; #join threads
foreach (@threads){$_->join;}
sub main {
        foreach my $count (0..$#VAR) {
		#my $namefile = "tadtmp/tmp_".$tmpname."-".$count.".zzz";
		#push $VAR[$count], $namefile;
		while(1) {
			if ($queue->pending() <100) {
				$queue->enqueue($VAR[$count]);
				last;
			}
		}
	}
	foreach(1..5) { $queue-> enqueue(undef); }
}

sub processor {
	my $query;
	while ($query = $queue->dequeue()){
		parseinput(@$query);
	}
}

sub parseinput{
        $i = 0;
	#my $file = pop @_;
	#open(OUT2, ">$file");
	foreach (@_){
                chomp; print $i++,"\n";
                $dbh = mysql();
                $syntax = "insert into genes_fpkm (library_id, gene_id, gene_short_name, chrom_no, chrom_start, chrom_stop, coverage, fpkm, fpkm_conf_low, fpkm_conf_high, fpkm_status ) values (?,?,?,?,?,?,?,?,?,?,?)";
                $sth = $dbh->prepare($syntax);
                my ($lib, $track, $class, $ref_id, $gene, $gene_name, $tss, $chrom, $start, $stop, $length, $coverage, $fpkm, $fpkm_low, $fpkm_high, $fpkm_stat ) = split /\t/;
                $sth ->execute($lib, $gene, $gene_name, $chrom, $start, $stop, $coverage, $fpkm, $fpkm_low, $fpkm_high, $fpkm_stat) or die "\nERROR:\t Complication in genes_fpkm table, consult documentation\n";
                $sth->finish;
		#sortposition($_);
	}
}

#foreach (@details){
#	#print $_;
#        chomp;
#	my ($lib, $track, $class, $ref_id, $gene, $gene_name, $tss, $chrom, $start, $stop, $length, $coverage, $fpkm, $fpkm_low, $fpkm_high, $fpkm_stat ) = split /\t/;
#        $sth ->execute($lib, $gene, $gene_name, $chrom, $start, $stop, $coverage, $fpkm, $fpkm_low, $fpkm_high, $fpkm_stat) or die "\nERROR:\t Complication in genes_fpkm table, consult documentation\n";
#}
