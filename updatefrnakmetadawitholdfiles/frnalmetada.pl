#!/usr/bin/perl
#convert the different files in a provided folder (standard input) to a tab-delimited output for import to frnak_metadaa table
use strict;
use Data::Dumper;

opendir(DIR, $ARGV[0]);
open (OUT, ">files.txt");
print OUT "library_id\tref_genome\tann_file\tann_file_ver\tstranded\tsequences\tuser\n";
my @files = readdir(DIR); close (DIR);
my %ALL;
`mkdir -p FINAL`;
foreach my $libs (@files){
  undef %ALL;
  if ($libs =~ /^([0-9].*)\.txt$/){ 
    my $lib_id = $1; 
    my $logfile = `head -n 1 $libs`;
    my %parsablegenomes = ("chicken" => 1, "alligator" => 2,"mouse" => 3, ); #genomes that work.
    my @allgeninfo = split('\s',$logfile);
  
    #also getting metadata info
    my ($str, $ann, $ref, $seq,$otherseq,$allstart, $allend) = 0; #defining calling variables
    #making sure I'm working on only the chicken files for now, need to find annotation of alligator
    my @checkerno = split('\/',$allgeninfo[$#allgeninfo]);
    my @numberz =split('_', $checkerno[$#checkerno]);
    my $no = 0;
    my ($refgenome, $temptest,$stranded, $sequences, $annotation, $annotationfile, $annfileversion);
    $numberz[0] =~ /^s([0-9].*)$/;
    if ($1 == $lib_id){ 
      my $number = 1;
      while ($number <= $#allgeninfo){
        unless ($allgeninfo[$number] =~ /-no-coverage-search/){
          if ($allgeninfo[$number] =~ /^\-/){
            my $old = $number++;
            $ALL{$allgeninfo[$old]} = $allgeninfo[$number];
          } else {
            unless (exists $ALL{$no}){
              $ALL{$no} = $allgeninfo[$number];
              $no++;
            }
          }
        }
        $number++;
      }
      unless (exists $ALL{"-G"}) {
        $annotation = undef;
        $annotationfile = undef;
        $annfileversion = undef;
      } else {
        $annotation = $ALL{"-G"};
        $annotationfile = uc ( (split('\.',((split("\/", $annotation))[-1])))[-1] ); #(annotation file)
        $annfileversion = substr(`head -n 1 $annotation`,2,-1); #annotation file version
      }
      unless (exists $ALL{"-library-type"}) { $stranded = undef; } else { $stranded = $ALL{"-library-type"}; }
      $refgenome = $ALL{0}; $seq = $ALL{1}; $otherseq = $ALL{2};
      #print Data::Dumper->Dump( [ \%ALL ], [ qw(*GENE) ] ); die;
      my ($one, $two);
      $one = ( ( split('\/', $seq) ) [-1]);
      unless ($one =~ /[fastq|gz|fq]$/){
        $one .= ".";
        $one = `locate $one | head -n 1`;
        chop $one;
      }
      unless(length($otherseq)<1){
        $two = ( ( split('\/', $otherseq) ) [-1]);
        unless ($two =~ /[fastq|gz|fq]$/){
          $two .= ".";
          $two = `locate $one | head -n 1`;
          chop $two;
        }
      }
      unless(length($otherseq)<1){ #sequences
        $sequences = ( ( split('\/', $one) ) [-1]).",". ( ( split('\/', $two) ) [-1]);
      } else {
        $sequences = ( ( split('\/', $one) ) [-1]);
      }
      unless ($sequences =~ /[fastq|gz|fq]$/){ print Data::Dumper->Dump( [ \%ALL ], [ qw(*INDEX) ] ); print "Sequence error: \'$lib_id \t $one \t $two \t$sequences\'\n"; next;}
      if (exists $parsablegenomes{$refgenome} || $refgenome =~ /Galgal4/){
        if ($refgenome =~ /Galgal4/){$refgenome = "chicken";}
        print OUT "$lib_id\t$refgenome\t$annotationfile\t$annfileversion\t$stranded\t$sequences\told update\n";
        `mv $libs FINAL/`;
      }
    }
  }
}
