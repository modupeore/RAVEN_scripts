#!/usr/bin/perl
use strict;
#extracting genes names and chromosomal location from GTF file into a json file.
#this is for genome browser view for the variants-genename webpage.
 
my $gtf = $ARGV[0];
my (%GENE,%POSstart, %POSend, %CHR);
my ($count,$sumofgenes) = (0,0);

GTF_FILE();

open(OUT,">jsongtfgenes.json") or die "can't open the output file";
print OUT "{\n";
foreach my $genename (sort keys %GENE){
  if ($genename =~ /\S+/){
    print OUT "\t\"$genename\" : [";
    my $counted = 0;
    #foreach my $bbb (sort {$a <=> $b} keys %{ $GINFO{$genename} }){
      #$counted++;
      $count++;
      #my @details = split("\t", $GINFO{$genename}{$bbb});
      print OUT "\n\t\t{\"chr\" : \"chr","$CHR{$genename}\",\"start\" : \"$POSstart{$genename}\",";
      print OUT "\"stop\" : \"$POSend{$genename}\" }";
    #  if ($counted < $Gno{$genename}) {
    #    print OUT ",";
    #  } 
    #}
    print OUT "\n\t]";
    if ($count < $sumofgenes) {
      print OUT ",\n";
    }
  }
}
print OUT "\n}\n";
print "a=>$sumofgenes\tb=>$count\n";
close (OUT);
#subroutine


sub GTF_FILE {
  open(GTF,"<",$gtf) or die "$gtf can't open";
  while (<GTF>) {
    chomp;
    my @all = split("\t", $_);
    if ($all[0] =~ /^\d/ || $all[0] =~ /^[L|M|Z|W]/) {
      #if ($#all > 8 || $#all < 8 ) { print $#all,"\n";}
      my @newall = split("\;", $all[8]);
      if($all[2] =~ /gene/){
        foreach my $abc (@newall){
          if ($abc =~ /gene_name/){
            my @bcd = split(" ",$abc,2);
          
            $bcd[1] = substr($bcd[1],1,-1);
            if (exists $GENE{$bcd[1]}){
              #$sumofgenes--;
              print "Already exists $bcd[1]\t$GENE{$bcd[1]}\n";
              #$GENE{$bcd[1]}=$GENE{$bcd[1]}+1;
            }
            else {
              $sumofgenes++;
              $GENE{$bcd[1]}=1;
            }
            $POSstart{$bcd[1]} = $all[3];
            $POSend{$bcd[1]} = $all[4];
            $CHR{$bcd[1]} = $all[0];
          }
        }
      }
    }
  }
  close (GTF);
}

