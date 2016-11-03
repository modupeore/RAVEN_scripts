#!/usr/bin/perl
use strict;
use File::Basename;

#VARIABLES
my (%GT, %TISSUE, %REF, %ALT, %QUAL, %CSQ, %DBSNP);
my (%ODACSQ,%number);
my (%NEWQUAL, %NEWCSQ, %NEWREF, %NEWDBSNP, %NEWALT,%NEWGT);
my (%subref, %subalt, %subgt);
my %MTD;

#PERFORM SUBROUTINES
PROCESS($ARGV[0]);
SORTER();
MTD();
our $headerinfo = HEADER();

#GET OUTPUT FILE NAME
my $out = fileparse($ARGV[0], qr/\.[^.]*(\.txt)?$/);
my $output = $out.".vcf";

#PRINT OUTPUT
open (OUT, ">", $output) or die "Can't open vcf file --weird\n";
print OUT $headerinfo;
foreach my $chrom (sort {$a <=> $b} keys %NEWREF) {
  foreach my $position (sort {$a<=> $b} keys %{$NEWREF{$chrom}}) {
    foreach my $ref (sort {$a<=> $b} keys %{$NEWREF{$chrom}{$position}}) {
      print OUT "chr",$chrom,"\t",$position,"\t",$NEWDBSNP{$chrom}{$position}{$ref},"\t",$NEWREF{$chrom}{$position}{$ref},"\t";
      print OUT $NEWALT{$chrom}{$position}{$ref},"\t",$NEWQUAL{$chrom}{$position}{$ref},"\tPASS\tCSQ=",$NEWCSQ{$chrom}{$position}{$ref};
      print OUT ";MTD=",$MTD{$chrom}{$position}{$ref},"\tGT\t",$NEWGT{$chrom}{$position}{$ref};
      print OUT "\n";
    }
  }
}
close (OUT);

#--------------------------------------------------------------------------------------------------------
sub HEADER {
#header information
  $headerinfo = <<'ENDOFFILE';
##fileformat=VCFv4.1
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##reference=file:///.GENOMES/chicken/chicken.fa
##VEP=v81 genebuild=2013-04 sift=sift5.2.2 dbSNP=140 assembly=Galgal4
##INFO=<ID=MTD,Number=.,Type=String,Description="Metadata information from Schmidt Lab. Format:Library|Tissue|Quality|Genotype">
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format:Consequence|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|Protein_position|Amino_acids|Codons|Existing_variation|VARIANT_CLASS">
#CHROM    POS    ID    REF    ALT    QUAL    FILTER    INFO    FORMAT    Label
ENDOFFILE
  return $headerinfo;
}

sub PROCESS {
#processt the input file
  open (IN, $_[0]);
  while (<IN>){
    chomp;
    my @line = split("\t",$_,19);
    if ($line[0] =~ /^\d/){
      my $joint = "$line[6]|$line[7]|$line[8]|$line[9]|$line[10]|$line[11]|$line[12]|$line[13]|$line[14]|$line[15]|$line[16]";
      $line[1] = substr($line[1],3);
      $TISSUE{$line[0]} = lc($line[18]);
      $REF{$line[1]}{$line[2]}{$line[3]}{$line[0]} = $line[3];
      $ALT{$line[1]}{$line[2]}{$line[3]}{$line[0]} = $line[4];
      $QUAL{$line[1]}{$line[2]}{$line[3]}{$line[0]} = $line[5];
      if (exists $CSQ{$line[1]}{$line[2]}{$line[3]}{$line[0]}) {
        $number{$line[1]}{$line[2]}{$line[3]}{$line[0]} = $number{$line[1]}{$line[2]}{$line[3]}{$line[0]}++;
        $ODACSQ{$line[1]}{$line[2]}{$line[3]}{$line[0]}{$number{$line[1]}{$line[2]}{$line[3]}{$line[0]}} = $joint;
        $CSQ{$line[1]}{$line[2]}{$line[3]}{$line[0]} =  "$CSQ{$line[1]}{$line[2]}{$line[3]}{$line[0]},$joint";
      }
      else {
        $number{$line[1]}{$line[2]}{$line[3]}{$line[0]}= 1;
        $ODACSQ{$line[1]}{$line[2]}{$line[3]}{$line[0]}{1} = $joint;
        $CSQ{$line[1]}{$line[2]}{$line[3]}{$line[0]} =  $joint;
      }
      my $verdict = undef;
      if ($line[17] =~ /^homozygous/){
        $verdict = '1/1';
      }
      elsif ($line[17] =~ /alternate/){
        $verdict = '1/2';
      }
      elsif ($line[17] =~ /^heterozygous$/){
        $verdict = '0/1';
      }
      else {die "zygosity is blank\n";}
      if (exists $GT{$line[1]}{$line[2]}{$line[3]}{$line[0]}) {
        unless ($GT{$line[1]}{$line[2]}{$line[3]}{$line[0]} =~ $verdict){
          die "Genotype information is different: $verdict is not ",$GT{$line[1]}{$line[2]}{$line[3]}{$line[0]},"\n";
        }
      }
      else {
        $GT{$line[1]}{$line[2]}{$line[3]}{$line[0]} = $verdict;
      }
      if (length($line[15])<1){$line[15]='.';}
      if (exists $DBSNP{$line[1]}{$line[2]}{$line[3]}{$line[0]}) {
        unless ($DBSNP{$line[1]}{$line[2]}{$line[3]}{$line[0]} =~ $line[15]){
          die "DBSNPs information is different: $line[15] is not ",$DBSNP{$line[1]}{$line[2]}{$line[3]}{$line[0]},"\n";
        }
      }
      else {
        $DBSNP{$line[1]}{$line[2]}{$line[3]}{$line[0]} = $line[15];
      }
    }
  }
  close(IN);
}

sub SORTER {
  #SORT ALLELES
  foreach my $chrom (sort {$a <=> $b} keys %REF) {
    foreach my $position (sort {$a <=> $b} keys %{$REF{$chrom}}) {
      foreach my $ref (sort {$a <=> $b} keys %{$REF{$chrom}{$position}}) {
        foreach my $library (sort {$a <=> $b} keys %{$REF{$chrom}{$position}{$ref}}) {
          if (exists $subref{$chrom}{$position}{$ref}){
            unless ($subref{$chrom}{$position}{$ref} =~ $REF{$chrom}{$position}{$ref}{$library}){
              $subref{$chrom}{$position}{$ref}= $subref{$chrom}{$position}{$ref}.",".$REF{$chrom}{$position}{$ref}{$library};
            }
          }
          else {
            $subref{$chrom}{$position}{$ref}= $REF{$chrom}{$position}{$ref}{$library};
          }
        }
        foreach my $library (sort {$a<=> $b} keys %{$ALT{$chrom}{$position}{$ref}}) {
          if (exists $subalt{$chrom}{$position}{$ref}){
            unless ($subalt{$chrom}{$position}{$ref} =~ $ALT{$chrom}{$position}{$ref}{$library}){
              $subalt{$chrom}{$position}{$ref} = $subalt{$chrom}{$position}{$ref}.",".$ALT{$chrom}{$position}{$ref}{$library};
            }
          }
          else {
            $subalt{$chrom}{$position}{$ref}= $ALT{$chrom}{$position}{$ref}{$library};
          }
        }
      }
    }
  }
  
  #sub sort REF & ALT alleles
  foreach my $chrom (sort {$a <=> $b} keys %subref) {
    foreach my $position (sort {$a<=> $b} keys %{$subref{$chrom}}) {
      foreach my $ref (sort {$a<=> $b} keys %{$subref{$chrom}{$position}}) {
        my %refhash = ''; my %althash = ''; my $refkey = undef; my $altkey = undef;
        my @refarray = split(",", $subref{$chrom}{$position}{$ref});
        foreach (sort {$a cmp $b} @refarray) {$refhash{$_} = $_;}
        foreach (sort {$a cmp $b} keys %refhash){ $refkey .= $_.","; } 
        $NEWREF{$chrom}{$position}{$ref} = substr ($refkey, 1, -1); 
        
        my @altarray = split(",", $subalt{$chrom}{$position}{$ref});
        foreach (sort {$a cmp $b} @altarray) {$althash{$_} = $_;}
        foreach (sort {$a cmp $b} keys %althash){ $altkey .= $_.","; }
        $NEWALT{$chrom}{$position}{$ref} = substr ($altkey, 1,-1);
      }
    }
  }
  
  #SORT CONSEQUENCE
  foreach my $chrom (sort {$a <=> $b} keys %CSQ) {
    foreach my $position (sort {$a<=> $b} keys %{$CSQ{$chrom}}) {
      foreach my $ref (sort {$a<=> $b} keys %{$CSQ{$chrom}{$position}}) {
        foreach my $library (sort {$a <=> $b} keys %{$CSQ{$chrom}{$position}{$ref}}) {
          if (exists $NEWCSQ{$chrom}{$position}{$ref}){
            unless ($NEWCSQ{$chrom}{$position}{$ref} =~ $CSQ{$chrom}{$position}{$ref}{$library}){
              die "Consequence should be the same:  $chrom $position not $NEWCSQ{$chrom}{$position}{$ref} equals $CSQ{$chrom}{$position}{$ref}{$library}\n";
            }
          }  
          else {
            $NEWCSQ{$chrom}{$position}{$ref} = $CSQ{$chrom}{$position}{$ref}{$library};
          }
        }
      }
    }
  }
  
  #SORT QUALITY
  foreach my $chrom (sort {$a <=> $b} keys %QUAL) {
    foreach my $position (sort {$a<=> $b} keys %{$QUAL{$chrom}}) {
      foreach my $ref (sort {$a<=> $b} keys %{$QUAL{$chrom}{$position}}) {
        my @quality = undef;
        foreach my $library (sort {$a<=> $b} keys %{$QUAL{$chrom}{$position}{$ref}}) {
          push @quality, $QUAL{$chrom}{$position}{$ref}{$library};
        }
        @quality = sort {$a <=> $b} @quality;
        $NEWQUAL{$chrom}{$position}{$ref} = $quality[$#quality];
      }
    }
  }
  
  #SORT DBSNP
  foreach my $chrom (sort {$a <=> $b} keys %DBSNP) {
    foreach my $position (sort {$a<=> $b} keys %{$DBSNP{$chrom}}) {
      foreach my $ref (sort {$a<=> $b} keys %{$DBSNP{$chrom}{$position}}) {
        foreach my $library (sort {$a <=> $b} keys %{$DBSNP{$chrom}{$position}{$ref}}) {
          if (exists $NEWDBSNP{$chrom}{$position}{$ref}){
            unless ($NEWDBSNP{$chrom}{$position}{$ref} =~ $DBSNP{$chrom}{$position}{$ref}{$library}){
              $NEWDBSNP{$chrom}{$position}{$ref} = '.';
            }
          }  
          else {
            $NEWDBSNP{$chrom}{$position}{$ref} = $DBSNP{$chrom}{$position}{$ref}{$library};
          }
        }
      }
    }
  }
  
  #SORT GENOTYPE
  foreach my $chrom (sort {$a <=> $b} keys %GT) {
    foreach my $position (sort {$a<=> $b} keys %{$GT{$chrom}}) {
      foreach my $ref (sort {$a<=> $b} keys %{$GT{$chrom}{$position}}) {
        foreach my $library (sort {$a <=> $b} keys %{$GT{$chrom}{$position}{$ref}}) {
          if (exists $NEWGT{$chrom}{$position}{$ref}){
            unless ($NEWGT{$chrom}{$position}{$ref} =~ $GT{$chrom}{$position}{$ref}{$library}){
              $subgt{$chrom}{$position}{$ref}{$GT{$chrom}{$position}{$ref}{$library}}++;
              #die "Genotype should be the same:  $chrom $position not $NEWGT{$chrom}{$position}{$ref} equals $GT{$chrom}{$position}{$ref}{$library}\n";
            }
          }  
          else {
            $subgt{$chrom}{$position}{$ref}{$GT{$chrom}{$position}{$ref}{$library}} = 1;
            $NEWGT{$chrom}{$position}{$ref} = $GT{$chrom}{$position}{$ref}{$library};
          }
        }
      }
    }
  }
  
  #order genotype
  my %odagt;
  foreach my $chrom (sort {$a <=> $b} keys %subgt) {
    foreach my $position (sort {$a<=> $b} keys %{$subgt{$chrom}}) {
      foreach my $ref (sort {$a<=> $b} keys %{$subgt{$chrom}{$position}}) {
        if ( (exists $subgt{$chrom}{$position}{$ref}{'0/1'}) && (exists $subgt{$chrom}{$position}{$ref}{'1/2'}) ){
          if ( $subgt{$chrom}{$position}{$ref}{'0/1'} > $subgt{$chrom}{$position}{$ref}{'1/2'} ) {
            $subgt{$chrom}{$position}{$ref}{'0/1'} =  $subgt{$chrom}{$position}{$ref}{'0/1'} + $subgt{$chrom}{$position}{$ref}{'1/2'};
          }
          if ( $subgt{$chrom}{$position}{$ref}{'0/1'} < $subgt{$chrom}{$position}{$ref}{'1/2'} ) {
            $subgt{$chrom}{$position}{$ref}{'1/2'} =  $subgt{$chrom}{$position}{$ref}{'0/1'} + $subgt{$chrom}{$position}{$ref}{'1/2'};
          }
          if ( $subgt{$chrom}{$position}{$ref}{'0/1'} == $subgt{$chrom}{$position}{$ref}{'1/2'} ) {
            $subgt{$chrom}{$position}{$ref}{'0/1'} =  $subgt{$chrom}{$position}{$ref}{'0/1'} + $subgt{$chrom}{$position}{$ref}{'1/2'};
          }
        }
        foreach my $geno (sort {$a <=> $b} keys %{$subgt{$chrom}{$position}{$ref}}){
          $odagt{$chrom}{$position}{$ref}{$subgt{$chrom}{$position}{$ref}{$geno}} = $geno;
        }
      }
    }
  }
  foreach my $chrom (sort {$a <=> $b} keys %odagt) {
    foreach my $position (sort {$a<=> $b} keys %{$odagt{$chrom}}) {
      foreach my $ref (sort {$a<=> $b} keys %{$odagt{$chrom}{$position}}) {
        my $newpost = (sort {$a <=> $b} keys %{$odagt{$chrom}{$position}{$ref}})[0];
        $NEWGT{$chrom}{$position}{$ref} = $odagt{$chrom}{$position}{$ref}{$newpost};
      }
    }
  }
}

sub MTD {
  #get metadata information
  foreach my $chrom (sort {$a <=> $b} keys %QUAL) {
    foreach my $position (sort {$a<=> $b} keys %{$QUAL{$chrom}}) {
      foreach my $ref (sort {$a<=> $b} keys %{$QUAL{$chrom}{$position}}) {
        foreach my $library (sort {$a<=> $b} keys %{$QUAL{$chrom}{$position}{$ref}}) {
          if (exists $MTD{$chrom}{$position}{$ref}) {
            $MTD{$chrom}{$position}{$ref} = $MTD{$chrom}{$position}{$ref}.",$library|$TISSUE{$library}|$QUAL{$chrom}{$position}{$ref}{$library}|$GT{$chrom}{$position}{$ref}{$library}";
          }
          else {
            $MTD{$chrom}{$position}{$ref} = "$library|$TISSUE{$library}|$QUAL{$chrom}{$position}{$ref}{$library}|$GT{$chrom}{$position}{$ref}{$library}";
          }
        }
      }
    }
  }
}