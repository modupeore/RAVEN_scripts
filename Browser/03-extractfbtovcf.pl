#!/usr/bin/perl
use strict;
use DBI;
use File::Basename;
use Sort::Key::Natural qw(natsort);
use lib '/home/modupe/SCRIPTS/SUB';
use passw;
use routine;


#This is to extract variants and all the pertinent information into a txt file based on the Chicken line
#Output is stored in the folder specified

# DATABASE ATTRIBUTES
my $basepath = int(rand(10000));
`mkdir $basepath`;
`rm -rf $basepath/*`;
my $statusfile = "$basepath/mystatus.txt";

open(STATUS,">>$statusfile") or die "$basepath does not exist\n"; close (STATUS);
# - - - - - - - - - - - - G L O B A L  V A R I A B L E S- - - - - - - - - - - - - -
my ($dbh, $sth);
our ($chickenpath, $mousepath, $alligatorpath) = FBPATHS();

#EXTRACT FROM FB VARIABLES
my (%HashDirectory, @thelines, %LineHash);
my ($ibis, $syntax);

#This is convert the extracted tab-delimited text file into a vcf file to be visualized on the genome browser
#VARIABLES
my (%GT, %TISSUE, %REF, %ALT, %QUAL, %CSQ, %DBSNP);
my (%NEWQUAL, %NEWCSQ, %NEWREF, %NEWDBSNP, %NEWALT, %NEWGT);
my (%subref, %subalt, %subgt);
my ($headerinfo, %MTD);

# - - - - - - - - - - - - - - - - M A I N - - - - - - - - - - - - - - - - - - - - - -

#GETTING ALL THE LIBRARIES FROM THE DATABASE.
@thelines = undef;
$ibis = "/home/modupe/.bin/bin/ibis -d $chickenpath -q '";
$syntax = "select line,countdistinct(line) where line != \"NULL\"' -v -o temp.txt";

print "$ibis$syntax";
`$ibis$syntax`; #execute ibis
print "\nDone with ibis " ,scalar (localtime), "\n";
open(IN, "<temp.txt");
while (<IN>){
	chomp;
	my $linename = (split /,/)[0];
	$linename =~ s/^"|"$//g;
	my $newlinename = $linename; $newlinename =~ s/ /_/g;$newlinename =~ s/&//g;
	$LineHash{$linename} = $newlinename;
}
close (IN);

print "\nDone with parse dictionary " ,scalar (localtime), "\n";

`rm -rf temp.txt`; #removing temp file

#PARSING EACH FILE OUT OF THE DATABASE. 
foreach my $sublibraries (keys %LineHash){
  opendir (DIR, $basepath) or die "Folder \"$basepath\" doesn't exist\n"; 
  my @Directory = readdir(DIR); close(DIR);
  #pushing each subfolder
#	print $basepath;
  foreach (@Directory) { #print "\t\t\t$_\n";
    if ($_ =~ /^(\w.*)\.txt$/) {
      $HashDirectory{$1}= $1;
    }
  }
  unless (exists $HashDirectory{$LineHash{$sublibraries}}){
		if (length($sublibraries) > 1) {
			open(STATUS, ">>$statusfile");
			print STATUS "Processing: $sublibraries\tStarting: ",scalar (localtime),"\n";
			close(STATUS);
			open(OUT,">$basepath/$LineHash{$sublibraries}\.txt") or die "cant open output\n";
			$syntax= "select library,chrom,position,ref,alt,quality,consequence,genename,geneid,feature,transcript,genetype,proteinposition,aachange,codon,dbsnp,class,zygosity,tissue where line = \"$sublibraries\" order by chrom,position' -o newtemp.txt";
	
			`$ibis$syntax`; #running ibis
			
			open(STATUS, ">>$statusfile");
			print STATUS "Midway: $sublibraries\tIBIS complete: ",scalar (localtime),"\n";
			close (STATUS);
    
			print OUT "library_id\tchrom\tposition\tref_allele\talt_allele\tquality\tconsequence\tgene_name\tgene_id\tfeature\ttranscript\tgene_type\tprotein_position\taminoacid_change\tcodon_change\texisting_variant\tvariant_class\tzygosity\ttissue\n";
    
			my $i = 0;
		
			open(IN2, "<newtemp.txt");
			while (<IN2>){
				chomp;
				$i++;
			
				my ($library, $chrom, $position,$ref, $alt, $quality,$consequence, $genename, $geneid,$feature, $transcript, $genetype,$pposition, $aachange, $codon,$dbsnp, $class, $zygosity,$tissue) = split(/\, /, $_, 19);
				#removing the quotation marks from the words
				$chrom =~ s/"//g; $ref =~ s/"//g; $alt =~ s/"//g;
				$consequence =~ s/"//g; $genename =~ s/"//g; $geneid =~ s/"//g;
				$feature =~ s/"//g; $transcript =~ s/"//g; $genetype =~ s/"//g;
				$aachange =~ s/"//g; $codon =~ s/"//g; $dbsnp =~ s/"//g;
				$class =~ s/"//g; $zygosity =~ s/"//g; $tissue =~ s/"//g;
				$quality = sprintf("%.5g", $quality);

				if ($genename =~ /NULL/) { $genename = "-"; }
				if ($dbsnp =~ /NULL/) { $dbsnp = " "; }
			
				print OUT "$library\t$chrom\t$position\t$ref\t$alt\t$quality\t$consequence\t$genename\t$geneid\t$feature\t$transcript\t$genetype\t$pposition\t$aachange\t$codon\t$dbsnp\t$class\t$zygosity\t$tissue\n";
			}
		
			open(STATUS, ">>$statusfile");
			print STATUS "End: $sublibraries\t$i\tDone: ",scalar (localtime),"\n";
			close (OUT); close (STATUS);
		}
	}
}
`rm -rf newtemp.txt`; #remove newtemp file

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
print "\n\n*********DONE*********\n\n";
# - - - - - - - - - - - - - - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -

#CONVERT TO VCF FILE
undef %HashDirectory;

foreach my $sublibraries (keys %LineHash){ # print $sublibraries,"\n";
  #print "dddd ",$sublibraries,"\n";
	(%GT, %TISSUE, %REF, %ALT, %QUAL, %CSQ, %DBSNP) = '';
  (%NEWQUAL, %NEWCSQ, %NEWREF, %NEWDBSNP, %NEWALT,%NEWGT) = '';
  (%subref, %subalt, %subgt) = '';
  ($headerinfo, %MTD) = '';
  opendir (DIR, $basepath) or die "Folder \"$basepath\" doesn't exist\n"; 
  my @Directory = readdir(DIR); close(DIR);
  #pushing each subfolder
  foreach (@Directory) {
    if ($_ =~ /^(\w.*)\.vcf$/) { #print "ccc ",$_;
      $HashDirectory{$1}= $1;
    }
  }
  unless (exists $HashDirectory{$LineHash{$sublibraries}}){
		print "yes ",$sublibraries,"\n\n";
		open(STATUS, ">>$statusfile") or die "Can not open statusfile\n";
		print STATUS "VCF: $sublibraries\tProcessing: ",scalar (localtime),"\n";
		close (STATUS);
		my $out = $basepath."/".$LineHash{$sublibraries}.".txt";
		PROCESS($out);
		SORTER();
		MTD();
		$headerinfo = HEADER();

		#GET OUTPUT FILE NAME
#		my $out = $sublibraries.".txt";
#fileparse($ARGV[0], qr/\.[^.]*(\.txt)?$/);
		my $output = $basepath."/".$LineHash{$sublibraries}.".vcf";

		#PRINT VCF OUTPUT
		open (OUT, ">", $output) or die "Can't open vcf file --weird\n";
		print OUT $headerinfo;
		foreach my $chrom (natsort keys %NEWREF) {
		  foreach my $position (sort {$a<=> $b} keys %{$NEWREF{$chrom}}) {
		    foreach my $ref (sort {$a<=> $b} keys %{$NEWREF{$chrom}{$position}}) {
		      print OUT "chr",$chrom,"\t",$position,"\t",$NEWDBSNP{$chrom}{$position}{$ref},"\t",$NEWREF{$chrom}{$position}{$ref},"\t";
		      print OUT $NEWALT{$chrom}{$position}{$ref},"\t",$NEWQUAL{$chrom}{$position}{$ref},"\tPASS\tCSQ=",$NEWCSQ{$chrom}{$position}{$ref};
		      print OUT ";MTD=",$MTD{$chrom}{$position}{$ref},"\tGT\t",$NEWGT{$chrom}{$position}{$ref};
		      print OUT "\n";
		    }
		  }
		}
		open(STATUS, ">>$statusfile");
		print STATUS "VCF: $sublibraries\tDone: ",scalar (localtime),"\n";
		close (STATUS); close (OUT);
	
	
		# STAGE 3: CONVERT TO VCF.GZ & TBI
		`bgzip $output`;
		`tabix -p vcf $output.gz`;
	}
}
	#STAGE 4: MOVE ALL THE FILES TO UCSC LOCATION
	#`rm -rf /home/modupe/public_html/UCSC/`;
	`mv $basepath/*vcf.gz* /home/modupe/public_html/UCSC`;
	`rm -rf $basepath`;
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
print "\n\n*********DONE*********\n\n";
# - - - - - - - - - - - - - - - - -S U B R O U T I N E S - - - - - - - - - - - - - -

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
      $CSQ{$line[1]}{$line[2]}{$line[3]}{$joint} =  $joint;
			
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
        foreach my $link (sort {$a <=> $b} keys %{$CSQ{$chrom}{$position}{$ref}}) {
          if (exists $NEWCSQ{$chrom}{$position}{$ref}){
            $NEWCSQ{$chrom}{$position}{$ref} .= ",$link";
					} else {
            $NEWCSQ{$chrom}{$position}{$ref} = $link;
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
					elsif ( $subgt{$chrom}{$position}{$ref}{'0/1'} < $subgt{$chrom}{$position}{$ref}{'1/2'} ) {
						$subgt{$chrom}{$position}{$ref}{'1/2'} =  $subgt{$chrom}{$position}{$ref}{'0/1'} + $subgt{$chrom}{$position}{$ref}{'1/2'};
					}
					elsif ( $subgt{$chrom}{$position}{$ref}{'0/1'} == $subgt{$chrom}{$position}{$ref}{'1/2'} ) {
						$subgt{$chrom}{$position}{$ref}{'0/1'} =  $subgt{$chrom}{$position}{$ref}{'0/1'} + $subgt{$chrom}{$position}{$ref}{'1/2'};
					}
					else{die "something is wrong";}
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

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
exit;

