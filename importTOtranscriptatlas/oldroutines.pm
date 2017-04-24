
sub DBVARIANTSdeprecated{	
  print "\n\tINSERTING VARIANTS INTO THE DATABASE\n\n";
  #disconnecting and connecting again to database just incase
  $dbh->disconnect(); 
  $dbh = mysql();

  $_[1] =~ /^library_(\d*)$/;
  my $libnumber = "$1";
  my $folder = undef;

  #VEP file
  my @splitinput = split('\/', $_[0]);
  foreach my $i (0..$#splitinput-1){$folder.="$splitinput[$i]/";$i++;}
  my $information = fileparse($_[0], qr/(\.vcf)?$/);
  
  #variant_metadata
  my $gatk_version = ( ( split('\/',$GATKDIR)) [-2] ); my $vep_version = ( ( split('\/',$VEP)) [-4] ); my $picard_version = ( ( split('\/',$PICARDDIR)) [-2] );

  #OUTPUT file
  my $output = "$folder$information".".v1table";
  my $output2 = "$folder$information".".v2table";
  open(OUTDBVAR,">$output");
  open(OUTDBVAR2,">$output2");
  
  print OUTDBVAR "#CHR\tPOS\tREF\tALT\tQUAL\tCLASS\t";
  print OUTDBVAR "ZYGOSITY\tdbSNP\n";
  
  print OUTDBVAR2 "#CHR\tPOS\tCONQ\t";
  print OUTDBVAR2 "GENE\tSYMBOL\tTRANSCRIPT\t";
  print OUTDBVAR2 "FEATURE\tTYPE OF GENE\tPROTEIN POSITION\t";
  print OUTDBVAR2 "AA CHANGE\tCODON CHANGE\n";
  my $date = `date +%Y-%m-%d`;

  #initializing the hash tables . . .
  %VCFhash = ();
  %VEPhash = ();
  %ExtraWork = ();
  #running through subroutines . . . 
  VEPVARIANT($_[0]);
  my ($itsnp,$itindel,$itvariants) = 0;
  #printing to output table & variant_results table
  #VARIANT_SUMMARY
  $sth = $dbh->prepare("insert into variants_summary ( library_id, ANN_version, Picard_version, GATK_version, date ) values (?,?,?,?,?)");
  $sth ->execute($libnumber, $vep_version, $picard_version, $gatk_version, $date);
  #VARIANT_RESULTS
  foreach my $abc (sort keys %VCFhash) {
    foreach my $def (sort {$a <=> $b} keys %{ $VCFhash{$abc} }) {
      my @vcf = split('\|', $VCFhash{$abc}{$def});
      my @ext = split('\|', $ExtraWork{$abc}{$def});
      if ($vcf[3] =~ /,/){
        my $first = split(",",$vcf[1]);
        if (length $vcf[0] == length $first){
          $itvariants++; $itsnp++;
        }
        else {
          $itvariants++; $itindel++;
        }
      }
      elsif (length $vcf[0] == length $vcf[1]){
        $itvariants++; $itsnp++;
      }
      else {
        $itvariants++; $itindel++;
      }
              
      print OUTDBVAR "$abc\t$def\t$vcf[0]\t$vcf[1]\t$vcf[2]\t$ext[0]\t";
      print OUTDBVAR "$vcf[3]\t$ext[1]\n";
              
      #to variant_result
      $sth = $dbh->prepare("insert into variants_result ( library_id, chrom, position, ref_allele, alt_allele, quality, variant_class,
                           zygosity, existing_variant ) values (?,?,?,?,?,?,?,?,?)");
      $sth ->execute($libnumber, $abc, $def, $vcf[0], $vcf[1], $vcf[2], $ext[0], $vcf[3], $ext[1]);
    }
  }	
  close (OUTDBVAR);	
  foreach my $abc (sort keys %VEPhash) {
    foreach my $def (sort {$a <=> $b} keys %{ $VEPhash{$abc} }) {
      foreach my $ghi (sort keys %{ $VEPhash{$abc}{$def} }) {
        foreach my $jkl (sort keys %{ $VEPhash{$abc}{$def}{$ghi} }) {
          foreach my $mno (sort keys %{ $VEPhash{$abc}{$def}{$ghi}{$jkl} }){
            my @vep = split('\|', $VEPhash{$abc}{$def}{$ghi}{$jkl}{$mno});
            if(length($vep[4]) < 1){$vep[4] = "-";}
            if (length($vep[1]) < 1) {$vep[1] = "-";}
            print OUTDBVAR2 "$abc\t$def\t$ghi\t";
            print OUTDBVAR2 "$vep[1]\t$vep[8]\t$vep[2]\t";
            print OUTDBVAR2 "$vep[0]\t$vep[7]\t$vep[4]\t";
            print OUTDBVAR2 "$vep[5]\t$vep[6]\n";
            
            #to variants_annotation
            $sth = $dbh->prepare("insert into variants_annotation ( library_id, chrom, position, consequence, gene_id, gene_name,
                                transcript, feature, gene_type,protein_position, aminoacid_change, codon_change ) values
                                (?,?,?,?,?,?,?,?,?,?,?,?)");
            $sth ->execute($libnumber, $abc, $def, $ghi, $vep[1], $vep[8], $vep[2], $vep[0], $vep[7], $vep[4], $vep[5], $vep[6]);
          }
        }
      }
    }
  }
  close (OUTDBVAR2);
  
  #VARIANT_SUMMARY
  $syntax = "update variants_summary set total_VARIANTS = $itvariants, total_SNPS = $itsnp, 
  total_INDELS = $itindel, status = \'done\' where library_id like \"$libnumber\"";
  $sth = $dbh->prepare($syntax);
  $sth ->execute();
}
sub VEPVARIANTdeprecated { #import vep variants annotation to database
  #working on VEP variants
  my %Geneinfo = ''; my %Transcriptinfo = ''; my %Conqinfo = ''; my %Varclass = '';
  my %Featureinfo = ''; my %Proinfo = ''; my %Prochangeinfo = ''; my %Codoninfo = '';
  my %dbSNPinfo = ''; my %locinfo = ''; my %Entrezinfo = ''; my %GENEtype = '';
  my %GENEname = ''; my $position; my %location = '';
  unless(open(FILE,$_[0])){print "File \'$_[0]\' doesn't exist\n";exit;}
  my $verd;
  my @file = <FILE>;
  chomp @file;
  close (FILE);
  foreach my $chr (@file){
    unless ($chr =~ /^#/){
      my @chrdetails = split('\t', $chr);
      
      #removing the random chromosomes (scaffolds) - because no vital information can be found for them.
      my @morechrsplit = split(';', $chrdetails[7]);
      if (((split(':', $chrdetails[9]))[0]) eq '0/1'){$verd = "heterozygous";}
      elsif (((split(':', $chrdetails[9]))[0]) eq '1/1'){$verd = "homozygous";}
      elsif (((split(':', $chrdetails[9]))[0]) eq '1/2'){$verd = "heterozygous alternate";}
      
      #Processing the VEP section
      my @finalchrsplit = split("\,",(((split('=',((split(';',$chrdetails[7]))[-1])))[1]))); 
      foreach my $FCR (0..$#finalchrsplit){
        my @vepdetails = split('\|', $finalchrsplit[$FCR]);
        if ($vepdetails[1] !~ /WITHIN_NON_CODING_GENE/){
      
          #VCFhash information
          $VCFhash{$chrdetails[0]}{$chrdetails[1]} = "$chrdetails[3]|$chrdetails[4]|$chrdetails[5]|$verd";
      
          if ($vepdetails[14] < 1) {
            $vepdetails[14] = "-";
          }
          $location{"$chrdetails[0]|$chrdetails[1]|$vepdetails[1]|$vepdetails[14]|$vepdetails[4]"} = "$chrdetails[0]|$chrdetails[1]|$vepdetails[1]|$vepdetails[14]|$vepdetails[4]"; #specifying location with consequence.
          
          #GENE - 1
          if (exists $Geneinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]}){
            unless ($Geneinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} eq $vepdetails[4]){
              print "something is wrong with the code >> GENE";
              print "\na $chrdetails[0]\tb $chrdetails[1]\tc $vepdetails[14]\td $vepdetails[4]\te $vepdetails[4]\tf $Geneinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]}\n";		
            }
          }
          else {$Geneinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} = $vepdetails[4];}
          
          #TRANSCRIPT - 2
          if (exists $Transcriptinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]}){
            unless ($Transcriptinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} eq $vepdetails[6]){
              my $temp1details = "$Transcriptinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]},$vepdetails[6]";
              $Transcriptinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} = $temp1details;
            }
          }
          else {$Transcriptinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} = $vepdetails[6];}
          
          #CONSEQUENCE - 3
          if (exists $Conqinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} ){
            unless ($Conqinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} eq $vepdetails[1]){
              print "something is wrong with the code >> consequence";
              print "\na $chrdetails[0]\tb $chrdetails[1]\tc $vepdetails[14]\td $vepdetails[4]\te $vepdetails[1]\tf $Conqinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]}\n";		
              #exit;
            }
          }
          else {$Conqinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} = $vepdetails[1];}
          
          #VARIANT CLASS - 4
          if (exists $Varclass{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]}){
            unless ($Varclass{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} eq $vepdetails[20]){
              print "something is wrong with the code >> variantclass";
              print "\na $chrdetails[0]\tb $chrdetails[1]\tc $vepdetails[14]\td $vepdetails[4]\te $vepdetails[20]\tf $Varclass{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]}\n";		
              #exit;
            }
          }
          else {$Varclass{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} = $vepdetails[20];}
         
          #FEATURE TYPE - 5
          if (exists $Featureinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]}){
            unless ($Featureinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} eq $vepdetails[5]){
              print "something is wrong with the code >> feature type";
              print "\na $chrdetails[0]\tb $chrdetails[1]\tc $vepdetails[14]\td $vepdetails[4]\te $vepdetails[5]\tf $Featureinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]}\n";		
            }
          }
          else {$Featureinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} = $vepdetails[5];}
          
          #PROTEIN POSITION - 6
          if (exists $Proinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]}){
            unless ($Proinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} eq $vepdetails[14]){
              print "something is wrong with the code >> protein position";
              print "\na $chrdetails[0]\tb $chrdetails[1]\tc $vepdetails[14]\td $vepdetails[4]\te $vepdetails[14]\tf $Proinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]}\n";		
            }		
          }
          else {$Proinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} = $vepdetails[14];}
         
          #PROTEIN CHANGE - 7 (or Amino Acid)
          if (exists $Prochangeinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]}){
            unless ($Prochangeinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} eq $vepdetails[15]){
              my $temp2details = "$Prochangeinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]},$vepdetails[15]";
              $Prochangeinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} = $temp2details;
            }
          }
          else {$Prochangeinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} = $vepdetails[15];}
          
          #CODON CHANGE - 8
          if (exists $Codoninfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]}){
            unless ($Codoninfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} eq $vepdetails[16]){
              my $temp3details = "$Codoninfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]},$vepdetails[16]";
              $Codoninfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} = $temp3details;
            }		
          }
          else {$Codoninfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} = $vepdetails[16];}
          
          #dbSNP - 9
          if (exists $dbSNPinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]}){
            unless ($dbSNPinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} eq $vepdetails[17]){
              print "something is wrong with the code >> dbSNP";
              print "\na $chrdetails[0]\tb $chrdetails[1]\tc $vepdetails[14]\td $vepdetails[4]\te $vepdetails[17]\tf $dbSNPinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]}\n";		
            }		
          }
          else {$dbSNPinfo{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} = $vepdetails[17];}
          
          #GENE name - 10
          if (exists $GENEname{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]}){
            unless ($GENEname{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} eq $vepdetails[3]){
              print "something is wrong with the code >> genename";
              print "\na $chrdetails[0]\tb $chrdetails[1]\tc $vepdetails[14]\td $vepdetails[4]\te $vepdetails[3]\tf $GENEname{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]}\n";		
            }	
          }
          else {$GENEname{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} = $vepdetails[3];}
          
          #GENE type - 11
          if (exists $GENEtype{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]}){
            unless ($GENEtype{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} eq $vepdetails[7]){
              my $temp4details = "$GENEtype{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]},$vepdetails[7]";
              $GENEtype{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} = $temp4details;
            }	
          }
          else {$GENEtype{$chrdetails[0]}{$chrdetails[1]}{$vepdetails[1]}{$vepdetails[14]}{$vepdetails[4]} = $vepdetails[7];}
        }
      }
    }
  }

  foreach my $alldetails (keys %location){
        my ($chrdetails1, $chrdetails2, $consequences, $pposition, $genename) = split('\|', $alldetails);
        #cleaning up the text
        my $clean1 = CLEANUP($Varclass{$chrdetails1}{$chrdetails2}{$consequences}{$pposition}{$genename});
        my $clean2 = CLEANUP($Featureinfo{$chrdetails1}{$chrdetails2}{$consequences}{$pposition}{$genename});
        my $clean3 = CLEANUP($Geneinfo{$chrdetails1}{$chrdetails2}{$consequences}{$pposition}{$genename});
        my $clean4 = CLEANUP($Transcriptinfo{$chrdetails1}{$chrdetails2}{$consequences}{$pposition}{$genename});
        my $clean5 = CLEANUP($Conqinfo{$chrdetails1}{$chrdetails2}{$consequences}{$pposition}{$genename});
        my $clean6 = CLEANUP($Proinfo{$chrdetails1}{$chrdetails2}{$consequences}{$pposition}{$genename});
        my $clean7 = CLEANUP($Prochangeinfo{$chrdetails1}{$chrdetails2}{$consequences}{$pposition}{$genename});
        my $clean8 = CLEANUP($Codoninfo{$chrdetails1}{$chrdetails2}{$consequences}{$pposition}{$genename});
        my $clean9 = CLEANUP($dbSNPinfo{$chrdetails1}{$chrdetails2}{$consequences}{$pposition}{$genename});
        my $clean10 = CLEANUP($GENEtype{$chrdetails1}{$chrdetails2}{$consequences}{$pposition}{$genename});
        my $clean11 = CLEANUP($GENEname{$chrdetails1}{$chrdetails2}{$consequences}{$pposition}{$genename});
        $VEPhash{$chrdetails1}{$chrdetails2}{$consequences}{$pposition}{$genename} = "$clean2|$clean3|$clean4|$clean5|$clean6|$clean7|$clean8|$clean10|$clean11";
        $ExtraWork{$chrdetails1}{$chrdetails2} = "$clean1|$clean9";
  }   
}
sub DBSNPDATVARIANTS{	
  print "\n\tINSERTING SNPDAT VARIANTS INTO THE DATABASE\n\n";
  #disconnecting and connecting again to database just incase
  $dbh->disconnect(); 
  $dbh = mysql();

  $_[2] =~ /^library_(\d*)$/;
  my $libnumber = $1;
  my $folder = undef;

  #VCF file
  my @splitinput = split('\/', $_[1]);
  foreach my $i (0..$#splitinput-1){$folder.="$splitinput[$i]/";$i++;}
  my $information = fileparse($_[0], qr/(\.vcf)?$/);
  
  #variant_metadata
   my $gatk_version = ( ( split('\/',$GATKDIR)) [-2] ); my $ann_version = ( ( split('\/',$SNPdat)) [-2] ); my $picard_version = ( ( split('\/',$PICARDDIR)) [-2] );

  #OUTPUT file
  my $output = "$folder$information".".v1table";
  my $output2 = "$folder$information".".v2table";
  
  open(OUTDBVAR,">$output");
  open(OUTDBVAR2,">$output2");
  
  print OUTDBVAR "#CHR\tPOS\tREF\tALT\tQUAL\tCLASS\t";
  print OUTDBVAR "ZYGOSITY\tdbSNP\n";
  
  print OUTDBVAR2 "#CHR\tPOS\tCONQ\t";
  print OUTDBVAR2 "GENE\tSYMBOL\tTRANSCRIPT\t";
  print OUTDBVAR2 "FEATURE\tTYPE OF GENE\tPROTEIN POSITION\t";
  print OUTDBVAR2 "AA CHANGE\tCODON CHANGE\n";
  my $date = `date +%Y-%m-%d`;

  #initializing the hash tables . . .
  %VCFhash = ();
  %VEPhash = ();
  %ExtraWork = ();
  #running through subroutines . . . 
  SNPdatVARIANT($_[0], $_[1]);
  my ($itsnp,$itindel,$itvariants) = 0;
  #printing to output table & variant_results table
  #VARIANT_SUMMARY
  $sth = $dbh->prepare("insert into variants_summary ( library_id, ANN_version, Picard_version, GATK_version, date ) values (?,?,?,?,?)");
  $sth ->execute($libnumber, $ann_version, $picard_version, $gatk_version, $date);

  #VARIANT_RESULTS
  foreach my $abc (sort keys %VCFhash) {
    foreach my $def (sort {$a <=> $b} keys %{ $VCFhash{$abc} }) {
      my @vcf = split('\|', $VCFhash{$abc}{$def});
      my $variantclass;
      my $first = split(",",$vcf[1]);
      if ($vcf[3] =~ /,/){                   
        if (length($vcf[0]) == length($first)){
          $itvariants++; $itsnp++; $variantclass = "SNV";
        }
        elsif (length($vcf[0]) > length($first)) {
          $itvariants++; $itindel++; $variantclass = "deletion";
        }
        elsif (length($vcf[0]) < length($first)) {
          $itvariants++; $itindel++; $variantclass = "insertion";
        }
        else {
          print "first False variant, cross-check asap"; exit; 
        }
      }
      elsif (length $vcf[0] == length $vcf[1]){
        $itvariants++; $itsnp++; $variantclass = "SNV";
      }
      elsif (length $vcf[0] > length ($vcf[1])) {
        $itvariants++; $itindel++; $variantclass = "deletion";
      }
      elsif (length $vcf[0] < length ($vcf[1])) {
        $itvariants++; $itindel++; $variantclass = "insertion";
      }
      else {
        print "second False variant, cross-check asap\n"; exit;
      }
      print OUTDBVAR "$abc\t$def\t$vcf[0]\t$vcf[1]\t$vcf[2]\t$variantclass\t";
      print OUTDBVAR "$vcf[3]\t$ExtraWork{$abc}{$def}\n";
              
      $sth = $dbh->prepare("insert into variants_result ( library_id, chrom, position, ref_allele, alt_allele, quality, variant_class,
                          zygosity, existing_variant ) values (?,?,?,?,?,?,?,?,?)");
      $sth ->execute($libnumber, $abc, $def, $vcf[0], $vcf[1], $vcf[2], $variantclass, $vcf[3], $ExtraWork{$abc}{$def});                  
    }
  }
  close (OUTDBVAR);
  
  foreach my $abc (sort keys %VEPhash) {
    foreach my $def (sort {$a <=> $b} keys %{ $VEPhash{$abc} }) {
      foreach my $ghi (sort keys %{ $VEPhash{$abc}{$def} }) {
          foreach my $jkl (sort keys %{ $VEPhash{$abc}{$def}{$ghi} }) {
            my @vep = split('\|', $VEPhash{$abc}{$def}{$ghi}{$jkl});
            
            if(length($vep[4]) < 1){
                $vep[4] = "-";
            }
            if (length($vep[1]) < 1) {
                $vep[1] = "-";
            }
              print OUTDBVAR2 "$abc\t$def\t$ghi\t";
              print OUTDBVAR2 "$vep[1]\t$vep[8]\t$vep[2]\t";
              print OUTDBVAR2 "$vep[0]\t$vep[7]\t$vep[4]\t";
              print OUTDBVAR2 "$vep[5]\t$vep[6]\n";
              
              #to variants_annotation
              $sth = $dbh->prepare("insert into variants_annotation ( library_id, chrom, position, consequence, gene_id, gene_name,
                                 transcript, feature, gene_type,protein_position, aminoacid_change, codon_change ) values
                                (?,?,?,?,?,?,?,?,?,?,?,?)");
              $sth ->execute($libnumber, $abc, $def, $ghi, $vep[1], $vep[8], $vep[2], $vep[0], $vep[7], $vep[4], $vep[5], $vep[6]);
          }
      }
    }
  }
  close (OUTDBVAR2);
  
  #VARIANT_SUMMARY
  $syntax = "update variants_summary set total_VARIANTS = $itvariants, total_SNPS = $itsnp, 
  total_INDELS = $itindel, status = \'done\' where library_id like \"$libnumber\"";
  $sth = $dbh->prepare($syntax);
  $sth ->execute();
}
sub SNPdatVARIANT {
  #working on VEP variants
  my %Geneinfo = ''; my %Transcriptinfo = ''; my %Conqinfo = '';
  my %Featureinfo = ''; my %Proinfo = ''; my %Prochangeinfo = ''; my %Codoninfo = '';
  my %dbSNPinfo = ''; my %locinfo = ''; my %Entrezinfo = ''; my %GENEtype = '';
  my %GENEname = ''; my $position; my %location = ''; my ($consqofgene,$idofgene, $codonchange, $aminoacid);
	
  unless(open(VAR,$_[0])){print "File \'$_[0]\' doesn't exist\n";exit;}
  my $verd; my @variantfile = <VAR>; chomp @variantfile; close (VAR);
  foreach my $chr (@variantfile){
    unless ($chr =~ /^#/){
      my @chrdetail = split('\t', $chr);
      #removing the random chromosomes (scaffolds) - because no vital information can be found for them.
      my @morechrsplit = split(';', $chrdetail[7]);
      if (((split(':', $chrdetail[9]))[0]) eq '0/1'){$verd = "heterozygous";}
      elsif (((split(':', $chrdetail[9]))[0]) eq '1/1'){$verd = "homozygous";}
      elsif (((split(':', $chrdetail[9]))[0]) eq '1/2'){$verd = "heterozygous almyternate";}
      #VCFhash information
      $VCFhash{$chrdetail[0]}{$chrdetail[1]} = "$chrdetail[3]|$chrdetail[4]|$chrdetail[5]|$verd";
    }
  }
  unless(open(FILE,$_[1])){print "File \'$_[0]\' doesn't exist\n";exit;}
  my @file = <FILE>; chomp @file; close (FILE); shift(@file);
  foreach my $chr2 (@file){
    #Processing the ANNotation section
    my @chrdetails = split('\t', $chr2);
    $chrdetails[0] = "\L$chrdetails[0]";
    if ($chrdetails[10] =~ /ID=AMISG/){
      #getting gene id and name of gene
      my @geneidlist = split("\;", $chrdetails[10]);
      foreach (@geneidlist){
        if($_ =~ /^ID=AMISG/){$idofgene = substr($_,3);}
        last;
      }
              
      #change codon format to be consistent
      unless ($chrdetails[19] eq "NA"){
        $chrdetails[19] = "\L$chrdetails[19]";
    
        my $aq = ((split('\[', $chrdetails[19],2))[0]); my $bq = ((split('\]', $chrdetails[19],2))[1]); 
        my $dq = uc ((split('\[', ((split('\/', $chrdetails[19],2))[0]),2)) [1]); 
        my $eq = uc ((split('\]',((split('\/', $chrdetails[19],2))[1]),2)) [0]);
        $codonchange = "$aq$dq$bq/$aq$eq$bq";
        $aminoacid = substr(substr($chrdetails[20],1),0,-1);
      }else {
        $codonchange = undef; $aminoacid = undef;
      }
      #determine the consequence of the SNP change
      if($chrdetails[21] eq "N") {
        if ($aminoacid =~ /^M/){$consqofgene = "START_LOST";}
        elsif($aminoacid =~ /M$/){$consqofgene = "START_GAINED";}
        elsif($aminoacid =~/^\-/){$consqofgene = "STOP_LOST";}
        elsif($aminoacid =~ /\-$/){$consqofgene = "STOP_GAINED";}	
        else {$consqofgene = "NON_SYNONYMOUS_CODING";}
      }
      elsif ($chrdetails[21] eq "Y"){
        $consqofgene = "SYNONYMOUS_CODING";
        $aminoacid = ((split('\/', $aminoacid,2))[0]);
      }
      elsif($chrdetails[21] eq "NA") {
        if ($chrdetails[3] =~ /Intergenic/){$consqofgene = "INTERGENIC";}
        elsif($chrdetails[3] =~ /Intronic/){$consqofgene = "INTRONIC";}
        elsif($chrdetails[3] =~ /Exonic/){$consqofgene = "EXON";}
        else {$consqofgene = "NA";}
      }
      else {
        print "There's a problem with the consequence of the gene\n\n"; exit;
      }
      foreach (0..$#chrdetails){if ($chrdetails[$_] eq "NA"){ $chrdetails[$_] = undef;}} #converting NA to undef
              
      $location{"$chrdetails[0]|$chrdetails[1]|$consqofgene|$idofgene"} = "$chrdetails[0]|$chrdetails[1]|$consqofgene|$idofgene"; #specifying location with consequence.
              
      #GENE - 1
      if (exists $Geneinfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene}){
        unless ($Geneinfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} eq $idofgene){
          print "something is wrong with the code >> Geneinfo\n";
          exit;
        }
      }
      else {	
        $Geneinfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} = $idofgene;
      }
      
      #TRANSCRIPT - 2
      if (exists $Transcriptinfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene}){
        unless ($Transcriptinfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} eq $AMIST{$idofgene}){
          print "something is wrong with the code >> transcript\n";
          exit;
        }
      }
      else {	
        $Transcriptinfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} = $AMIST{$idofgene};
      } 	
      
      #CONSEQUENCE - 3
      if (exists $Conqinfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} ){
        unless ($Conqinfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} eq $consqofgene){
          print "something is wrong with the code >> consequence\n";
          exit;
        }			
      }
      else {	
        $Conqinfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} = $consqofgene;
      }	
        
      #FEATURE TYPE - 5 #SNPdat
      if (exists $Featureinfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene}){
        unless ($Featureinfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} eq $chrdetails[5]){
          print "something is wrong with the code >> feature\n";
          exit;
        }			
      }
      else {	
        $Featureinfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} = $chrdetails[5];
      }

      #PROTEIN POSITION - 6 #SNPdat this isn't provided
      $Proinfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} = undef;

      #PROTEIN CHANGE - 7 (or Amino Acid) #SNPdat
      if (exists $Prochangeinfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene}){
        unless ($Prochangeinfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} eq $aminoacid){
          my $temp1details = "$Prochangeinfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene},$aminoacid";
          $Prochangeinfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} = $temp1details;
        }			
      }
      else {	
        $Prochangeinfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} = $aminoacid;
      }
	
      #CODON CHANGE - 8 #SNPdat
      if (exists $Codoninfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene}){
        unless ($Codoninfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} eq $codonchange){
          my $temp2details = "$Codoninfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene},$codonchange";
          $Codoninfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} = $temp2details;
        }			
      }
      else {	
        $Codoninfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} = $codonchange;
      }
      
	#dbSNP - 9 #SNPdat # no dbSNP info for alligator
	$dbSNPinfo{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} = undef;			
	
	#GENE name - 10 #name of gene from the alligator gene list. #SNPdat
      if (exists $GENEname{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene}){
        unless ($GENEname{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} eq $AMISG{$idofgene}){
          print "something is wrong with the code >> gene name\n";
          exit;
        }
	}
	else {	
	  $GENEname{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} = $AMISG{$idofgene};
	}
	
	#GENE type - 11 #SNPdat
      if (exists $GENEtype{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene}){
        unless ($GENEtype{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} eq $chrdetails[3]){
          print "something is wrong with the code >> gene type\n";
          exit;
        }
	}
	else {	
	  $GENEtype{$chrdetails[0]}{$chrdetails[1]}{$consqofgene}{$idofgene} = $chrdetails[3];
	}
    }
  }

  foreach my $alldetails (keys %location){
    my ($chrdetails1, $chrdetails2, $consequences, $genename) = split('\|', $alldetails);
    #cleaning up the text
    my $clean2 = CLEANUP($Featureinfo{$chrdetails1}{$chrdetails2}{$consequences}{$genename});
    my $clean3 = CLEANUP($Geneinfo{$chrdetails1}{$chrdetails2}{$consequences}{$genename});
    my $clean4 = CLEANUP($Transcriptinfo{$chrdetails1}{$chrdetails2}{$consequences}{$genename});
    my $clean5 = CLEANUP($Conqinfo{$chrdetails1}{$chrdetails2}{$consequences}{$genename});
    my $clean6 = CLEANUP($Proinfo{$chrdetails1}{$chrdetails2}{$consequences}{$genename});
    my $clean7 = CLEANUP($Prochangeinfo{$chrdetails1}{$chrdetails2}{$consequences}{$genename});
    my $clean8 = CLEANUP($Codoninfo{$chrdetails1}{$chrdetails2}{$consequences}{$genename});
    my $clean9 = CLEANUP($dbSNPinfo{$chrdetails1}{$chrdetails2}{$consequences}{$genename});
    my $clean10 = CLEANUP($GENEtype{$chrdetails1}{$chrdetails2}{$consequences}{$genename});
    my $clean11 = CLEANUP($GENEname{$chrdetails1}{$chrdetails2}{$consequences}{$genename});
    $VEPhash{$chrdetails1}{$chrdetails2}{$consequences}{$genename} = "$clean2|$clean3|$clean4|$clean5|$clean6|$clean7|$clean8|$clean10|$clean11";
    $ExtraWork{$chrdetails1}{$chrdetails2} = $clean9;
  } 
}


sub CLEANUP {
  #cleaning up the VEP variants so that it doesn't have repetitions in the output
  my @unclean = split(',', $_[0]);
  my ($cleansyntax, %Hashdetails, %Hash2) = undef;
  foreach (0..$#unclean){
    if ($unclean[$_] =~ /^[a-zA-Z0-9]/){
      $Hashdetails{$unclean[$_]} = $_;
    }
  }
  foreach my $unique (keys %Hashdetails){
    if ($unique =~ /^[a-zA-Z0-9]/){
      $Hash2{$Hashdetails{$unique}} = $unique;
    }
  }
  foreach my $final (sort keys %Hash2){
    if ($Hash2{$final} =~ /^[a-zA-Z0-9]/){
      $cleansyntax .= ",$Hash2{$final}";
    }
  }
  my $returnclean = substr($cleansyntax, 1);
  return $returnclean;
}
sub ALLIGATOR {
  my ($theparent, $theid); 
  unless(open(FILE,$_[0])){exit "File \'$_[0]\' doesn't exist\n";}
  my @file = <FILE>;
  chomp @file; shift (@file);
  close (FILE);
  foreach my $line (@file){
    my @details = split(',', $line);
    $AMISG{$details[0]} = $details[1];
  } 
  unless(open(FILER,$_[1])){exit "File \'$_[1]\' doesn't exist\n";}
  my @filer = <FILER>;
  chomp @filer; shift (@filer);
  close (FILER);
  foreach (@filer){
    my @newdetails = split('\t', $_);
    if ($newdetails[8] =~ /Parent=AMISG/){
	my @whatIwant = split("\;", $newdetails[8]);
	foreach (@whatIwant){
	  if ($_ =~ /ID=/){
	    $theid = substr($_, 3);
        }
	  elsif ($_ =~ /Parent=/){
	    $theparent = substr($_, 7);
        }
	}
	$AMIST{$theparent} = $theid;
    }
  }    
}
sub SUMMARYstmts{
  #transcripts_summary from the database
  $dbh = mysql();
  print "\n\tEXECUTING SELECT STATEMENT ON THE DATABASE TABLES \n";
  $syntax = "select count(*) from transcripts_summary";
  $sth = $dbh->prepare($syntax);
  $sth->execute or die "SQL Error: $DBI::errstr\n";
  while (@row = $sth->fetchrow_array() ) {print "Number of rows in \"transcripts_summary\" table \t:\t @row\n";}
  #genes_fpkm
  $syntax = "select count(*) from genes_fpkm";
  $sth = $dbh->prepare($syntax);
  $sth->execute or die "SQL Error: $DBI::errstr\n";
  while (@row = $sth->fetchrow_array() ) {print "Number of rows in \"genes_fpkm\" table \t\t:\t @row\n";}
  #isoforms_fpkm
  $syntax = "select count(*) from isoforms_fpkm";
  $sth = $dbh->prepare($syntax);
  $sth->execute or die "SQL Error: $DBI::errstr\n";
  while (@row = $sth->fetchrow_array() ) {print "Number of rows in \"isoforms_fpkm\" table \t:\t @row\n";}
  #variant_summary
  $syntax = "select count(*) from variants_summary";
  $sth = $dbh->prepare($syntax);
  $sth->execute or die "SQL Error: $DBI::errstr\n";
  while (@row = $sth->fetchrow_array() ) {print "Number of rows in \"variants_summary\" table \t:\t @row\n";}
  #variant_list
  $syntax = "select count(*) from variants_result";
  $sth = $dbh->prepare($syntax);
  $sth->execute or die "SQL Error: $DBI::errstr\n";
  while (@row = $sth->fetchrow_array() ) {print "Number of rows in \"variants_results\" table \t:\t @row\n";}
  #variant_annotion
  $syntax = "select count(*) from variants_annotation";
  $sth = $dbh->prepare($syntax);
  $sth->execute or die "SQL Error: $DBI::errstr\n";
  while (@row = $sth->fetchrow_array() ) {print "Number of rows in \"variants_annotation\" table \t:\t @row\n";}
  #frnak_metadata
  $syntax = "select count(*) from frnak_metadata";
  $sth = $dbh->prepare($syntax);
  $sth->execute or die "SQL Error: $DBI::errstr\n";
  while (@row = $sth->fetchrow_array() ) {print "Number of rows in \"frnak_metadata\" table \t:\t @row\n";}
}

