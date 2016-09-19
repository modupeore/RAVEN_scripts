#!/usr/bin/perl
#table of the different output
$newdir = $ARGV[0]; %NEWKEY = ();
opendir(DIR, $newdir);
my @Directory = readdir(DIR); close(DIR);
foreach my $folder (@Directory){
  if ($folder =~ /[A-Z]/){
    print "$folder\n";
    $newfolder = "$newdir$folder/";
    opendir(FDIR, "$newfolder");
    my @FDirectory = readdir(FDIR); close(FDIR); 
    foreach my $FILE (@FDirectory) { 
      if ($FILE =~ /(\d+).*.vcf$/){
	#total
        $question = "grep \"^chr\" ".$folder."/".$1."-".$folder.".vcf | wc -l";
        $answer = `$question`; chomp $answer;
	$NEWKEY{$1}{"total"} = $answer;
        #not random
	$question = "grep \"^chr\" ".$newfolder.$1."-".$folder."_PASS.vcf | wc -l";
        $answer = `$question`; chomp $answer;
        $NEWKEY{$1}{"random"} = $answer;
        #count
	$question = "grep \"^chr\" ".$newfolder.$1."-".$folder."-raw_snp.vcf | wc -l";
        $answer = `$question`; chomp $answer;
        $NEWKEY{$1}{"snps"} = $answer;
        $question = "grep \"^chr\" ".$newfolder.$1."-".$folder."-raw_indel.vcf | wc -l";
        $answer = `$question`; chomp $answer;
        $NEWKEY{$1}{"indels"} = $answer;
        #count fil
        $question = "grep \"PASS\" ".$newfolder.$1."-".$folder."-filtered_snp.vcf | wc -l";
        $answer = `$question`; chomp $answer;
        $NEWKEY{$1}{"snpsfil"} = $answer;
        $question = "grep \"PASS\" ".$newfolder.$1."-".$folder."-filtered_indel.vcf | wc -l";
        $answer = `$question`; chomp $answer;
        $NEWKEY{$1}{"indelsfil"} = $answer;
        #count fil custom
        $question = "grep \"PASS\" ".$newfolder.$1."-".$folder."-filtered_snp30.vcf | wc -l";
        $answer = `$question`; chomp $answer;
        $NEWKEY{$1}{"snpsfil30"} = $answer;
        $question = "grep \"PASS\" ".$newfolder.$1."-".$folder."-filtered_indel60.vcf | wc -l";
        $answer = `$question`; chomp $answer;
        $NEWKEY{$1}{"indelsfil60"} = $answer;
      }
    }
    foreach $a (keys %NEWKEY){
		print $a."\t".$NEWKEY{$a}{"total"}."\t".$NEWKEY{$a}{"random"}."\t".$NEWKEY{$a}{"snps"}."\t".$NEWKEY{$a}{"snpsfil"}."\t".
				$NEWKEY{$a}{"snpsfil30"}."\t".$NEWKEY{$a}{"indels"}."\t".$NEWKEY{$a}{"indelsfil"}."\t".$NEWKEY{$a}{"indelsfil60"}."\n";
		delete $NEWKEY{$a};
	}
  }
}
