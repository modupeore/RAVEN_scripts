#!/usr/bin/perl

open(FILE,$ARGV[0]) or die "Cant open $ARGV[0]\n";
while(<FILE>){
 chomp;
 unless ($_ =~ /^ID/){
  @newlist = split("\t", $_);
  @ggalist = split("gga", $newlist[3]);
  foreach (1..$#ggalist){
   print $newlist[0]."\tgga".$ggalist[$_]."\n";
#if (exists $FILER{$newlist[0]}){
 #   $FILER{$newlist[0]} = $FILER{$newlist[0]}.",gga".$ggalist[$_];
 #  } else {
  #  $FILER{$newlist[0]} = "gga".$ggalist[$_];
   #}
   }
  }
}
#foreach (sort keys %FILER) {
 #print "$_\t$FILER{$_}]n";
