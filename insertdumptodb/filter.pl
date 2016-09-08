#!/usr/bin/perl

open(VAR,"<", $ARGV[0]);
open(OUT1, ">bbbbb.txt");
open(OUT2, ">ccccc.txt");
while (<VAR>) {
  chomp;
  @values = split('VALUES \(', $_);
  @indit = split(",",$values[1]);
  if ($indit[0] >= 1522){
	print OUT1 "$_\n";
  }
  else {print OUT2 "$_\n";}
}
