#!/usr/bin/perl
open(DIR,"<", $ARGV[0]);
open (BB, ">birdlibraries");
open (FR, ">frnkmetadata");
open(GENE, ">genefpkm");
open(ISO, ">isoformfpkm");
open(VAR, ">varisum");
open (VRR, ">varires");
open (VAA, ">variann");
open (TRR, ">transum");
open (OUT2, ">restofthem");
while(<DIR>){
if($_ =~ /genes_fpkm/) {
	print GENE $_;
}
elsif($_ =~ /isoforms_fpkm/) {
        print ISO $_;
}
elsif($_ =~ /bird_libraries/) {
        print BB $_;
}
elsif($_ =~ /frnak_metadata/) {
        print FR $_;
}
elsif($_ =~ /transcripts_summary/) {
        print TRR $_;
}
elsif($_ =~ /variants_annotation/) {
        print VAA $_;
}
elsif($_ =~ /variants_result/) {
        print VRR $_;
}
elsif($_ =~ /variants_summary/) {
        print VAR $_;
}
else {
print OUT2 $_;
}
}
