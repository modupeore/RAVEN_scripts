#!/usr/bin/perl
use Sort::Key::Natural qw(natsort);
#Order chromosome count file

open (IN, $ARGV[0]) or die "add chromosome distribution file";
while (<IN>) {
	chomp;
	@all = split /\t/;
	unless ($all[0] =~ /random/) {$GENOME{$all[0]} = $all[1];}
}
foreach $key (natsort keys %GENOME) { print $key,"\t",$GENOME{$key},"\n"; }
