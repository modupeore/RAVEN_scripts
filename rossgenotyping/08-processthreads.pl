#!/usr/bin/perl
use strict;
use threads;
use Thread::Queue;
#convert the 07-htrzygty_window to multiple threads

my $files = `for i in \$(ls -1 ../RE*/hetero*txt);
do
c=\$(echo \$i | awk -F'-' '{print \$2}');
d=\$(echo \$i | awk -F'/' '{print \$2}');
echo perl ~/SCRIPTS/rossgenotyping/07-htrzygty_window.pl -t \$i -f ~/.big_ten/Galgal5/Galgal-5.fa -w 10000 -s 2000 -o \$d\\_\$c-10_20.txt;
echo perl ~/SCRIPTS/rossgenotyping/07-htrzygty_window.pl -t \$i -f ~/.big_ten/Galgal5/Galgal-5.fa -w 10000 -s 5000 -o \$d\\_\$c-10_50.txt;
echo perl ~/SCRIPTS/rossgenotyping/07-htrzygty_window.pl -t \$i -f ~/.big_ten/Galgal5/Galgal-5.fa -w 20000 -s 2000 -o \$d\\_\$c-20_20.txt;
echo perl ~/SCRIPTS/rossgenotyping/07-htrzygty_window.pl -t \$i -f ~/.big_ten/Galgal5/Galgal-5.fa -w 20000 -s 5000 -o \$d\\_\$c-20_50.txt;
echo perl ~/SCRIPTS/rossgenotyping/07-htrzygty_window.pl -t \$i -f ~/.big_ten/Galgal5/Galgal-5.fa -w 30000 -s 2000 -o \$d\\_\$c-30_20.txt;
echo perl ~/SCRIPTS/rossgenotyping/07-htrzygty_window.pl -t \$i -f ~/.big_ten/Galgal5/Galgal-5.fa -w 30000 -s 5000 -o \$d\\_\$c-30_50.txt;
echo perl ~/SCRIPTS/rossgenotyping/07-htrzygty_window.pl -t \$i -f ~/.big_ten/Galgal5/Galgal-5.fa -w 40000 -s 2000 -o \$d\\_\$c-40_20.txt;
echo perl ~/SCRIPTS/rossgenotyping/07-htrzygty_window.pl -t \$i -f ~/.big_ten/Galgal5/Galgal-5.fa -w 40000 -s 5000 -o \$d\\_\$c-40_50.txt;
echo perl ~/SCRIPTS/rossgenotyping/07-htrzygty_window.pl -t \$i -f ~/.big_ten/Galgal5/Galgal-5.fa -w 40000 -s 20000 -o \$d\\_\$c-40_20000.txt;
done`;

my @VAR = split(/\n/, $files);
my @threads;
my $queue = new Thread::Queue();
my $builder=threads->create(\&main); #create thread for each subarray into a thread
push @threads, threads->create(\&processor) for 1..5; #execute 5 threads
$builder->join; #join threads
foreach (@threads){$_->join;}
				
				
sub main {
  foreach my $count (0..$#VAR) {
		while(1) {
			if ($queue->pending() < 100) {
				$queue->enqueue($VAR[$count]);
				last;
			}
		}
	}
	foreach(1..5) { $queue-> enqueue(undef); }
}

sub processor {
	my $query;
	while ($query = $queue->dequeue()){
		print $query,"\n";
		system($query);
	}
}


