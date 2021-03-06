#!usr/bin/perl
use strict;
use Pod::Usage;
use Getopt::Long;
use Sort::Key::Natural qw(natsort);

my ($split, $count, %REALGENOME,$header);

my $usage ="
=> Reference file (.fa) manipulation
Argument needed for '$0'
\t[ -split ] (split reference genome to individual chromosomes)
\t[ -count ] (measure individual chromosomal length)
<Reference genome file>
";

GetOptions('split|s'=>\$split, 'count|c'=>\$count  ) or die $usage;

die $usage unless ($split || $count);
@ARGV == 1 || die $usage;

if ($split){ SPLITGENOME($ARGV[0]); }
if ($count) { GENOMECOUNT($ARGV[0]); }

sub GENOME {
	$/ = ">";
        open (REF, "<", $_[0]) or die "Can't open $_[0]\n";
        my @genomefile = <REF>; close REF;
        shift @genomefile;
        foreach (@genomefile){
                my @pieces = split /\n/;
                my $total;
                foreach my $num (1..$#pieces){
                        $total .= $pieces[$num];
                }
                $pieces[0] =~ /(\S*).*/; $header = $1;
		unless ($header =~ /^chr/) {$header = "chr".$header; } 
		$total = substr($total,0,-1);
                $REALGENOME{$header} = $total;
        }
	$/ = "\n";
}
sub SPLITGENOME{
	open(OUT2, ">", "anotherall.fa");
	GENOME($_[0]);
	foreach my $key (natsort keys %REALGENOME){
		`mkdir -p other`;
		my $name = "other/$key".'.fa';
		open (OUT, ">", $name ) or die "Can't open $name\n";
		print OUT ">chr".$key."\n";
		print OUT "$REALGENOME{$key}\n";
		print OUT2 ">chr".$key."\n";
                print OUT2 "$REALGENOME{$key}\n";
		close OUT;
	}
}

sub GENOMECOUNT{
        GENOME($_[0]);
	my $name = "genomecount.txt";
        open (OUT, ">", $name ) or die "Can't open $name\n";
        foreach my $key (natsort keys %REALGENOME){
                my $name = "genomecount.txt"; 
                print OUT $key,"\t",length($REALGENOME{$key}),"\n";
        }
	close OUT;
}
