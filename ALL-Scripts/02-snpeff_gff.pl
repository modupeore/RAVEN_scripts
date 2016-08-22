#!usr/bin/perl


my $shell = "#!/bin/sh
#PBS -N SNPeff-gff
#PBS -S /bin/bash
#PBS -V
#PBS -l walltime=300:00:00,cput=300:00:00,nodes=1
#PBS -q long
#PBS -d ";

my $folder="/home/amodupe/CARL/PERMANENT/snpEFFgff";
my $WKDIR = $ARGV[0];

open(OUT,">z-snpeffgff.sh");
print OUT $shell;
print OUT "$folder\n";

print OUT "time java -jar /usr/local/snpEff/snpEff.jar eff -v -i vcf -o gatk -s $WKDIR/all_UnifiedGenotyper_summary-gff.html Galgal4_76 $WKDIR/all_Unified_DP5.vcf > $WKDIR/all_Unified_DP5_snpEff-gff.vcf\n";

close(OUT);

print `qsub z-snpeffgff.sh`;
system "rm -rf z-snpeffgff.sh";
