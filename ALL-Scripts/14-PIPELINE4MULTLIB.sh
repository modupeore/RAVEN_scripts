#!/bin/bash
#JUST SPECIFY THE FOLDER THE MERGED BAM OUTPUT IS IN THIS SCRIPT

#PBS -N pipeline4mul-true
#PBS -S /bin/bash
#PBS -V
#PBS -l walltime=300:00:00,cput=3000:00:00,nodes=1:ppn=8
#PBS -d /home/amodupe/CARL/TORQUEoutput

echo start JOB_ID `date`

#changeable VARIABLES
WKDIR=/home/amodupe/CARL/FINAL/

#-----------------------------------

#permanent VARIABLES
PICARDDIR=/usr/local/picard-tools-1.67
GATKDIR=/usr/local/GATK-3.1.1
REF=/home/amodupe/CARL/Galgal4_genome.fa

#------------

##Process: 
#markduplicates
java -jar $PICARDDIR/MarkDuplicates.jar INPUT=$WKDIR/newbwasorted.bam OUTPUT=$WKDIR/newbwasorted_mdup.bam M=$WKDIR/newbwassorted_mdup.metrics CREATE_INDEX=true

#create sequence dictionary
java -jar /usr/local/picard-tools-1.67/CreateSequenceDictionary.jar R=$REF O=/home/amodupe/CARL/Galgal4_genome.dict

#FAIDX
#/usr/local/samtools-0.1.18/samtools faidx $REF

##GATK
java -jar $GATKDIR/GenomeAnalysisTK.jar -T UnifiedGenotyper -nt 6 -R $REF -I $WKDIR/newbwasorted_mdup.bam -o $WKDIR/all_Unified.vcf

#what i used for fasta files.
#####java -jar /usr/local/GATK-3.1.1/GenomeAnalysisTK.jar -T UnifiedGenotyper --defaultBaseQualities 40 -nt 6 -R Galgal4_genome.fa -I newbwasorted_mdup.bam -o all.vcf

#perl to select DP > 5 & get header information
perl /home/amodupe/CARL/Scripts/01-filteringDP5.pl $WKDIR $WKDIR/all_Unified.vcf

#snpEFF-gff
echo "Running snpEFF-GFF"
echo ""
perl /home/amodupe/CARL/Scripts/02-snpeff_gff.pl $WKDIR
echo "done"
echo ""
echo ""


#snpEFF-gtf
echo "Running snpEFF-GTF"
echo ""
perl /home/amodupe/CARL/Scripts/03-snpeff_gtf.pl $WKDIR
echo "done"
echo ""
echo end JOB_ID `date`
