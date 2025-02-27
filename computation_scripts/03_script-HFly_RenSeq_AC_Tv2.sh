#!/bin/bash

###########################
## Author: Narinder Singh
###########################

export PATH=$PATH:/homes/usr/bin:/homes/usr/bin/bin

#$ -l h_rt=72:00:00
#$ -l mem=15G
#$ -N HFly_RenSeq_AC_Tv2_parts
#$ -pe single 6
#$ -cwd

currDir=`pwd`

## setting up paths
cd $currDir
dbPath=/bulk/genomes/CS_NRGene/pseudomolecules_v1.0/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta
seq_files=/bulk/data/hf_KU2147_RenSeq/clean_data/sickle_trimmed/samples_Tv2

## aligning reads to ref genome using bwa mem
ls ${seq_files} | grep -v "R2.fastq" | sed -e 's/_R1.fastq//g' | \
	parallel bwa mem -t 6 $dbPath ${seq_files}/'{}'_R1.fastq ${seq_files}/'{}'_R2.fastq '>' '{}'.sam

ls | grep '.sam' | sed -e 's/.sam//g' | parallel samtools view -f2 -hb -o '{}'.bam '{}'.sam
ls | grep '.bam' | sed -e 's/.bam//g' | parallel samtools sort -T '{}'_sorted -o '{}'_sorted.bam '{}'.bam
ls | grep '_sorted.bam' | sed -e 's/_sorted.bam//g' | parallel samtools rmdup -S '{}'_sorted.bam '{}'_rmdup.bam

mkdir variants

# calling SNPs using samtools
ls | grep "_rmdup.bam" > rmdup.bam.files
samtools mpileup -IBg -b rmdup.bam.files -f $dbPath -t AD,DP,SP,ADF,ADR,DV -o variants/HFly_RenSeq_AC_Tv2_parts.bcf

cd variants

# creating vcf file with potential SNP sites
bcftools call -c -v HFly_RenSeq_AC_Tv2_parts.bcf > HFly_RenSeq_AC_Tv2_parts.vcf

# extracting allele count data
perl /homes/scripts/vcf2AC_v2.pl HFly_RenSeq_AC_Tv2_parts.vcf > HFly_RenSeq_AC_Tv2_parts.txt

