#!/bin/sh
#PBS -N bmiedema_1000g_vcf_fill_err_log
#PBS -q devel
#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
#PBS -l mem=2gb
#PBS -l file=1gb
#PBS -e /gcc/groups/gcc/tmp01/bmiedema/data/1000G/filtered/err_log/err_log.err
#PBS -o /gcc/groups/gcc/tmp01/bmiedema/data/1000G/filtered/err_log/err_log.out

export PERL5LIB="/gcc/tools/Perl/lib/perl5/:${PERL5LIB}"

module load vcftools

cd /gcc/groups/gcc/tmp01/bmiedema/data/1000G/filtered

vcftools --gzvcf /gcc/groups/gcc/tmp01/bmiedema/data/1000G/unfiltered/gtc.ALL.chrY.phase3_integrated.20130502.sites.vcf.gz.gz --bed /gcc/groups/gcc/tmp01/bmiedema/scripts/vcf-filter-batch/canonical_exon_sequences_filtered_headered_sorted.bed --out filtered.gtc.ALL.chrY.phase3_integrated.20130502.sites.vcf.gz.gz --recode --remove-indels --recode-INFO-all > /gcc/groups/gcc/tmp01/bmiedema/data/1000G/filtered/filtered.gtc.ALL.chrY.phase3_integrated.20130502.sites.vcf.gz.gz.log
bgzip -c ./filtered.gtc.ALL.chrY.phase3_integrated.20130502.sites.vcf.gz.gz.recode.vcf > ./filtered.gtc.ALL.chrY.phase3_integrated.20130502.sites.vcf.gz.gz.recode.vcf.gz
tabix -s 1 -b 2 -e 3 ./filtered.gtc.ALL.chrY.phase3_integrated.20130502.sites.vcf.gz.gz.recode.vcf.gz
