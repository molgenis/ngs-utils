#!/bin/sh
#PBS -N bmiedema_gonl_vcf_fill_err_log
#PBS -q devel
#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
#PBS -l mem=2gb
#PBS -l file=1gb
#PBS -e /gcc/groups/gcc/tmp01/bmiedema/data/gonl/filtered/err_log/err_log.err
#PBS -o /gcc/groups/gcc/tmp01/bmiedema/data/gonl/filtered/err_log/err_log.out

export PERL5LIB="/gcc/tools/Perl/lib/perl5/:${PERL5LIB}"

module load vcftools

cd /gcc/groups/gcc/tmp01/bmiedema/data/gonl/filtered

vcftools --gzvcf /gcc/groups/gcc/tmp01/bmiedema/data/gonl/unfiltered/gonl.chr13.snps_indels.r5.vcf.gz --bed /gcc/groups/gcc/tmp01/bmiedema/scripts/vcf-filter-batch/canonical_exon_sequences_filtered_headered_sorted.bed --out filtered.gonl.chr13.snps_indels.r5.vcf.gz --recode --remove-indels --recode-INFO-all > /gcc/groups/gcc/tmp01/bmiedema/data/gonl/filtered/filtered.gonl.chr13.snps_indels.r5.vcf.gz.log
bgzip -c ./filtered.gonl.chr13.snps_indels.r5.vcf.gz.recode.vcf > ./filtered.gonl.chr13.snps_indels.r5.vcf.gz.recode.vcf.gz
tabix -s 1 -b 2 -e 3 ./filtered.gonl.chr13.snps_indels.r5.vcf.gz.recode.vcf.gz
