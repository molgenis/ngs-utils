#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use Cwd;

###To remove incorrect header from VCF files use: sed '2,25d' <filename> > <new file>
###Only the first line describing the VCF version and the line containing column headers should be kept

my $lines;
my @array;
my $element;
my $length;
my @snps;
my $count = 0;
my $linecount = 0;

open (OUTPUTchr1, ">/data/p255198/resources/gatk_resources/resources/1kg.pilot_release.merged.indels.sites.hg19.chr1.vcf");
open (OUTPUTchr2, ">/data/p255198/resources/gatk_resources/resources/1kg.pilot_release.merged.indels.sites.hg19.chr2.vcf");
open (OUTPUTchr3, ">/data/p255198/resources/gatk_resources/resources/1kg.pilot_release.merged.indels.sites.hg19.chr3.vcf");
open (OUTPUTchr4, ">/data/p255198/resources/gatk_resources/resources/1kg.pilot_release.merged.indels.sites.hg19.chr4.vcf");
open (OUTPUTchr5, ">/data/p255198/resources/gatk_resources/resources/1kg.pilot_release.merged.indels.sites.hg19.chr5.vcf");
open (OUTPUTchr6, ">/data/p255198/resources/gatk_resources/resources/1kg.pilot_release.merged.indels.sites.hg19.chr6.vcf");
open (OUTPUTchr7, ">/data/p255198/resources/gatk_resources/resources/1kg.pilot_release.merged.indels.sites.hg19.chr7.vcf");
open (OUTPUTchr8, ">/data/p255198/resources/gatk_resources/resources/1kg.pilot_release.merged.indels.sites.hg19.chr8.vcf");
open (OUTPUTchr9, ">/data/p255198/resources/gatk_resources/resources/1kg.pilot_release.merged.indels.sites.hg19.chr9.vcf");
open (OUTPUTchr10, ">/data/p255198/resources/gatk_resources/resources/1kg.pilot_release.merged.indels.sites.hg19.chr10.vcf");
open (OUTPUTchr11, ">/data/p255198/resources/gatk_resources/resources/1kg.pilot_release.merged.indels.sites.hg19.chr11.vcf");
open (OUTPUTchr12, ">/data/p255198/resources/gatk_resources/resources/1kg.pilot_release.merged.indels.sites.hg19.chr12.vcf");
open (OUTPUTchr13, ">/data/p255198/resources/gatk_resources/resources/1kg.pilot_release.merged.indels.sites.hg19.chr13.vcf");
open (OUTPUTchr14, ">/data/p255198/resources/gatk_resources/resources/1kg.pilot_release.merged.indels.sites.hg19.chr14.vcf");
open (OUTPUTchr15, ">/data/p255198/resources/gatk_resources/resources/1kg.pilot_release.merged.indels.sites.hg19.chr15.vcf");
open (OUTPUTchr16, ">/data/p255198/resources/gatk_resources/resources/1kg.pilot_release.merged.indels.sites.hg19.chr16.vcf");
open (OUTPUTchr17, ">/data/p255198/resources/gatk_resources/resources/1kg.pilot_release.merged.indels.sites.hg19.chr17.vcf");
open (OUTPUTchr18, ">/data/p255198/resources/gatk_resources/resources/1kg.pilot_release.merged.indels.sites.hg19.chr18.vcf");
open (OUTPUTchr19, ">/data/p255198/resources/gatk_resources/resources/1kg.pilot_release.merged.indels.sites.hg19.chr19.vcf");
open (OUTPUTchr20, ">/data/p255198/resources/gatk_resources/resources/1kg.pilot_release.merged.indels.sites.hg19.chr20.vcf");
open (OUTPUTchr21, ">/data/p255198/resources/gatk_resources/resources/1kg.pilot_release.merged.indels.sites.hg19.chr21.vcf");
open (OUTPUTchr22, ">/data/p255198/resources/gatk_resources/resources/1kg.pilot_release.merged.indels.sites.hg19.chr22.vcf");
open (OUTPUTchr23, ">/data/p255198/resources/gatk_resources/resources/1kg.pilot_release.merged.indels.sites.hg19.chr23.vcf");
open (OUTPUTchr24, ">/data/p255198/resources/gatk_resources/resources/1kg.pilot_release.merged.indels.sites.hg19.chr24.vcf");
open (OUTPUTchr25, ">/data/p255198/resources/gatk_resources/resources/1kg.pilot_release.merged.indels.sites.hg19.chr25.vcf");

print OUTPUTchr1 "##fileformat=VCFv4.0\n" . "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n";
print OUTPUTchr2 "##fileformat=VCFv4.0\n" . "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n";
print OUTPUTchr3 "##fileformat=VCFv4.0\n" . "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n";
print OUTPUTchr4 "##fileformat=VCFv4.0\n" . "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n";
print OUTPUTchr5 "##fileformat=VCFv4.0\n" . "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n";
print OUTPUTchr6 "##fileformat=VCFv4.0\n" . "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n";
print OUTPUTchr7 "##fileformat=VCFv4.0\n" . "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n";
print OUTPUTchr8 "##fileformat=VCFv4.0\n" . "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n";
print OUTPUTchr9 "##fileformat=VCFv4.0\n" . "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n";
print OUTPUTchr10 "##fileformat=VCFv4.0\n" . "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n";
print OUTPUTchr11 "##fileformat=VCFv4.0\n" . "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n";
print OUTPUTchr12 "##fileformat=VCFv4.0\n" . "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n";
print OUTPUTchr13 "##fileformat=VCFv4.0\n" . "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n";
print OUTPUTchr14 "##fileformat=VCFv4.0\n" . "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n";
print OUTPUTchr15 "##fileformat=VCFv4.0\n" . "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n";
print OUTPUTchr16 "##fileformat=VCFv4.0\n" . "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n";
print OUTPUTchr17 "##fileformat=VCFv4.0\n" . "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n";
print OUTPUTchr18 "##fileformat=VCFv4.0\n" . "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n";
print OUTPUTchr19 "##fileformat=VCFv4.0\n" . "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n";
print OUTPUTchr20 "##fileformat=VCFv4.0\n" . "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n";
print OUTPUTchr21 "##fileformat=VCFv4.0\n" . "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n";
print OUTPUTchr22 "##fileformat=VCFv4.0\n" . "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n";
print OUTPUTchr23 "##fileformat=VCFv4.0\n" . "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n";
print OUTPUTchr24 "##fileformat=VCFv4.0\n" . "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n";
print OUTPUTchr25 "##fileformat=VCFv4.0\n" . "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n";

open (indelVCF, "</data/p255198/resources/gatk_resources/resources/1kg.pilot_release.merged.indels.sites.hg19.vcf");
while ($lines=<indelVCF>){
    chomp $lines;
    $linecount++;
    @array = split('\t', $lines);
    $length = @array;
    my $chr = $array[0];
    $chr =~ s/chr//mg;
    $chr =~ s/X/23/mg;
    $chr =~ s/Y/24/mg;
    $chr =~ s/M/25/mg;
    if ($chr =~ m/#.+?/gs){
        #skip line
    }
    elsif ($chr == 1){
        print OUTPUTchr1 $lines . "\n";
    }
    elsif ($chr == 2){
        print OUTPUTchr2 $lines . "\n";
    }
    elsif ($chr == 3){
        print OUTPUTchr3 $lines . "\n";
    }
    elsif ($chr == 4){
        print OUTPUTchr4 $lines . "\n";
    }
    elsif ($chr == 5){
        print OUTPUTchr5 $lines . "\n";
    }
    elsif ($chr == 6){
        print OUTPUTchr6 $lines . "\n";
    }
    elsif ($chr == 7){
        print OUTPUTchr7 $lines . "\n";
    }
    elsif ($chr == 8){
        print OUTPUTchr8 $lines . "\n";
    }
    elsif ($chr == 9){
        print OUTPUTchr9 $lines . "\n";
    }
    elsif ($chr == 10){
        print OUTPUTchr10 $lines . "\n";
    }
    elsif ($chr == 11){
        print OUTPUTchr11 $lines . "\n";
    }
    elsif ($chr == 12){
        print OUTPUTchr12 $lines . "\n";
    }
    elsif ($chr == 13){
        print OUTPUTchr13 $lines . "\n";
    }
    elsif ($chr == 14){
        print OUTPUTchr14 $lines . "\n";
    }
    elsif ($chr == 15){
        print OUTPUTchr15 $lines . "\n";
    }
    elsif ($chr == 16){
        print OUTPUTchr16 $lines . "\n";
    }
    elsif ($chr == 17){
        print OUTPUTchr17 $lines . "\n";
    }
    elsif ($chr == 18){
        print OUTPUTchr18 $lines . "\n";
    }
    elsif ($chr == 19){
        print OUTPUTchr19 $lines . "\n";
    }
    elsif ($chr == 20){
        print OUTPUTchr20 $lines . "\n";
    }
    elsif ($chr == 21){
        print OUTPUTchr21 $lines . "\n";
    }
    elsif ($chr == 22){
        print OUTPUTchr22 $lines . "\n";
    }
    elsif ($chr == 23){
        print OUTPUTchr23 $lines . "\n";
    }
    elsif ($chr == 24){
        print OUTPUTchr24 $lines . "\n";
    }
    elsif ($chr == 25){
        print OUTPUTchr25 $lines . "\n";
    }
    else{
        print "Error:" . $lines . "\n";
    }
}