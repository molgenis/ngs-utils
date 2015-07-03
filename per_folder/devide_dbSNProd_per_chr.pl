#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use Cwd;

my $lines;
my @array;
my $element;
my $length;
my @snps;
my $count = 0;
my $linecount = 0;

open (OUTPUTchr1, ">/data/p255198/resources/gatk_resources/resources/dbsnp_129_b37_chr1.rod");
open (OUTPUTchr2, ">/data/p255198/resources/gatk_resources/resources/dbsnp_129_b37_chr2.rod");
open (OUTPUTchr3, ">/data/p255198/resources/gatk_resources/resources/dbsnp_129_b37_chr3.rod");
open (OUTPUTchr4, ">/data/p255198/resources/gatk_resources/resources/dbsnp_129_b37_chr4.rod");
open (OUTPUTchr5, ">/data/p255198/resources/gatk_resources/resources/dbsnp_129_b37_chr5.rod");
open (OUTPUTchr6, ">/data/p255198/resources/gatk_resources/resources/dbsnp_129_b37_chr6.rod");
open (OUTPUTchr7, ">/data/p255198/resources/gatk_resources/resources/dbsnp_129_b37_chr7.rod");
open (OUTPUTchr8, ">/data/p255198/resources/gatk_resources/resources/dbsnp_129_b37_chr8.rod");
open (OUTPUTchr9, ">/data/p255198/resources/gatk_resources/resources/dbsnp_129_b37_chr9.rod");
open (OUTPUTchr10, ">/data/p255198/resources/gatk_resources/resources/dbsnp_129_b37_chr10.rod");
open (OUTPUTchr11, ">/data/p255198/resources/gatk_resources/resources/dbsnp_129_b37_chr11.rod");
open (OUTPUTchr12, ">/data/p255198/resources/gatk_resources/resources/dbsnp_129_b37_chr12.rod");
open (OUTPUTchr13, ">/data/p255198/resources/gatk_resources/resources/dbsnp_129_b37_chr13.rod");
open (OUTPUTchr14, ">/data/p255198/resources/gatk_resources/resources/dbsnp_129_b37_chr14.rod");
open (OUTPUTchr15, ">/data/p255198/resources/gatk_resources/resources/dbsnp_129_b37_chr15.rod");
open (OUTPUTchr16, ">/data/p255198/resources/gatk_resources/resources/dbsnp_129_b37_chr16.rod");
open (OUTPUTchr17, ">/data/p255198/resources/gatk_resources/resources/dbsnp_129_b37_chr17.rod");
open (OUTPUTchr18, ">/data/p255198/resources/gatk_resources/resources/dbsnp_129_b37_chr18.rod");
open (OUTPUTchr19, ">/data/p255198/resources/gatk_resources/resources/dbsnp_129_b37_chr19.rod");
open (OUTPUTchr20, ">/data/p255198/resources/gatk_resources/resources/dbsnp_129_b37_chr20.rod");
open (OUTPUTchr21, ">/data/p255198/resources/gatk_resources/resources/dbsnp_129_b37_chr21.rod");
open (OUTPUTchr22, ">/data/p255198/resources/gatk_resources/resources/dbsnp_129_b37_chr22.rod");
open (OUTPUTchr23, ">/data/p255198/resources/gatk_resources/resources/dbsnp_129_b37_chrX.rod");
open (OUTPUTchr24, ">/data/p255198/resources/gatk_resources/resources/dbsnp_129_b37_chrY.rod");
open (OUTPUTchr25, ">/data/p255198/resources/gatk_resources/resources/dbsnp_129_b37_chrM.rod");

open (dbSNProd, "</data/p255198/resources/gatk_resources/resources/dbsnp_129_b37.rod");
while ($lines=<dbSNProd>){
    chomp $lines;
    $linecount++;
    @array = split('\t', $lines);
    $length = @array;
    my $nm = $array[0];
    my $chr = $array[1];
    $chr =~ s/chr//mg;
    $chr =~ s/X/23/mg;
    $chr =~ s/Y/24/mg;
    $chr =~ s/M/25/mg;
    if ($chr == 1){
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