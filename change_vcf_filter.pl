#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;

###
#Created by FvD: freerk.van.dijk@gmail.com
###

###
#Usage: perl change_vcf_filter.pl <input.vcf> <output.vcf> <DP> <QUALITY>
#Example: perl /data/fvandijk/concordance/change_vcf_filter.pl /data/fvandijk/inhouse_concordance/CD/Cytochip/VCF/0605-041.0411_L8.HSpe22.variant_annotator.ftl.human_g1k_v37.2011_05_25_09_58.snps.variantannotated.vcf /data/fvandijk/inhouse_concordance/CD/Cytochip/VCF/0605-041_q20_dp10.vcf 10 20
###

my $input = $ARGV[0];
my $output = $ARGV[1];
my $depth = $ARGV[2];
my $quality = $ARGV[3];
my $lines;
my @columns;

chomp $input;
chomp $output;
chomp $depth;
chomp $quality;

open (INPUT, "<$input") or die 'cannot open file $!';
open (OUTPUT, ">$output") or die 'cannot open file $!';

while ($lines = <INPUT>){
    chomp $lines;
    if ($lines =~ m/^#.+/gs){
        print OUTPUT "$lines\n";
    }else{
        my @columns = split('\t', $lines);
        my $chr = $columns[0];
        my $pos = $columns[1];
        my $id = $columns[2];
        my $ref = $columns[3];
        my $alt = $columns[4];
        my $qual = $columns[5];
        my $filter = $columns[6];
        my $info = $columns[7];
        my $format = $columns[8];
        my $sample = $columns[9];
        if ($info =~ m/.+DP=([0-9]{1,});.+/gs || $info =~ m/DP=([0-9]{1,});.+/gs || $info =~ m/DP=([0-9]{1,})/gs){
            my $dp = $1;
            #print "$dp\n";
            if ($qual ne "."){
                if ($qual >= $quality && $dp >= $depth){
                    $filter = "PASS";
                }else{
                    $filter = "Filtered";
                }
                print OUTPUT "$chr\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t$info\t$format\t$sample\n";
            }
        }
    }
}