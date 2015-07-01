#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use File::Glob ':glob';


#############################
##USAGE: perl create_interval_list.pl
#############################

#Variables#
my $dictlines;
my $exonlines;
my @exonarray;
my $baitlines;
my @baitarray;

#Filehandlers#
open (EXONS,"</gcc/resources/b37/intervals/Sureselect_MP_Capture_Library_Design_0679001_exons_b37_human_g1k_v37.bed") or die "Can't open $_$!";
open (BAITS,"</gcc/resources/b37/intervals/Sureselect_MP_Capture_Library_Design_0679001_baits_b37_human_g1k_v37.bed") or die "Can't open $_$!";
#open (BAITS,"</data/fvandijk/resources/hg18/intervals/test.bed") or die "Can't open $_$!";
open (DICT,"</gcc/resources/b37/indices/human_g1k_v37.dict") or die "Can't open $_$!";
open (EXONS_OUT,">/gcc/resources/b37/intervals/Sureselect_MP_Capture_Library_Design_0679001_exons_b37_human_g1k_v37.interval_list") or die "Can't open $_$!";
open (BAITS_OUT,">/gcc/resources/b37/intervals/Sureselect_MP_Capture_Library_Design_0679001_baits_b37_human_g1k_v37.interval_list") or die "Can't open $_$!";

#Read .dict file of used reference genome and extract the @HD and @SQ lines. Also change Sort Order#
while ($dictlines=<DICT>){
    chomp $dictlines;
    if ($dictlines =~ m/\@HD.+/gs || $dictlines =~ m/\@SQ.+/gs){
        $dictlines =~ s/SO:unsorted/SO:coordinate/mg;
        print EXONS_OUT $dictlines . "\n";
        print BAITS_OUT $dictlines . "\n";
    }
}

#Read Agilent Exon file#
while ($exonlines=<EXONS>){
    chomp $exonlines;
    if ($exonlines =~ m/[0-9XYMT]{1,2}.+/gs){
        @exonarray = split('\t', $exonlines);#Split lines by tab
        my $chrom = $exonarray[0];
        my $start = $exonarray[1];
        my $end = $exonarray[2];
        my $target = $exonarray[3];
        #my $strand = $exonarray[5];
        $chrom =~ s/chr//mg;#Remove chr, Un_, _random from chromosome name
        #print EXONS_OUT "$chrom\t$start\t$end\t$strand\t$target\n";#Print chromosome and other information to output file
        print EXONS_OUT "$chrom\t$start\t$end\t" . "+" . "\t$target\n";#Print chromosome and other information to output file
    }
}

#Read Agilent Bait file#
while ($baitlines=<BAITS>){
    chomp $baitlines;
    if ($baitlines =~ m/[0-9XYMT]{1,2}.+/gs){
        @baitarray = split('\t', $baitlines);#Split lines by tab
        my $chrbait = $baitarray[0];
        my $startbait = $baitarray[1];
        my $endbait = $baitarray[2];
        my $targetbait = $baitarray[3];
        #my $strandbait = $baitarray[5];
        $chrbait =~ s/chr//mg;#Remove chr from chromosome name
        print BAITS_OUT "$chrbait\t$startbait\t$endbait\t" . "+" . "\t$targetbait\n";#Print output to file
    }
}
