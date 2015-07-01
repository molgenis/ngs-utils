#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use File::Glob ':glob';

#####################################################################
# Changed tool to only create intervals for chromosome 1-22, X and Y#
#####################################################################


#############################
##USAGE: perl create_interval_list.pl
#############################

#Variables#
my $dictlines;
my $exonlines;
my @exonarray;
my $baitlines;
my @baitarray;

print "WARNING: BAITS files not correctly implemented!! Contact Freerk! :-)\n";

#Filehandlers#
open (EXONS,"</target/gpfs2/gcc/home/mdijkstra/selectgenes/microcephaly.bed") or die "Can't open $_$!";
open (BAITS,"</target/gpfs2/gcc/home/mdijkstra/selectgenes/microcephaly.bed") or die "(<- Warning, this file is no BAITS file!! Can't open $_$!";
#open (BAITS,"</data/fvandijk/resources/hg18/intervals/test.bed") or die "Can't open $_$!";
open (DICT,"</target/gpfs2/gcc/resources/hg19/indices/human_g1k_v37.dict") or die "Can't open $_$!";
open (EXONS_OUT,">/target/gpfs2/gcc/home/mdijkstra/selectgenes/microcephaly.interval_list") or die "Can't open $_$!";
open (BAITS_OUT,">/target/gpfs2/gcc/home/mdijkstra/test.interval_list_REMOVE") or die "Can't open $_$!";

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
    if ($exonlines =~ m/chr.+/gs){
        @exonarray = split('\t', $exonlines);#Split lines by tab
        my $chr = $exonarray[0];
        my $start = $exonarray[1];
        my $end = $exonarray[2];
        my $target = $exonarray[3];
        my $strand = $exonarray[5];
        $chr =~ s/chr//mg;#Remove chr, Un_, _random from chromosome name
        #$chr =~ s/Un_//mg;
        #$chr =~ s/_random//mg;
        #$chr =~ s/gl/GL/mg;
        if ($chr !~ m/.+_.+/gs){
        if ($chr =~ m/^[0-9XYMT]{1,2}/gs){
            print EXONS_OUT "$chr\t$start\t$end\t$strand\t$target\n";#Print chromosome and other information to output file
        }
        }
    }
}

#Read Agilent Bait file#
while ($baitlines=<BAITS>){
    chomp $baitlines;
    if ($baitlines =~ m/chr.+/gs){
        @baitarray = split('\t', $baitlines);#Split lines by tab
        my $chrbait = $baitarray[0];
        my $startbait = $baitarray[1];
        my $endbait = $baitarray[2];
        my $targetbait = $baitarray[3];
        my $strandbait = $baitarray[5];
        $chrbait =~ s/chr//mg;#Remove chr, Un_, _random from chromosome name
        #$chrbait =~ s/Un_//mg;
        #$chrbait =~ s/_random//mg;
        #$chrbait =~ s/gl/GL/mg;
        if ($chrbait !~ m/.+_.+/gs){
        if ($chrbait =~ m/^[0-9XYMT]{1,2}/gs){
            print BAITS_OUT "$chrbait\t$startbait\t$endbait\t$strandbait\t$targetbait\n";#Print output to file
        }
        }
    }
}
