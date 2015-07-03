#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use File::Glob ':glob';

#############################
##USAGE: perl collect_GoNL_sample_coverage_stats.pl <working_directory> <outputfile_prefix>
#############################

my $input_dir = $ARGV[0];
my $output_dir = $ARGV[1];
chomp $input_dir;
chomp $output_dir;
my $line;
my $lines;

##Retrieve files containing QC information
#/data/gcc/projects/gonl/results/first_batch
my @coverage_prop = glob "$input_dir/*/*.coverage.*proportions";
my @coverage_counts = glob "$input_dir/*/*.coverage.*counts";

open (OUTPUTPROP, ">$output_dir.proportions.txt");
open (OUTPUTCOUNTS, ">$output_dir.counts.txt");

while ($_ = shift @coverage_prop) {
    #print $_ . "\n";
    chomp $_;
    open (COVPROP, "<$_") or die "Can't open $_$!";
    while ($line = <COVPROP>){
        chomp $line;
        if ($line =~ m/[A-Za-z0-9]{1,10}.+1\.00.+/gs){
            print OUTPUTPROP "$line\n";
        }
    }
}

while ($_ = shift @coverage_counts) {
    chomp $_;
    if ($_ =~ m/.+\/([A-Za-z0-9]{1,}).vc03.(.+).sample_cumulative_coverage_counts/gs){
        my $sample = $1;
        open (COVCOUNTS, "<$_") or die "Can't open $_$!";
        while ($lines = <COVCOUNTS>){
            chomp $lines;
            if ($lines =~ m/(NSamples_1)\t(.+)/gs){
                my $counts = $2;
                print OUTPUTCOUNTS "$sample\t$counts\n";
            }
        }
    }
}