#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use Benchmark;
use File::Glob ':glob';

#############################
##USAGE: perl collect_bamIndex_stats.pl <working_directory> <outputfile_prefix>
#############################

my $input_dir = $ARGV[0];
my $output_dir = $ARGV[1];
my $avg_readlength = 90;
chomp $input_dir;
chomp $output_dir;
my $line;
my $lines;
my @results;
my $elements;

##Retrieve files containing QC information
#/data/gcc/projects/gonl/results/first_batch
my @bamindexfiles = glob "$input_dir/*/*BamIndexStats";

open (OUTPUT, ">$output_dir.bamIndex_stats.txt");

while ($_ = shift @bamindexfiles) {
    chomp $_;
    my @chrs;
    if ($_ =~ m/([A-Za-z0-9]{1,6})\/([A-Za-z0-9]{1,10}).human(.+).(.+).BamIndexStats/gs){
        my $sample = $1;
        print OUTPUT "$sample";
        #print "SAMPLE\t";
        open (BAMINDEXFILE, "<$_") or die "Can't open $_$!";
        while ($lines = <BAMINDEXFILE>){
            if ($lines =~ m/(.+) length=\t([0-9]{1,}).+Aligned= ([0-9]{1,}).+/gs){
                my $chr = $1;
                my $chrlength = $2;
                my $aligned = ($3 * $avg_readlength);
                my $result = ($aligned/$chrlength);
                print OUTPUT "\t$result";
            }
        }
    }print OUTPUT "\n";
}