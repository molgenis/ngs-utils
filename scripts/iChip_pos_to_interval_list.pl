#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use File::Glob ':glob';
use File::stat;
use Time::localtime;

##########################
#### This script reads a VCF, retrieves the
#### positions of all SNPs and converts
#### these to an inteval_list
##########################

#perl iChip_pos_to_interval_list <VCF file> <output.interval_list>

my $input = $ARGV[0];
my $output = $ARGV[1];

open (INPUT, "<$input");
open (OUTPUT, ">$output");

my $line;
while ($line = <INPUT>){
    chomp $line;
    if ($line !~ m/^#.+/gs){
        my @array = split('\t', $line);
        my $chr = $array[0];
        my $pos = $array[1];
        print OUTPUT "$chr:$pos-$pos\n";
    }
}