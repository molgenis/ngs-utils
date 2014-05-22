#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;

############################################################
##USAGE: createBWAandSAMtoolsIndices.pl <path to .fa file>##
############################################################

##Variables
my $file;
my $bwa_dir;
my $samtools_dir;
my $bwa_index;
my $samtools_index;

##RETRIEVE INPUT PARAMETERS##
$file = $ARGV[0];
chomp $file;
$file =~ s/.fai//mg;

##BWA and SAMtools directory##
$bwa_dir = "/data/fvandijk/tools/bwa_45_patched";
$samtools_dir="/data/fvandijk/tools/samtools-0.1.8";

##Create indices##
$bwa_index="$bwa_dir/bwa index -p $file -a is $file";
system($bwa_index); #Execute command
    
$samtools_index="$samtools_dir/samtools faidx $file";
system($samtools_index);

##Check to see if command generated is correct##
#print $bwa_index . "\n";
#print $samtools_index . "\n";