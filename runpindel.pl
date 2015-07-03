#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use File::Glob ':glob';

#############################
##USAGE: perl runpindel.pl <path_to_pindel> <reference_genome.fa> <path_to_output4pindel.txt> <output> <path_to_breakdancer_result> <indel_size> <cores??> <chromosome>
#############################

##Retrieve input parameters
my $path_to_pindel = $ARGV[0];
my $ref_genome = $ARGV[1];
my $output4pindel = $ARGV[2];
my $outputdir = $ARGV[3];
my $brkdncresult = $ARGV[4];
my $indel_size = $ARGV[5];
my $cores = $ARGV[6];

##Variables
my $line;
my @chromosomes;
my $chr;
my $execute_pindel;
my $chrom;

##Read reference genome
open (REFGENOME, "<$ref_genome");

##Extract chromosomes from reference genome fasta file
while ($line=<REFGENOME>){
    chomp $line;
    if ($line =~ m/>(.+)\sdna:.+/gs || $line =~ m/>(.+)\sgi.+/gs || $line =~ m/>([0-9A-Za-z]{1,2})/gs){
#        $chr = $1;
        print "$chr" . "\n";
        push (@chromosomes, $chr);                          #Push chromosomes in array
    }
}

##Execute statement for pindel
#$execute_pindel = "$path_to_pindel $ref_genome $output4pindel $outputdir $brkdncresult $indel_size $cores $chrom";

##Execute pindel foreach chromosome in array
foreach $chrom (@chromosomes){
    $execute_pindel = "$path_to_pindel $ref_genome $output4pindel $outputdir $brkdncresult $indel_size $cores $chrom";
    system($execute_pindel);
}