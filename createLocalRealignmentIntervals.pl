#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use Cwd;

####################################################################################################
##USAGE: createLocalRealignmentIntervals.pl <index.fa> <dbsnp.rod> <indels.vcf> <output.intervals>##
####################################################################################################

##Variables##
my $index;
my $dbsnp;
my $indels;
my $file;
my $gatk_dir;
my $command;

##RETRIEVE INPUT PARAMETERS##
$index = $ARGV[0];
$dbsnp = $ARGV[1];
$indels = $ARGV[2];
$file = $ARGV[3];
chomp $index;
chomp $dbsnp;
chomp $indels;
chomp $file;

##GATK directory##
$gatk_dir = "data/fvandijk/tools";

##Create indel file##
$command = "java -jar $gatk_dir/Sting/dist/GenomeAnalysisTK.jar -T RealignerTargetCreator" .
    " -R $index" .
    " -D $dbsnp" .
    " -B:indels,VCF $indels" .
    " -o $file";
system($command);
#print $command . "\n";