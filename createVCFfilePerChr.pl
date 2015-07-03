#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;

####################################################
##USAGE: createVCFfilePerChr.pl <VCF file per chr>##
####################################################

###Only the first line describing the VCF version and the line containing column headers should be kept
###To remove incorrect header from VCF files use: sed '2,25d' <filename> > <new file>

##Variables##
my $file;
my $genome;
my $chr_num;
my $lines;
my @array;

##RETRIEVE INPUT PARAMETERS##
$file = $ARGV[0];
chomp $file;

##Create filehandlers##
open (OUTPUT, ">$file") or die $!;

##Check for used genome build (hg18/hg19)##
if ($file =~ m/.*\/resources\/(hg[0-9]{1,2})\/indels\/1kg.*chr([0-9]{1,2}).vcf/gs){
    $genome = $1;
    $chr_num = $2;
    #print "genome $genome and chr is $chr_num\n";
}

##Open VCF file corresponding to right build name##
open (INDELVCF, "</data/fvandijk/resources/$genome/indels/1kg.pilot_release.merged.indels.sites.$genome.vcf") or die $!;

##Print VCF header into output##
print OUTPUT "##fileformat=VCFv4.0\n" . "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n";

##Read line, determine if chr of interest is in it, if so, write to output file##
while ($lines=<INDELVCF>){
    chomp $lines;
    if ($lines =~ m/#.+/gs){
        #skip line
    }else{
        @array = split('\t', $lines);
        my $chr = $array[0];
        $chr =~ s/chr//mg;
        $chr =~ s/X/23/mg;
        $chr =~ s/Y/24/mg;
        $chr =~ s/M/25/mg;
        if ($chr == $chr_num){
            print OUTPUT $lines . "\n";
        }
    }
}