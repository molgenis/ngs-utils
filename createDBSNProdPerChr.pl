#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;

#####################################################
##USAGE: createDBSNProdPerChr.pl <dbsnprod.chr.rod>##
#####################################################

##Variables
my $lines;
my $file;
my $dbsnpbuild;
my $chr_num;
my $genome;
my $ncbibuild;
my $chr;
my @array;

##RETRIEVE INPUT PARAMETERS##
$file = $ARGV[0];
chomp $file;

##Open filehandlers##
open (OUTPUT, ">$file") or die $!;
#open (OUTPUT, ">/data/fvandijk/resources/$genome/dbsnp/dbsnp_$dbsnpbuild" . "_b" . $ncbibuild . "_chr" . $chr_num . ".rod") or die $!;

##Retrieve genome, dbsnpbuild, ncbibuild and chr number##
if ($file =~ m/.*\/resources\/(hg[0-9]{1,2})\/dbsnp\/dbsnp_([0-9]{1,3})_b([0-9]{1,2})_chr([0-9]{1,2}).rod/gs){
    $genome = $1;
    $dbsnpbuild = $2;
    $ncbibuild = $3;
    $chr_num = $4;
}


##Open correct rod file corresponding to genome, dbsnpbuild and chr##
open (dbSNProd, "</data/fvandijk/resources/$genome/dbsnp/dbsnp_$dbsnpbuild" . "_b" . $ncbibuild . ".rod") or die $!;

##Search per line if chr number is in file, if so, write to output##
while ($lines=<dbSNProd>){
    chomp $lines;
    @array = split('\t', $lines);
    $chr = $array[1];
    $chr =~ s/chr//mg;
    $chr =~ s/X/23/mg;
    $chr =~ s/Y/24/mg;
    $chr =~ s/M/25/mg;
    if ($chr == $chr_num){
        print OUTPUT $lines . "\n";
    }else{
        #skip line
    }
}