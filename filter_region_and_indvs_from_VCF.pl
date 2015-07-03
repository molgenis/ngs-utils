#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use List::Util qw(first);

my ($help, $inputvcf, $outputvcf, $indv, $region);

#### get options
GetOptions(
                "h"                         => \$help,
                "inputVCF=s"                => \$inputvcf,
                "outputVCF=s"               => \$outputvcf,
                "individuals=s"             => \$indv,
                "region=s"                  => \$region  
          );
usage() and exit(1) if $help;
# mandatory args
usage() and exit(1) unless $inputvcf;
usage() and exit(1) unless $outputvcf;
usage() and exit(1) unless $indv;
usage() and exit(1) unless $region;

chomp $inputvcf;
chomp $outputvcf;
chomp $indv;
chomp $region;

#Open input and output files
open (INPUT, "<", $inputvcf ) or die $!;
open (OUTPUT, ">", $outputvcf ) or die $!;

my $chr;
my $first;
my $last;
if ($region =~ m/(.+):(.+)-(.+)/gs){
    $chr = $1;
    $chr =~ s/chr//gs;
    $first = $2;
    $last = $3;
}else{
    print "ERROR: Specified region has incorrect format!\n";
    exit(1);
}

my $lines;
my @columnindices;

while ($lines = <INPUT>){
    chomp $lines;
    if ($lines =~ m/^##/gs) {
        #Print header lines
        print OUTPUT "$lines\n";
    }elsif ($lines =~ m/^#CHROM.+/gs){
        #Sample line, extract indices of parents
        my @array = split("\t", $lines);
        print OUTPUT $array[0];
        for (my $i=1; $i<=8; $i++){
            push(@columnindices, $i);
            print OUTPUT "\t" . $array[$i];
        }
        foreach my $elem (@array){
            #$index++;
            if ($elem =~ m/$indv/gs) {
                my $index = first { $array[$_] eq $elem } 0..$#array;
                push(@columnindices, $index);
                print OUTPUT "\t" . $elem;
            }   
        }
        print OUTPUT "\n";
        shift(@columnindices);
    }else{
        #Searchline
        my @array = split("\t", $lines);
        my $arrchr = $array[0];
        $arrchr =~ s/chr//gs;
        my $pos = $array[1];
        if ("$arrchr" eq "$chr" && $pos >= $first && $pos <= $last){
            print OUTPUT "$arrchr\t$pos";
            foreach my $ele (@columnindices){
                print OUTPUT "\t" . $array[$ele];
            }
            print OUTPUT "\n";
        }
    }
    #print OUTPUT "\n";
}



sub usage {
        print <<EOF;
#########################################################################################
This script filters a VCF file based on individuals (devined by a regular expression) and
a region of interest.
#########################################################################################
Usage: ./filter_region_and_indvs_from_VCF.pl
\t-inputVCF                 Input VCF file to extract individuals from.
\t-outputVCF                Output VCF file.
\t-individuals              Regular expression to filter on, between quotes.
                                    Example: To filter on only GoNL parents use
                                    "[A-Za-z]{1,1}[0-9]{1,3}[ABab]"
\t-region                   Region to filter (chr:start-stop).
                                    Example: 21:44471301-44498472
#########################################################################################
EOF
 
}