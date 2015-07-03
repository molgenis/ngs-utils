#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use File::Glob ':glob';
use Getopt::Long;
use List::Util qw(first);
use POSIX;



#############################
##USAGE: perl create_interval_list.pl -Ref <FILE> -Exons <FILE> -Baits <FILE> -Strand TRUE|FALSE
#############################

my ($help, $HumanRef, $ExonBed, $BaitsBed);
my $Strand = "FALSE";

#### get options
GetOptions(
                "h"                             => \$help,
                "Ref=s"                      => \$HumanRef,
                "Exons=s"                       => \$ExonBed,
                "Baits=s"                       => \$BaitsBed,
                "Strand=s"                       => \$Strand,
                );
usage() and exit(1) if $help;
# mandatory args
usage() and exit(1) unless $HumanRef;
usage() and exit(1) unless $ExonBed;
usage() and exit(1) unless $BaitsBed;

chomp $HumanRef;
chomp $ExonBed;
chomp $BaitsBed;


#Variables#
my $dictlines = $HumanRef;
my $exonlines = $ExonBed;
my @exonarray;
my $baitlines = $BaitsBed;
my @baitarray;

#Filehandlers#
print $ExonBed . "\n";
print $BaitsBed . "\n";
print $HumanRef . "\n";

open (EXONS,"<$ExonBed.bed") or die "Can't open $_$!";
open (BAITS,"<$BaitsBed.bed") or die "Can't open $_$!";
open (DICT,"<$HumanRef") or die "Can't open $_$!";
open (EXONS_OUT,">$ExonBed.interval_list") or die "Can't open $_$!";
open (BAITS_OUT,">$BaitsBed.interval_list") or die "Can't open $_$!";


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
    if ($exonlines =~ m/[0-9XYMT]{1,2}.+/gs){
        @exonarray = split('\t', $exonlines);#Split lines by tab
        my $chrom = $exonarray[0];
        my $start = $exonarray[1];
        my $end = $exonarray[2];
        my $target = $exonarray[3];
        my $strandexon = "";

	if ($Strand =~ /TRUE|True|true/){
        	$strandexon = $exonarray[5];
        	}
        	elsif ($Strand =~ /FALSE|False|false/){
        	$strandexon = "+";
        	}
		else{
		usage() and exit(1)
		}

        $chrom =~ s/chr//mg;#Remove chr, Un_, _random from chromosome name
        print EXONS_OUT "$chrom\t$start\t$end\t$strandexon"."\t$target\n";#Print chromosome and other information to output file
    }
}

#Read Agilent Bait file#
while ($baitlines=<BAITS>){
    chomp $baitlines;
    if ($baitlines =~ m/[0-9XYMT]{1,2}.+/gs){
        @baitarray = split('\t', $baitlines);#Split lines by tab
        my $chrbait = $baitarray[0];
        my $startbait = $baitarray[1];
        my $endbait = $baitarray[2];
        my $targetbait = $baitarray[3];
	my $strandbait = "";

        if ($Strand =~ /TRUE|True|true/){
        	$strandbait = $baitarray[5];
        	}
		elsif ($Strand =~ /FALSE|False|false/){
                $strandbait = "+";
                }
        	else{
		usage() and exit(1)
        	}

        $chrbait =~ s/chr//mg;#Remove chr from chromosome name
        print BAITS_OUT "$chrbait\t$startbait\t$endbait\t$strandbait"."\t$targetbait\n";#Print output to file
    }
}

sub usage {
        print <<EOF;
#########################################################################################
This script merges two coverage interval summary files for one sample and creates a *.bed
file with coverage profiles, which can be displayed in the UCSC Genome Browser.
#########################################################################################
Usage: ./create_interval_listv4.pl
\t-Ref                 Human referencee. (Example: /gcc/resources/b37/indices/human_g1k_v37.dict)
\t-Exons                    Exons Bedfile.
\t-Baits                    Baits Bedfile.

(optional)

\t-Strand                   Is a stand present in the inputBEDfile? TRUE|FALSE (Default FALSE)

output: filename ".interval_list" is added as extention for both inputfiles
#########################################################################################
EOF
 
}
