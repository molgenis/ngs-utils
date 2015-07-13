#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use File::Glob ':glob';

my ($help, $input, $output, $filter, $format);

#### get options
GetOptions(
                "h"             => \$help,
                "vcf=s"         => \$input,
                "output=s"      => \$output,
                "filter=s"      => \$filter,
                "format=s"      => \$format
          );
usage() and exit(0) if $help;
# mandatory args
usage() and exit(0) unless $input;
usage() and exit(0) unless $output;
usage() and exit(0) unless $filter;
usage() and exit(0) unless $format;

print "Starting conversion...\n\n";

#Open input and output file
open (INPUT, "<", $input) or die $!;
open (OUTPUT, ">", $output) or die $!;

print "Output written to: $output\n\n";

#Print header line for tabdelim file
print OUTPUT "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\t";
#print OUTPUT "Chromosome\tPosition\tdbSNPid\tReference_allele\tAlternative_allele\tQuality\tFilter\t";

#Substitute commas from cmdline input with tabs and add to header line
chomp $filter;
my @filterline = split(",",$filter);
$filter =~ s/,/\t/gs;
print OUTPUT "$filter\t";

#Substitute commas from cmdline input with tabs for FORMAT line
chomp $format;
my @formatline = split(",",$format);
$format =~ s/,/\t/gs;
print OUTPUT "$format\n";

#Read VCF file
my $line;
while ($line=<INPUT>){
    #Read SNPs from VCF file
    chomp $line;
    if ($line !~ m/^#.+/gs){
        #print "$line\n";
        my @full_line = split("\t", $line);
        my $chr = $full_line[0];
        my $pos = $full_line[1];
        my $id = $full_line[2];
        my $ref = $full_line[3];
        my $alt = $full_line[4];
        my $qual = $full_line[5];
        my $filter = $full_line[6];
        print OUTPUT "$chr\t$pos\t$id\t$ref\t$alt\t$qual\t$filter";
        my $info = $full_line[7];
        my $format = $full_line[8];
        my $sampleinfo = $full_line[9];
        my @info_line = split(";", $info);
        my @format_line = split(":", $format);
        my @sample_line = split(":", $sampleinfo);
        foreach my $input_elem (@filterline){
            my $count = 0;
            my $value = "N/A";
            #print "$input_elem BLAAT";
            foreach my $info_elem (@info_line){
                #$count++;
                #print "$info_elem\t$input_elem\n";
                if ($info_elem =~ m/$input_elem\=(.+)/gs){
                    $value = $1;
                }
            }
            print OUTPUT "\t$value";
        }
        #Filter on format and sample field
        foreach my $format_elem (@formatline){
            my $search = "$format_elem";
            my $formatvalue = "N/A";
            my %index;
            @index{@format_line} = (0..$#format_line);
            my $index = $index{$search};
            $formatvalue = $sample_line[$index];
            print OUTPUT "\t$formatvalue";
        }
    print OUTPUT "\n";
    }
}

print "Finished conversion!\n";

sub usage {
        print <<EOF;
##############################################################################################
This script takes a VCF file containing ONE sample as input and produces a tab separated
output file. The user can provide a list of annotations (words) as found in the INFO and
FORMAT column from which the script filters the values.
For questions please contact: freerk.van.dijk\@gmail.com
##############################################################################################
Usage: ./vcf2tab.pl
\t-vcf          The input VCF file to filter, this file should contain only ONE sample.
\t-output       The output file. (tab separated)
\t-filter       The annotations to filter on, comma separated list.
                      Example: AC,AF,AN,BaseCounts
\t-format       The annotations in the format field to filter on, comma separated list.
                      Example: GT,DP
##############################################################################################
EOF
 
}
