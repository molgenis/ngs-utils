#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use Getopt::Long;


my ($help, $input, $output, $delimiter);

#### get options
GetOptions(
                "h"             => \$help,
                "inputvcf=s"    => \$input,
                "outputvcf=s"   => \$output,
                "delimiter=s"   => \$delimiter
          );
usage() and exit(1) if $help;
# mandatory args
usage() and exit(1) unless $input;
usage() and exit(1) unless $output;
usage() and exit(1) unless $delimiter;


open (INPUT, "<", $input) or die $!;
open (OUTPUT, ">", $output) or die $!;

print "\nStarting ID conversion..\n\n";

my $lines;
while ($lines = <INPUT>){
    chomp $lines;
    if ($lines =~ m/^#.+/gs){
        print OUTPUT "$lines\n";
    }else{
        if ($lines =~ m/([0-9XYMT]{1,2})\t([0-9]+)\t([^\t]+)(\t.+)/gs){
            print OUTPUT "$1\t$2\t$1$delimiter$2$4\n";
        }
    }
}

print "Converting finished!\n";

close(INPUT);
close(OUTPUT);


sub usage {
        print <<EOF;
#########################################################################################
convert_snpIDs (version 2)

This script converts the snpIDs in a VCF file to a format used by Impute2. The delimiter
for a snpID needs to be specified as well. 
#########################################################################################
Usage: ./convert_snpIDs.pl
\t-inputvcf     VCF file to convert IDs from
\t-outputvcf    Output VCF containing updated IDs
\t-delimiter    The delimiter to be used in the snpID
#########################################################################################
EOF
 
}
