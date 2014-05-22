#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use Getopt::Long;


my ($help, $input, $output, $filtervalue);
my $count = 0;

#### get options
GetOptions(
                "h"                     => \$help,
                "inputBeagleR2file=s"   => \$input,
                "outputSnpIDfile=s"     => \$output,
                "R2threshold=s"                  => \$filtervalue
          );
usage() and exit(1) if $help;
# mandatory args
usage() and exit(1) unless $input;
usage() and exit(1) unless $output;
usage() and exit(1) unless $filtervalue;


open (INPUT, "<", $input) or die $!;
open (OUTPUT, ">", $output) or die $!;

print "\nStarting filtering..\n\n";

my $lines;
while ($lines = <INPUT>){
    chomp $lines;
    my @array = split("\t", $lines);
    my $snpID = $array[1];
    my $pos = $array[2];
    my $refall = $array[3];
    my $altall = $array[4];
    my $r2value = $array[5];
    $r2value =~ s/NaN/0.00/g;
    if ($r2value >= $filtervalue){
        print OUTPUT "$snpID\n";
        $count++;
    }
}   



print "Filtering finished!\n";
print "Number of SNPs left: $count\n";

close(INPUT);
close(OUTPUT);


sub usage {
        print <<EOF;
#########################################################################################
filter_BeagleR2 (version 1)

This script filters a beagleR2 file on user specified R2 value. Output is a file
containing a list of snpIDs.
#########################################################################################
Usage: ./convert_snpIDs.pl
\t-inputBeagleR2file    Input beagleR2 file to filter
\t-outputSnpIDfile      Output file
\t-R2threshold          R2 threshold value to filter on (>=)
#########################################################################################
EOF

}