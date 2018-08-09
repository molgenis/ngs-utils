#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use List::Util qw(first);
use POSIX;

my ($help, $inputbed, $output, $outputfolder);

#### get options
GetOptions(
                "h"                             => \$help,
                "input=s"                       => \$inputbed,
                "output=s"                      => \$output,
		"outputfolder=s"                => \$outputfolder,
          );
usage() and exit(1) if $help;
# mandatory args
usage() and exit(1) unless $inputbed;
usage() and exit(1) unless $output;
usage() and exit(1) unless $outputfolder;

chomp $inputbed;
chomp $output;
chomp $outputfolder;


#Open input and output files
open (INPUT, "<", $inputbed ) or die $!;

#Read bed file
my $number=0;
my $regnum = 1;
my $binSize= 1;
while (my $lines=<INPUT>){
        chomp $lines;
        if ($lines !~ m/^track.+/gs) {
                #print $lines . "\n";
                #Remove chr before chrNumber and substitute M with MT
                $lines =~ s/^chr//i;
		$lines =~ s/^M\t/MT\t/i;
                #Split line
                my @array = split("\t", $lines);
                my $chr = $array[0];
                my $start = $array[1];
                my $stop = $array[2];
		my $gene = $array[3];
                my $region = ($stop-$start);
                #Iterate over region and create bins
		open (OUTPUT, ">>", "$outputfolder/$output.per_base.intervals" ) or die $!;
                for (my $i=($start+1); $i<=$stop; $i=($i+$binSize)){
                    print OUTPUT "$chr\t" . $i . "\t" . $i . "\t$gene\n";
                }
                close(OUTPUT);

                $regnum++; 
        }else{
                #Negative check
                #print "$lines\n";
        }
}
sub usage {
        print <<EOF;
#########################################################################################
This script splits a bed file in regions of a length specified by the user.
#########################################################################################
Usage: ./create_per_base_intervals.pl
\t-input\t\t\tInput bed file.
\t-output\t\t\tOutput prefix
\t-outputfolder\t\tOutputfolder	
Example usage: perl create_per_base_intervals.pl -input target_exons.bed -output exonIntervals
#########################################################################################
EOF

}
