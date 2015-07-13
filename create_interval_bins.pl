#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use List::Util qw(first);
use POSIX;

my ($help, $dict, $binSize, $inputbed, $output);



#### get options
GetOptions(
                "h"                             => \$help,
                "dict=s"                        => \$dict,
                "binSize=s"                     => \$binSize,
                "input=s"                       => \$inputbed,
                "output=s"                      => \$output,
          );
usage() and exit(1) if $help;
# mandatory args
usage() and exit(1) unless $dict;
usage() and exit(1) unless $binSize;
usage() and exit(1) unless $inputbed;
usage() and exit(1) unless $output;

chomp $dict;
chomp $inputbed;
chomp $output;
chomp $binSize;

#Copy dict file to OUTPUT
`cp $dict $output.1.interval_list`;
`cp $dict $output.2.interval_list`;

#Open input and output files
open (INPUT, "<", $inputbed ) or die $!;



#Read bed file
my $number=0;
my $regnum = 1;
while (my $lines=<INPUT>){
        chomp $lines;
	if ($lines eq ""){
		
	}
        elsif ($lines !~ m/^track.+/gs) {
                #print $lines . "\n";
                #Remove chr before chrNumber and substitute M with MT
                $lines =~ s/^chr//gs;
                $lines =~ s/^M/MT/gs;
                my $idx = 1;
                #Split line
                my @array = split("\t", $lines);
                my $chr = $array[0];
                my $start = $array[1];
                my $stop = $array[2];
                my $name = $array[3];
                my $region = ($stop-$start);
                #Iterate over region and create bins
                for (my $i=$start; $i<$stop; $i+=$binSize){
                    #Check if number is even, else write away in file2
                    if ($number % 2 == 0) {
                        open (OUTPUT, ">>", "$output.1.interval_list" ) or die $!;
                    }else{
                        open (OUTPUT, ">>", "$output.2.interval_list" ) or die $!;
                    }
                    #Check for all, but first line
                    if ($idx>1){
                        if (($i+$binSize+$binSize+$idx)>$stop){
                            print OUTPUT "$chr\t" . ($i+$idx-1) . "\t$stop\t+\tregion$regnum\_chunk$idx\n";
                            last;
                            $number++;
                        }else {
                            print OUTPUT "$chr\t" . ($i+$idx-1) . "\t" . ($i+$binSize+$idx-1) . "\t+\tregion$regnum\_chunk$idx\n";
                            $idx++;
                            $number++;
                        }
                        
                    }else {
                    	if (($i+$binSize+$binSize+$idx)>$stop){
                        	print OUTPUT "$chr\t" . ($i+$idx-1) . "\t$stop\t+\tregion$regnum\_chunk$idx\n";
                            last;
                            $number++;
                            }
                            else {
                        	                      	
                        	print OUTPUT "$chr\t$i\t" . ($i+$binSize) . "\t+\tregion$regnum\_chunk$idx\n";
                        	$idx++;
                        	$number++;
                        	}
                    }
                    #$number++;
                    close(OUTPUT);
                }
                $regnum++; 
        }else{
                #Negative check
                #print "$lines\n";
        }
}



sub usage {
        print <<EOF;
#########################################################################################
This script splits a bed file in regions of a length specified by the user, writing away
each bin to two seperate files by taking turns.
#########################################################################################
Usage: ./create_interval_bins.pl
\t-dict                     Fasta reference dictory file (*.dict).
\t-binSize                  Size of output bins.
\t-input                    Input bed file.
\t-output                   Output prefix, resulting in <output>.1.interval_list and
                                                              <output>.2.interval_list
Example usage: perl create_bins.pl -dict human_g1k_v37.dict -binSize 49 \\
-input target_exons.bed -output exonIntervals
#########################################################################################
EOF
 
}