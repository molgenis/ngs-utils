#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use File::Glob ':glob';

my ($help, $sample, $file1, $file2, $output);

#### get options
GetOptions(
                "h"                             => \$help,
                "sample=s"                      => \$sample,
                "file1=s"                       => \$file1,
                "file2=s"                       => \$file2,
                "output=s"                      => \$output,
          );
usage() and exit(1) if $help;
# mandatory args
usage() and exit(1) unless $sample;
usage() and exit(1) unless $file1;
usage() and exit(1) unless $file2;
usage() and exit(1) unless $output;

chomp $sample;
chomp $file1;
chomp $file2;
chomp $output;

print "\nMerging the following two files:\n";
print "$file1\n";
print "$file2\n\n";
print "Writing output to:\n";
print "$output\n\n";

#Cat file 1 and 2, split regions and sort on chr,start,stop
`cat $file1 | awk 'NR>1' > $output.tmp1.txt ; awk 'NR>1' $file2 >> $output.tmp1.txt`;
`perl -pi -e 's/:/\t/g' $output.tmp1.txt`;
`perl -pi -e 's/-/\t/g' $output.tmp1.txt`;
`awk '{ if ( \$1 ~ /[X]/ ) {\$1="23"; print \$0}else if ( \$1 ~ /[Y]/ ) {\$1="24"; print \$0}else{ print \$0} }' OFS="\t" $output.tmp1.txt | sort -k1n -k2n -k3n > $output.tmp2.txt`;

#Open inputs and output files
open (INPUT, "<", "$output.tmp2.txt" ) or die $!;
open (OUTPUT, ">", "$output" ) or die $!;

#Write header (track line) to output file
print OUTPUT "track name=\"$sample\" description=\"Coverage sample $sample\" visibility=1 itemRgb=\"On\"\n";

#Read files into an array
my @data = <INPUT>;
my $lastidx = $#data;

my $chrCheck;
my $startCheck;
my $stopCheck;
my $tenxCheck;
my $twentyxCheck;
my $count = 0;
my $outputline;

my @array;
#Foreach line check coverage
foreach my $elem (@data){
    $count++;
    chomp $elem;
    my @line = split("\t", $elem);
    my $chr = $line[0];
    my $start = $line[1];
    my $stop = $line[2]+1;
    my $tenx = $line[10];
    my $twentyx = $line[11];
    my $color;
    #Check color
    if ($twentyx == 100){
        $color = "0,255,0"; #Green
    }elsif ($tenx == 100){
        $color = "255,180,0"; #Orange
    }else{ #color is red
        $color = "255,0,0"; #Red
    }
    #Push first line of file in array
    if ($count == 1){
        push(@array, $chr);
        push(@array, $start);
        push(@array, $stop);
        push(@array, $color);
        next;
    }
    #Check if chr in newline matches chr in array
    my $arrayChr = $array[0];
    my $arrayStart = $array[1];
    my $arrayStop = $array[2];
    my $arrayColor = $array[3];
    if ($arrayChr == $chr){
        #print "$color\t$arrayColor\n";
        $outputline = "chr$arrayChr\t$arrayStart\t$arrayStop\t$sample\t0\t+\t$arrayStart\t$arrayStop\t$arrayColor\n";
        #print "$outputline";
        #Continue checks, now for start en end positions
        if ($start == $arrayStop){ #region adjacent to previous regions, try to merge by comparing colors
            #compare colors
            if ($color eq $arrayColor){ #regions can be merged, update end of region
                splice(@array, 2, 1, $stop);
            }else { #regions don't have same color, write away info in array and push new line into array
                print OUTPUT $outputline;
                #empty array
                undef(@array);
                push(@array, $chr);
                push(@array, $start);
                push(@array, $stop);
                push(@array, $color);
            }
        }else {#region not adjacent, write away info in array and push new line into array
            print OUTPUT $outputline;
            #empty array
            undef(@array);
            push(@array, $chr);
            push(@array, $start);
            push(@array, $stop);
            push(@array, $color);
        }
    }else {
        #print "BLAAT\t$arrayChr\t$chr\n";
        #Chrs not equal, new region, print contents from array;
        print OUTPUT $outputline;
        undef(@array);
        push(@array, $chr);
        push(@array, $start);
        push(@array, $stop);
        push(@array, $color);
    }
    #If last line of file reached, write away all info from array
    if ($count == $lastidx+1){
        print OUTPUT "chr$chr\t$start\t$stop\t$sample\t0\t+\t$start\t$stop\t$color\n";
    }
}

close(INPUT);
close(OUTPUT);

#Substitute chr23 and chr24 with chrX and chrY respectively
`perl -pi -e 's/^chr23/chrX/g' $output`;
`perl -pi -e 's/^chr24/chrY/g' $output`;

#Cleanup temp files
`rm $output.tmp1.txt`;
`rm $output.tmp2.txt`;

print "\nFinished merging!\n";

sub usage {
        print <<EOF;
#########################################################################################
This script merges two coverage interval summary files for one sample and creates a *.bed
file with coverage profiles, which can be displayed in the UCSC Genome Browser.
#########################################################################################
Usage: ./create_coverage_bed.pl
\t-sample                   Sample name.
\t-file1                    First input coverage interval summary file.
\t-file2                    Second input coverage interval summary file.
\t-output                   Output bed file
#########################################################################################
EOF
 
}
