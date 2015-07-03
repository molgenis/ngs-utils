#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use File::Glob ':glob';

my ($help, $pindelvcf, $ugvcf, $outputvcf);

#### get options
GetOptions(
                "h"                         => \$help,
                "pindelVCF=s"               => \$pindelvcf,
                "unifiedGenotyperVCF=s"     => \$ugvcf,
                "outputVCF=s"               => \$outputvcf  
          );
usage() and exit(1) if $help;
# mandatory args
usage() and exit(1) unless $pindelvcf;
usage() and exit(1) unless $ugvcf;
usage() and exit(1) unless $outputvcf;

chomp $pindelvcf;
chomp $ugvcf;

#Open input and output files
open (PINDEL, "<", $pindelvcf ) or die $!;
open (UG, "<", $ugvcf ) or die $!;
open (OUTINDELVCF, ">", "$outputvcf" ) or die $!;

my %pindel;
my %ug;
my %ugannotation;
my $sample;
my %contigs;

print "Starting merge..\n\n";


#Read pindel VCF file
while (my $lines = <PINDEL>){ #Read by line
    chomp $lines;
    if ($lines =~ m/#CHROM.+/gs){ #Extract samplename
        my @array = split("\t",$lines);
        $sample = $array[9];
        #print "SAMPLE: $sample\n\n";
        #print OUTPUT "$sample";
    }
    if ($lines !~ m/^#.+/gs){ #Read SNPs plus info
        my @line = split("\t", $lines);
        my $chr = $line[0];

        my $pos = $line[1];
        my $ref = $line[3];
        my $alt = $line[4];
        my $filt = $line[6];
        my $info = $line[7];
        my $sampleinfo = $line[9];
        my @splitsampleinfo = split(":", $sampleinfo);
        my $genotype = $splitsampleinfo[0];

        my $key = $chr . "&&&" . $pos;
        my $value = $ref . "&&&" . $alt;# . "_" . $genotype;
        $pindel{ $key } = $value;
    }
}
close(PINDEL);

#Read unified genotyper VCF file
while (my $lines = <UG>){ #Read by line
    chomp $lines;
    if ($lines =~ m/^##.+/gs){
        print OUTINDELVCF "$lines\n";
    }
    if ($lines =~ m/#CHROM.+/gs){ #Extract samplename
        my @array = split("\t",$lines);
        $sample = $array[9];
    }
    if ($lines !~ m/^#.+/gs){ #Read SNPs plus info
        my @line = split("\t", $lines);
        my $chr = $line[0];

        my $pos = $line[1];
        my $id = $line[2];
        my $ref = $line[3];
        my $alt = $line[4];
        my $qual = $line[5];
        my $filt = $line[6];
        my $info = $line[7];
        my $format = $line[8];
        my $sampleinfo = $line[9];
        my @splitsampleinfo = split(":", $sampleinfo);
        my $genotype = $splitsampleinfo[0];

        my $key = $chr . "&&&" . $pos;
        my $value = $ref . "&&&" . $alt;# . "_" . $genotype;
        my $ugvalue = "$id\&&&$qual\&&&$filt\&&&$info\&&&$format\&&&$sampleinfo";
        $ug{ $key } = $value;
        $ugannotation{ $key } = $ugvalue;
    }
}
close(UG);

my %events;

#Compare keys
my @extra;
for ( sort keys %pindel ) { #check if chr_pos combination exists in both files
    unless ( exists $ug{$_} ) {
        push (@extra, $_);
        next;
    }

    if ($pindel{$_} eq $ug{$_}){
        my @first = split("&&&", $_);
        my @second = split("&&&", $ug{$_});
        my @third = split("&&&", $ugannotation{$_});
        
        print OUTINDELVCF $first[0] . "\t";
        print OUTINDELVCF $first[1] . "\t";
        print OUTINDELVCF $third[0] . "\t";
        print OUTINDELVCF $second[0] . "\t";
        print OUTINDELVCF $second[1] . "\t";
        print OUTINDELVCF $third[1] . "\t";
        print OUTINDELVCF $third[2] . "\t";
        print OUTINDELVCF $third[3] . "\t";
        print OUTINDELVCF $third[4] . "\t";
        print OUTINDELVCF $third[5];
        print OUTINDELVCF "\n";
        
        #Retrieve event length
        my $reflength = length($second[0]);
        my $altlength = length($second[1]);
        my $eventlength = ($altlength-$reflength);
        if (exists $events{$eventlength}) {
            $events{$eventlength}++;
        } else {
            $events{$eventlength} = 1;
        } 
    }

    if ($pindel{$_} ne $ug{$_}){ #Check if detected event differs between both outputs
        my @first = split("&&&", $_);
        my @second = split("&&&", $ug{$_});
        my @third = split("&&&", $ugannotation{$_});
        my @pindelsecond = split("&&&",$pindel{$_});
        
        print OUTINDELVCF $first[0] . "\t";
        print OUTINDELVCF $first[1] . "\t";
        print OUTINDELVCF $third[0] . "\t";
        print OUTINDELVCF $second[0] . "\t";
        print OUTINDELVCF $second[1] . "\t";
        print OUTINDELVCF $third[1] . "\t";
        print OUTINDELVCF $third[2] . "\t";
        print OUTINDELVCF $third[3] . ";PindelREF=" . $pindelsecond[0] . ";PindelALT=" . $pindelsecond[1] . "\t";
        print OUTINDELVCF $third[4] . "\t";
        print OUTINDELVCF $third[5];
        print OUTINDELVCF "\n";

        #Retrieve event length
        my $reflength = length($second[0]);
        my $altlength = length($second[1]);
        my $eventlength = ($altlength-$reflength);
        if (exists $events{$eventlength}) {
            $events{$eventlength}++;
        } else {
            $events{$eventlength} = 1;
        } 
    }
}

print "Merging finished succesfull!\n";

#Write event lengths
my @lengths;
my @counts;
foreach (sort {$a<=>$b} keys %events) {
    push(@lengths, $_);
    push(@counts, $events{$_});
}

########
#Uncomment this part to generate graph with counts of SV events
########

#use GD::Graph::bars;

#Create eventsize graph
#open (GRAPH, "> $outputvcf.events.png");

#my @data = ( 
#    [@lengths], 
#    [@counts],
#); 

#my $my_graph = GD::Graph::bars->new(1200, 700);
#$my_graph->set_x_label_font(GD::gdMediumBoldFont);
#$my_graph->set_x_axis_font(GD::gdMediumBoldFont);
#$my_graph->set_y_label_font(GD::gdMediumBoldFont);
#$my_graph->set_y_axis_font(GD::gdMediumBoldFont);

#$my_graph->set( 
#x_label => 'Event size', 
#y_label => 'Occurence', 
#title => 'Number of found events ordered by event length', 
#two_axes => 1, 
#y_min_value => 0, 
#y_tick_number => 5,
#long_ticks => 1, 
#x_ticks => 0, 
#legend_marker_width => 24, 
#line_width => 2,
#bar_spacing => 4, 
#transparent => 0, 
#) or warn $my_graph->error;
#my $image = $my_graph->plot(\@data) or die $my_graph->error; #print graph
#print GRAPH $image->png;

sub usage {
        print <<EOF;
#########################################################################################
This script merges the VCF files generated by Pindel with those generated by the Uni-
fied Genotyper (GATK). The script merges calls when both methods detected an event on the
same coordinates. When two events detected by both methods have different REF/ALT alleles
the call of UG will be used in the output and annotated with the Pindel call in the INFO
field. NOTE: The output VCF file is unsorted!
#########################################################################################
Usage: ./merge_SVs.pl
\t-pindelVCF                VCF file produced by Pindel to merge.
\t-unifiedGenotyperVCF      VCF file produced by GATK Unified Genotyper to merge. This
                                  VCF file only contains INDELS!
\t-outputVCF                Output unsorted VCF file.
#########################################################################################
EOF
 
}
