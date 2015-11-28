#!/usr/bin/perl -w
#use strict;
use warnings;
use diagnostics;
use File::Glob ':glob';

#############################
##USAGE: perl collect_stats_first_batch.pl <working_directory> <output_directory>
#############################

my $dir;
my @asummet;
my @dedupmet;
my @asummet_filtered;
my @dedupmet_filtered;
my @files;
my @hsmet;
my @hsmet_filtered;
my @ismet_filtered;
my $filename;
my $count=0;
my $asumline;
my $dedupline;
my $hsline;
my $isline;
my $input;
my @asuminformation;
my @dedupinformation;
my @hsinformation;
my @isinformation;
my $asheader;
my $ddheader;
my $isheader;

my $input_dir = $ARGV[0];
my $output_dir = $ARGV[1];
chomp $input_dir;
chomp $output_dir;

##Retrieve files containing QC information
@asummet = glob "$input_dir/outputs/*/*.AlignmentSummaryMetrics";
@dedupmet = glob "$input_dir/outputs/*/*dedup.metrics";
@hsmet = glob "$input_dir/outputs/*/*.HsMetrics";
@ismet = glob "$input_dir/outputs/*/*.CollectInsertSizeMetrics";

open (OUTPUT, ">$output_dir");

##Output headers
#my $asum_header="Sample\tTotal_reads\tAligned_reads\tQ>20_Aligned_reads\tQ>20_Median_mismatches\tBad_cycles\tStrand_balance\tPercentage_chimeras\tPercentage_adapters\n";
#my $dedup_header="Sample\tReads_without_mapped_mate_pair\tRead_pairs_mapped\tUnmapped_reads\tUnpaired_read_duplicates\tRead_pair_duplicates\tRead_pair_optical_duplicates\tPercent_duplication\n";
#my $hs_header="Sample\tGenome_size(bases)\tExome_target(bases)\tTotal_reads\tUnique_reads\tUnique_reads_aligned\tOn_target_bases\tMean_coverage_exome_targets\tPercentage_exome_targets_coverage<2\tFold_80_base_penalty\tPercentage_target_bases_>2X\tPercentage_target_bases_>10X\tPercentage_target_bases_>20X\tPercentage_target_bases_>30X\tHybrid_selection_penalty_10X\tHybrid_selection_penalty_20X\tHybrid_selection_penalty_30X\n";
#y $is_header="Sample\tMedian_insert_size\tMinimum_insert_size\tMaximum_insert_size\tMean_insert_size\tStandard_deviation\tNumber_of_read_pairs\tRead_pair_orientation\n";

foreach my $element (@asummet){
    #print $element . "\n";
    if ($element =~ m/$input_dir\/outputs\/.*\/(.+).AlignmentSummaryMetrics/gs){
        push (@asummet_filtered, $1);
        #print $1 . "\t";
        #print "hit!!\n";
    }
}

foreach my $element (@dedupmet){
    #print $element . "\n";
    if ($element =~ m/$input_dir\/outputs\/.*\/(.+).dedup.metrics/gs){
        push (@dedupmet_filtered, $1);
        #print $1 . "\t";
        #print "hit!!\n";
    }
}

foreach my $element (@ismet){
    #print $element . "\n";
    if ($element =~ m/$input_dir\/outputs\/.*\/(.+).CollectInsertSizeMetrics/gs){
        push (@ismet_filtered, $1);
        #print $1 . "\t";
        #print "hit!!\n";
    }
}

while ($_ = shift @asummet) {
    #print $_ . "\n";
    open (ASUM,"<$_") or die "Can't open $_$!";
    if ($_ =~ m/.+\/outputs\/.*\/[0-9]{1,6}_I.+_(FC.+)_HU.+.AlignmentSummaryMetrics/gs){
        $asummet_filtered = $1;
        #100803_I147_FC2023JABXX_L3_HUModqRACDIAAPE.human_g1k_v37
        #print $1 . "\t";
        #print "hit!!\n";
    }elsif ($_ =~ m/.+\/outputs\/.*\/(.+)\..+.AlignmentSummaryMetrics/gs){
        $asummet_filtered = $1;
    }
    while($asumline = <ASUM>){
        #print $line;
        if ($asumline =~ m/(CATEGORY\t.+)/gs){
            $asheader = "LANE\t$1";
        }
        # if line matches data push into one array!
        if ($asumline =~ m/^PAIR\t[0-9]+\t.+/gs){
            chomp $asumline;
            #print "SAMPLE\t" . $line . "\t";
            #print "hit asumline \n";
            #$input = $asumline . "\t";
            $input = $asummet_filtered . "\t" . $asumline;
            push (@asuminformation, $input);
            #print $input . "\n";
            $count++;
        }
    }
}

while ($_ = shift @dedupmet) {
    #print $_ . "\n";
    open (DEDUP,"<$_") or die "Can't open $_$!";
    if ($_ =~ m/.+\/outputs\/.*\/[0-9]{1,6}_I.+_(FC.+)_HU.+.dedup.metrics/gs){
        $dedupmet_filtered = $1;
        #print $1 . "\t";
        #print "hit!!\n";
    }
    while($dedupline = <DEDUP>){
        #print $line;
        if ($dedupline =~ m/(LIBRARY\t.+)/gs){
            $ddheader = "LANE\t$1";
        }
        # if line matches data push into one array!
        if ($dedupline =~ m/.+\t[0-9].+\t[0-9].+\t.+/gs){
            chomp $dedupline;
            #print "SAMPLE\t" . $line . "\t";
            #print "hit dedupline \n";
            #$input = $dedupline . "\t";
            $input = $dedupmet_filtered . "\t" . $dedupline;
            push (@dedupinformation, $input);
            #print $input . "\n";
            $count++;
        }
    }
}

while ($_ = shift @ismet) {
    #print $_ . "\n";
    open (IS,"<$_") or die "Can't open $_$!";
    if ($_ =~ m/.+\/outputs\/.*\/[0-9]{1,6}_I.+_(FC.+)_HU.+.CollectInsertSizeMetrics/gs){
        $ismet_filtered = $1;
        #print $1 . "\t";
        #print "hit!!\n";
    }elsif ($_ =~ m/.+\/outputs\/.*\/(.+)\..+.CollectInsertSizeMetrics/gs){
        $ismet_filtered = $1;
    }
    while($isline = <IS>){
        #print $line;
        if ($isline =~ m/(MEDIAN_INSERT_SIZE\t.+)/gs){
            $isheader = "LANE\t$1";
        }
        # if line matches data push into one array!
        if ($isline =~ m/.+\t[0-9].+\t[0-9].+\t.+/gs){
            chomp $isline;
            #print "SAMPLE\t" . $line . "\t";
            #print "hit hsline \n";
            $input = $ismet_filtered . "\t" . $isline;
            push (@isinformation, $input);
            #print $input . "\n";
            $count++;
        }
    }
}

print OUTPUT "###Alignment Information###\n";
print OUTPUT $asheader;

foreach my $asumelement (@asuminformation){
    my @asumarray = split('\t', $asumelement);
    print OUTPUT $asumelement . "\n";
}

print OUTPUT "\n\n###Duplicates Information###\n";
print OUTPUT $ddheader;

foreach my $dedupelement (@dedupinformation){
    print OUTPUT $dedupelement . "\n";
}

print OUTPUT "\n\n###Insert Size Information\n";
print OUTPUT $isheader;

foreach my $iselement (@isinformation){
    print OUTPUT $iselement . "\n";
}