#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use File::Glob ':glob';
use File::stat;
use Time::localtime;

#############################
##USAGE: perl collect_GoNL_lane_statsV2.pl <working_directory> <report_output_file>
##Created 2011-09-05
#############################


##perl collect_GoNL_lane_stats.pl /target/gpfs2/gcc/from_millipede/fourth_batch /data/gcc/tools/scripts/testreport.txt


my $dir;

my $input_dir = $ARGV[0];
my $reportoutput_dir = $ARGV[1];
chomp $input_dir;
my $logline;
my @logs;
my @filenames;
my $count=0;
my @array;
my $header;
my $folder;
my %loghash;
my $stepname;
my $node;
my $asumline;
my $isline;
my $dupline;
my $analyzed_lanes = 0;
my $finished_lanes = 0;

##Retrieve files containing QC information
my @files = glob "$input_dir/*/*.runtime.log";
my @steplogs = glob "$input_dir/*/*.log";

##Open output directory
open (REPORTOUTPUT, ">$reportoutput_dir");

##Print header for report
my $reportheader = "sample\tlane_analysis_step\tanalysis_finished\t\tpurity_filter_reads\tpercent_pf_reads_aligned\tpercent_pf_reads_aligned_in_pairs\tstrand_balance\t\tmedian_insert_size\tmean_insert_size\tstandard_deviation\t\tread_pair_duplicates\tread_pair_optical_duplicates\tpercent_duplication\n";
print REPORTOUTPUT $reportheader;

##Walk through log file paths
while ($_ = shift @files) {
    #print $_ . "\n";
    my %hash=();
    my $sample;
    my $file;
    my @asuminformation;
    my @isinformation;
    my @dupinformation;
    my $pf_reads;
    my $pct_pf_reads_aligned;
    my $pct_pf_reads_aligned_in_pairs;
    my $sb;
    my $median_insert_size;
    my $mean_insert_size;
    my $sd;
    my $read_pair_duplicates;
    my $read_pair_optical_duplicates;
    my $pct_duplication;
    my $date_string = ctime(stat($_)->mtime); #Retrieve date of file
    if ($_ =~ m/$input_dir\/(.+)\/(.+L[0-9]{1,1}).+v37.runtime.log/gs){                                         #Retrieve log files
        $sample = $1;
        $file = $2;
        print REPORTOUTPUT $file;
        open (FILE,"<$_") or die "Can't open $_$!";
        $analyzed_lanes++;

        while($logline = <FILE>){                                               #Retrieve all lines starting with "Completed" and push them into array.
            chomp $logline;                                                     #Bottomline is pushed in array first
            if($logline =~ m/Completed (.+).ftl.+at.+in ([0-9]{1,}) seconds/gs){
                unshift (@array, $logline);
            }
        }
        close (FILE);
        foreach my $element(@array){
            if($element =~ m/Completed (.+).ftl.+at.+in ([0-9]{1,}) seconds/gs){   #Retrieve stepname and elapsed time from lines
                my $step = $1;
                my $elapsedtime = $2;
                if (exists($hash{$step})){                                      #If step not exists push into hash
                    # if the key is found in the hash come here
                }else{
                # come here if the key is not found in the hash
                    $hash{ $step } = $elapsedtime;
                    }
                }
            }
            my $key;
            my $num_pipeline_steps = 0;
            my $qc_steps = 0;
            foreach $key (sort (keys(%hash))) {                                     #Loop through hash sorting by key
                if ($key =~ m/.+.GoNL[0-9]{2,2}/gs){                                #If key matches "GoNL" the 5 step template system is used..QC then in step 3
                    $num_pipeline_steps = 5;
                    $qc_steps = 3;
                }else{                                                              #Else the 15 step template system is used
                    $num_pipeline_steps = 15;
                    $qc_steps = 7;
                }
            }
        
            my $count_steps = scalar keys %hash;    #Count number of steps in hash
            print REPORTOUTPUT "\t$count_steps";

            if ($count_steps == $num_pipeline_steps){    #If all steps succeeded print timestamp of logfile
                print REPORTOUTPUT "\t$date_string\t";
                $finished_lanes++;                          #Count finished lanes
            }else{  #Else print -
                $date_string = "-";
                print REPORTOUTPUT "\t$date_string\t";
            }

            #Print metrics from picardQC step
            if ($count_steps > $qc_steps){  #If step picardQC is passed assume all QC files are created
            my $open_asum = "$input_dir/$sample/$file.al06.picardQC.ftl.human_g1k_v37.AlignmentSummaryMetrics";   #Open alignment summary metrics file
            my $open_ismet = "$input_dir/$sample/$file.al06.picardQC.ftl.human_g1k_v37.CollectInsertSizeMetrics";
            my $open_dupmet = "$input_dir/$sample/$file.al07.mark_duplicates.ftl.human_g1k_v37.dedup.metrics";
                if (-e $open_asum){
                    open (ASUMFILE,"<$open_asum") or die "Can't open $open_asum$!";
                    while($asumline = <ASUMFILE>){
                        # if line matches data push into one array!
                        if ($asumline =~ m/^PAIR\t[0-9]+\t.+/gs){
                            chomp $asumline;
                            @asuminformation = split('\t', $asumline);
                            my $pf_reads = $asuminformation[2];
                            my $pct_pf_reads_aligned = $asuminformation[6]*100;
                            my $pct_pf_reads_aligned_in_pairs = $asuminformation[14]*100;
                            my $sb = $asuminformation[16]*100;
                            print REPORTOUTPUT "\t$pf_reads\t$pct_pf_reads_aligned\t$pct_pf_reads_aligned_in_pairs\t$sb\t";
                        }
                    }
                }else{
                    $pf_reads = "-";
                    $pct_pf_reads_aligned = "-";
                    $pct_pf_reads_aligned_in_pairs = "-";
                    $sb = "-";
                }
                    #Open insert size metrics file
                if (-e $open_ismet){
                    open (ISFILE,"<$open_ismet") or die "Can't open $open_ismet$!";
                    while($isline = <ISFILE>){
                    # if line matches data push into one array!
                        if ($isline =~ m/.+\t[0-9].+\t[0-9].+\t.+/gs){
                            chomp $isline;
                            @isinformation = split('\t', $isline);
                            $median_insert_size = $isinformation[0];
                            $mean_insert_size = $isinformation[3];
                            $sd = $isinformation[4];
                            print REPORTOUTPUT "\t$median_insert_size\t$mean_insert_size\t$sd\t";
                        }
                    }
                }else{
                    $median_insert_size = "-";
                    $mean_insert_size = "-";
                    $sd = "-";
                }
                    #Open duplicates metrics
                if (-e $open_dupmet){
                    open (DUPFILE,"<$open_dupmet") or die "Can't open $open_dupmet$!";
                    while($dupline = <DUPFILE>){
                    # if line matches data push into one array!
                        if ($dupline =~ m/.+\t[0-9].+\t[0-9].+\t.+/gs){
                            chomp $dupline;
                            @dupinformation = split('\t', $dupline);
                            $read_pair_duplicates = $dupinformation[5];
                            $read_pair_optical_duplicates = $dupinformation[6];
                            $pct_duplication = $dupinformation[7]*100;
                            print REPORTOUTPUT "\t$read_pair_duplicates\t$read_pair_optical_duplicates\t$pct_duplication";
                        }
                    }
                }else{
                    $read_pair_duplicates = "-";
                    $read_pair_optical_duplicates = "-";
                    $pct_duplication = "-";
                }
            }
        print REPORTOUTPUT "\n";
        undef %hash;
        undef @array;
        #}
    }
}

print REPORTOUTPUT "\nLanes analyzed: $analyzed_lanes\n";
print REPORTOUTPUT "Lanes finished: $finished_lanes\n";