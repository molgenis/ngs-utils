#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use File::Glob ':glob';
use File::stat;
use Time::localtime;

#############################
##USAGE: perl check_output2.pl <working_directory> <output_file> <report_output_file>
#############################


##perl check_lane_output.pl /target/gpfs2/gcc/from_millipede/third_batch /data/gcc/tools/scripts/test.txt /data/gcc/tools/scripts/test.txt /data/gcc/tools/scripts/testreport.txt


my $dir;

my $input_dir = $ARGV[0];
my $output_dir = $ARGV[1];
my $reportoutput_dir = $ARGV[2];
chomp $input_dir;
chomp $output_dir;
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

##Retrieve files containing QC information
my @files = glob "$input_dir/*/*.runtime.log";
my @steplogs = glob "$input_dir/*/*.log";
my @filesizes = glob "$input_dir/*/*";

##Open output directory
open (OUTPUT, ">$output_dir");
open (REPORTOUTPUT, ">$reportoutput_dir");

##Print output header
$header = "sample\tal00.fastqc.ftl\tal01.bwa_align_pair1\tal02.bwa_align_pair2\tal03.bwa_sampe\tal04.sam_to_bam\tal05.sam_sort\tal06.picardQC\tal07.mark_duplicates\tal08.realign\tal09.fixmates\tal10.covariates_before\tal11.recalibrate\tal12.sam_sort\tal13.covariates_after\tal14.analyze_covariates\n";
print OUTPUT $header;

##Print header for report
my $reportheader = "sample\tlane_analysis_step (of 15)\tanalysis_finished\t\tpurity_filter_reads\tpercent_pf_reads_aligned\tpercent_pf_reads_aligned_in_pairs\tstrand_balance\t\tmedian_insert_size\tmean_insert_size\tstandard_deviation\t\tread_pair_duplicates\tread_pair_optical_duplicates\tpercent_duplication\n";
print REPORTOUTPUT $reportheader;

##Walk through log file paths
while ($_ = shift @files) {
    #print $_ . "\n";
    my %hash=();
    my $sample;
    my $file;
    my $date;
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
    open (FILE,"<$_") or die "Can't open $_$!";
    my $date_string = ctime(stat($_)->mtime); #Retrieve date of file
    if ($_ =~ m/$input_dir\/(.+)\/(.+L[0-9]{1,1}).+v37.(.+).runtime.log/gs){                                         #Retrieve log files
        $sample = $1;
        $file = $2;
        $date = $3;
        print OUTPUT $file;
        print REPORTOUTPUT $file;
    }
    while($logline = <FILE>){                                               #Retrieve all lines starting with Completed and push them into array.
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
        foreach $key (sort (keys(%hash))) {                                     #Loop through hash sorting by key
            print OUTPUT "\t$hash{$key}";                                       #Print all values per key
        }
        my $count_steps = scalar keys %hash;    #Count number of steps in hash
        print REPORTOUTPUT "\t$count_steps";
        
        if ($count_steps >= 15){    #If all 15 steps succeeded print timestamp of logfile
            print REPORTOUTPUT "\t$date_string\t";
        }else{  #Else print -
            $date_string = "-";
            print REPORTOUTPUT "\t$date_string\t";
        }
        #Print metrics from picardQC step
        if ($count_steps > 7){  #If step picardQC is passed assume all files are created
            my $open_asum = "$input_dir/$sample/$file.al06.picardQC.ftl.human_g1k_v37.$date.AlignmentSummaryMetrics";   #Open alignment summary metrics file
            #print "$open_asum\n";
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
                
                #Open insert size metrics file
                my $open_ismet = "$input_dir/$sample/$file.al06.picardQC.ftl.human_g1k_v37.$date.CollectInsertSizeMetrics";
                open (ISFILE,"<$open_ismet") or die "Can't open $open_ismet$!";
                while($isline = <ISFILE>){
                # if line matches data push into one array!
                    if ($isline =~ m/.+\t[0-9].+\t[0-9].+\t.+/gs){
                    chomp $isline;
                    @isinformation = split('\t', $isline);
                    }
                }
                $median_insert_size = $isinformation[0];
                $mean_insert_size = $isinformation[3];
                $sd = $isinformation[4];
                
                #Open duplicates metrics
                my $open_dupmet = "$input_dir/$sample/$file.al07.mark_duplicates.ftl.human_g1k_v37.$date.dedup.metrics";
                open (DUPFILE,"<$open_dupmet") or die "Can't open $open_dupmet$!";
                while($dupline = <DUPFILE>){
                # if line matches data push into one array!
                    if ($dupline =~ m/.+\t[0-9].+\t[0-9].+\t.+/gs){
                    chomp $dupline;
                    @dupinformation = split('\t', $dupline);
                    }
                }
                $read_pair_duplicates = $dupinformation[5];
                $read_pair_optical_duplicates = $dupinformation[6];
                $pct_duplication = $dupinformation[7]*100;
                
            }
        }else{
            $pf_reads = "-";
            $pct_pf_reads_aligned = "-";
            $pct_pf_reads_aligned_in_pairs = "-";
            $sb = "-";
            $median_insert_size = "-";
            $mean_insert_size = "-";
            $sd = "-";
            $read_pair_duplicates = "-";
            $read_pair_optical_duplicates = "-";
            $pct_duplication = "-";
        }
    
    print REPORTOUTPUT "\t$median_insert_size\t$mean_insert_size\t$sd\t";
    print REPORTOUTPUT "\t$read_pair_duplicates\t$read_pair_optical_duplicates\t$pct_duplication";
    
    print OUTPUT "\n";
    print REPORTOUTPUT "\n";
    undef %hash;
    undef @array;
}

print OUTPUT "\n\n\n";

#Print header corresponding to filesizes table
my $header_filesizes = "sample\tal01.bwa_align_pair1\tal02.bwa_align_pair2\tal03.bwa_sampe\tal04.sam_to_bam\tal05.sam_sort\tal07.mark_duplicates\tal08.realign\tal09.fixmates\tal11.recalibrate\tal12.sam_sort\tsam_sort_index\n";
print OUTPUT $header_filesizes;

#Retrieve all filesizes and print them to table
while ($_ = shift @filesizes) {
    #print $_ . "\n";
    my %sizehash=();
    my $filesize = "-";
    open (FILE,"<$_") or die "Can't open $_$!";
    if ($_ =~ m/$input_dir\/.+\/(.+).ftl.+.1.sai/gs || $_ =~ m/$input_dir\/.+\/(.+).ftl.+.2.sai/gs || $_ =~ m/$input_dir\/.+\/(.+).ftl.+.sam/gs || $_ =~ m/$input_dir\/.+\/(.+).ftl.+.bam/gs){  #If file with one of these extensions exist
        if ($_ =~ m/$input_dir\/.+\/(.+).ftl.+.1.sli/gs){    
        
        }elsif ($_ =~ m/$input_dir\/.+\/(.+).al01.+.ftl.+.1.sai/gs){                                         
            my $file = $1;
            print OUTPUT "\n$file";
            my $filename = $_;
            $filesize = -s $filename;
            #print "File: $file\t$filesize\n";
            print OUTPUT "\t$filesize";
        }elsif ($_ =~ m/$input_dir\/.+\/(.+).ftl.+.2.sai/gs){
            my $file = $1;
            my $filename = $_;
            $filesize = -s $filename;
            #print "File: $file\t$filesize\n";
            print OUTPUT "\t$filesize";
        }elsif ($_ =~ m/$input_dir\/.+\/(.+).ftl.+.sam/gs){
            my $file = $1;
            my $filename = $_;
            $filesize = -s $filename;
            #print "File: $file\t$filesize\n";
            print OUTPUT "\t$filesize";
        }elsif ($_ =~ m/$input_dir\/.+\/(.+).ftl.+.recal.sorted.bam/gs){
            my $file = $1;
            my $filename = $_;
            $filesize = -s $filename;
            #print "File: $file\t$filesize\n";
            print OUTPUT "\t$filesize";
        }elsif ($_ =~ m/$input_dir\/.+\/(.+).ftl.+.bam.bai/gs){
            #remove .bam.bai files
        }elsif ($_ =~ m/$input_dir\/.+\/(.+).ftl.+.bam/gs) {
            my $file = $1;
            my $filename = $_;
            $filesize = -s $filename;
            #print "File: $file\t$filesize\n";
            print OUTPUT "\t$filesize";
        }
    }
    #print OUTPUT "\n";
}