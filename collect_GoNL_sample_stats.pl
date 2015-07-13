#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use File::Glob ':glob';
use File::stat;
use Time::localtime;
use Spreadsheet::WriteExcel;


#############################
##USAGE: perl collect_GoNL_sample_stats.pl <working_directory> <report_output_file>
##Created 2011-09-26
#############################


##perl collect_GoNL_sample_stats.pl /target/gpfs2/gcc/from_millipede/fourth_batch /data/gcc/tools/scripts/testreport.txt


my $dir;

my $input_dir = $ARGV[0];
my $output_xls = $ARGV[1];
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
my $analyzed_samples = 0;
my $finished_lanes = 0;

##Retrieve files containing QC information
my @files = glob "$input_dir/*/*.runtime.log";

##Create and open output XLS file
my $workbook  = Spreadsheet::WriteExcel->new($output_xls);
my $worksheet = $workbook->add_worksheet();

##Set color variables
my $green = $workbook->add_format(fg_color => 'green');
my $orange = $workbook->add_format(fg_color => 'orange');
my $red = $workbook->add_format(fg_color => 'red');
my $blue = $workbook->add_format(fg_color => 'blue');

#Set italic font
my $italic = $workbook->add_format(italic => 1);

##Set and print header for report
my $heading = $workbook->add_format(bold => 1);
$worksheet->write(0, 0, "sample", $heading);
$worksheet->write(0, 1, "current_lane_analysis_step (of 26)", $heading);
$worksheet->write(0, 2, "analysis_finished", $heading);
$worksheet->write(0, 4, "purity_filter_reads", $heading);
$worksheet->write(0, 5, "percent_pf_reads_aligned", $heading);
$worksheet->write(0, 6, "percent_pf_reads_aligned_in_pairs", $heading);
$worksheet->write(0, 7, "strand_balance", $heading);
$worksheet->write(0, 9, "median_insert_size", $heading);
$worksheet->write(0, 10, "mean_insert_size", $heading);
$worksheet->write(0, 11, "standard_deviation", $heading);
$worksheet->write(0, 13, "number_of_SNPs", $heading);
$worksheet->write(0, 14, "percentage_dbSNP", $heading);
$worksheet->write(0, 15, "ti/tv_known", $heading);
$worksheet->write(0, 16, "ti/tv_novel", $heading);
$worksheet->write(0, 17, "percentage_all_comp_het_called_het", $heading);
$worksheet->write(0, 18, "percentage_known_comp_het_called_het", $heading);
$worksheet->write(0, 19, "percentage_nonref_sensitivity", $heading);
$worksheet->write(0, 20, "percentage_nonref_discrepancy", $heading);
$worksheet->write(0, 21, "percentage_overall_concordance", $heading);
$worksheet->write(0, 23, "contamination_check", $heading);

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
    my $num_snps;
    my $pct_dbsnp;
    my $titv_known;
    my $titv_novel;
    my $all_comp_het_call_het;
    my $known_comp_het_call_het;
    my $nonref_sens;
    my $nonref_disc;
    my $overall_conc;
    my @selfrg;
    my @bestrg;
    my $seq_sm;
    my $rg;
    my $selfibd;
    my $mix;
    my $ref_a2;
    my $exhet;
    my $hetseq;
    my $hetaf;
    my $hetdiff;
    my $best_rg;
    my $best_sm;
    my $bestibd;
    my $date_string = ctime(stat($_)->mtime); #Retrieve date of file
    if ($_ =~ m/$input_dir\/(.+)\/([AGR]{1,1}[0-9]{1,3}[ABCDabcd]{1,1}).human_g1k_v37.runtime.log/gs){                                         #Retrieve log files
        $sample = $1;
        $file = $2;
        #print REPORTOUTPUT $file;
        open (FILE,"<$_") or die "Can't open $_$!";
        $analyzed_samples++;
        $worksheet->write($analyzed_samples, 0, "$file", $italic);
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
        
            my $count_steps = scalar keys %hash;    #Count number of steps in hash
            $worksheet->write($analyzed_samples, 1, "$count_steps");

            if ($count_steps >= 4){    #If all steps succeeded print timestamp of logfile
                $worksheet->write($analyzed_samples, 2, "$date_string");
                $finished_lanes++;                          #Count finished lanes
            }else{  #Else print -
                $date_string = "-";
                $worksheet->write($analyzed_samples, 2, "$date_string");
            }

            #Print metrics from picardQC step
            #if ($count_steps > 2){  #If step picardQC is passed assume all QC files are created
            my $open_asum = "$input_dir/$sample/$file.vc02.picardQC.ftl.human_g1k_v37.AlignmentSummaryMetrics";   #Open alignment summary metrics file
            my $open_ismet = "$input_dir/$sample/$file.vc02.picardQC.ftl.human_g1k_v37.CollectInsertSizeMetrics";
            my $open_con_check = "$input_dir/$sample/$file.vc01.unified_genotyper.ftl.human_g1k_v37.qc_check_snps.concordance.txt";
            my $open_vbID_selfRG = "$input_dir/$sample/$file.vc03b.contamination_checker.human_g1k_v37.selfRG";
            my $open_vbID_bestRG = "$input_dir/$sample/$file.vc03b.contamination_checker.human_g1k_v37.bestRG";
                if (-e $open_asum){
                    open (ASUMFILE,"<$open_asum") or die "Can't open $open_asum$!";
                    while($asumline = <ASUMFILE>){
                        # if line matches data push into one array!
                        if ($asumline =~ m/^PAIR\t[0-9]+\t.+/gs){
                            chomp $asumline;
                            @asuminformation = split('\t', $asumline);
                            $pf_reads = $asuminformation[2];
                            $pct_pf_reads_aligned = $asuminformation[6]*100;
                            $pct_pf_reads_aligned_in_pairs = $asuminformation[14]*100;
                            $sb = $asuminformation[16]*100;
                            #print REPORTOUTPUT "\t$pf_reads\t$pct_pf_reads_aligned\t$pct_pf_reads_aligned_in_pairs\t$sb\t";
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
                            #print REPORTOUTPUT "\t$median_insert_size\t$mean_insert_size\t$sd\t";
                        }
                    }
                }else{
                    $median_insert_size = "-";
                    $mean_insert_size = "-";
                    $sd = "-";
                }
                
                my $conline;
                if (-e $open_con_check){
                    open (CONFILE,"<$open_con_check") or die "Can't open $open_con_check$!";
                    while ($conline = <CONFILE>){
                        chomp $conline;
                        if ($conline =~ m/.+,([0-9]{1,}),(.+),(.+),(.+),(.+),(.+),(.+),(.+),(.+)/gs){
                            $num_snps = $1;
                            $pct_dbsnp = $2;
                            $titv_known = $3;
                            $titv_novel = $4;
                            $all_comp_het_call_het = $5;
                            $known_comp_het_call_het = $6;
                            $nonref_sens = $7;
                            $nonref_disc = $8;
                            $overall_conc = $9;
                            #print REPORTOUTPUT "\t$num_snps\t$pct_dbsnp\t$titv_known\t$titv_novel\t$all_comp_het_call_het\t$known_comp_het_call_het\t$nonref_sens\t$nonref_disc\t$overall_conc\t";
                        }
                    }
                }else{
                    $num_snps = "-";
                    $pct_dbsnp = "-";
                    $titv_known = "-";
                    $titv_novel = "-";
                    $all_comp_het_call_het = "-";
                    $known_comp_het_call_het = "-";
                    $nonref_sens = "-";
                    $nonref_disc = "-";
                    $overall_conc = "-";
                }
        
                #print REPORTOUTPUT;
                $worksheet->write($analyzed_samples, 4, "$pf_reads");
                $worksheet->write($analyzed_samples, 5, "$pct_pf_reads_aligned");
                $worksheet->write($analyzed_samples, 6, "$pct_pf_reads_aligned_in_pairs");
                $worksheet->write($analyzed_samples, 7, "$sb");
                
                #Write insert size metrics
                $worksheet->write($analyzed_samples, 9, "$median_insert_size");
                $worksheet->write($analyzed_samples, 10, "$mean_insert_size");
                $worksheet->write($analyzed_samples, 11, "$sd");
                
                #Write concordance metrics
                $worksheet->write($analyzed_samples, 13, "$num_snps");
                $worksheet->write($analyzed_samples, 14, "$pct_dbsnp");
                $worksheet->write($analyzed_samples, 15, "$titv_known");
                $worksheet->write($analyzed_samples, 16, "$titv_novel");
                $worksheet->write($analyzed_samples, 17, "$all_comp_het_call_het");
                $worksheet->write($analyzed_samples, 18, "$known_comp_het_call_het");
                $worksheet->write($analyzed_samples, 19, "$nonref_sens");
                $worksheet->write($analyzed_samples, 20, "$nonref_disc");
                $worksheet->write($analyzed_samples, 21, "$overall_conc");
                
                undef %hash;
                undef @array;
                
                my $selfrgline;
                my $bestrgline;
                my $vbID_column = 23;
                if (-e $open_vbID_selfRG && -e $open_vbID_bestRG){
                    open (SELFRG, "<$open_vbID_selfRG") or die "Can't open $open_vbID_selfRG$!";
                    open (BESTRG, "<$open_vbID_bestRG") or die "Can't open $open_vbID_bestRG$!";
                    while ($selfrgline = <SELFRG>){
                        chomp $selfrgline;
                        if ($selfrgline !~ m/^SEQ_SM.+/gs){
                            #print "$selfrgline\n";
                            @selfrg = split('\t', $selfrgline);
                            $seq_sm = $selfrg[0];
                            $rg = $selfrg[1];
                            $selfibd = $selfrg[3];
                            $mix = $selfrg[24];
                            $ref_a2 = $selfrg[15];
                            $exhet = $selfrg[23];
                            $hetaf = $selfrg[21];
                            $hetseq = $selfrg[22];
                            my $self_output = "RG=$rg;SELF_IBD=$selfibd;%MIX=$mix;%REF_A2=$ref_a2;EXHET=$exhet;%HET_AF=$hetaf;%HET_SEQ=$hetseq";
                            #print "$self_output\n";
                            if ("$sample" ne "$seq_sm"){
                                die "Sample $seq_sm in selfRG file doesn't match with $sample! Comparing wrong files!";
                            }
                            #if ($selfibd < 0.90 ){
                            #    $worksheet->write($analyzed_samples, $vbID_column, "$self_output", $orange);
                            #    $vbID_column++;
                            #}elsif ($mix > 0.05){
                            #    
                            #}
                            if ($selfibd eq "N/A" || $mix eq "N/A" || $ref_a2 eq "N/A" || $exhet eq "N/A"){ #Check for values containing Not Available
                                $worksheet->write($analyzed_samples, $vbID_column, "$self_output", $blue);
                                $vbID_column++;
                            }elsif ($selfibd < 1 && $mix > 0 && $ref_a2 > 0.01){
                                #Possible contamination
                                $worksheet->write($analyzed_samples, $vbID_column, "$self_output", $red);
                                $vbID_column++;
                            }elsif ($selfibd < 0.95 && $mix > 0.05){
                                #Possible sample swap -> check bestRG and bestIBD
                                while ($bestrgline = <BESTRG>){
                                    @bestrg = split('\t', $bestrgline);
                                    $best_rg = $bestrg[1];
                                    $best_sm = $bestrg[2];
                                    $bestibd = $bestrg[3];
                                    if ("$best_rg" eq "$rg"){
                                        $self_output = "RG=$rg;SELF_IBD=$selfibd;%MIX=$mix;%REF_A2=$ref_a2;EXHET=$exhet;%HET_AF=$hetaf;%HET_SEQ=$hetseq;BEST_SM=$best_sm;BEST_IBD=$bestibd";
                                        $worksheet->write($analyzed_samples, $vbID_column, "$self_output", $red);
                                        $vbID_column++;
                                    }
                                }
                            }elsif ($exhet < 0.90 || $exhet > 1.10){
                                #Excessive heterozygosity 0.05 deviation from 1
                                $worksheet->write($analyzed_samples, $vbID_column, "$self_output", $orange);
                                $vbID_column++;
                            }elsif ($selfibd < 0.95 || $mix > 0.05 || $ref_a2 > 0.01){
                                #One of the filtered values below/above thresholds from combined filter
                                $worksheet->write($analyzed_samples, $vbID_column, "$self_output", $orange);
                                $vbID_column++;
                            }else{
                                #Print normal values
                                $worksheet->write($analyzed_samples, $vbID_column, "$self_output");
                                $vbID_column++;
                            }
                        }
                    }
                }
                
            #}
        
        #}
    }
}

my $analyzed_row = $analyzed_samples + 2;
my $finished_row = $analyzed_samples + 3;

$worksheet->write($analyzed_row, 0, "#Samples analyzed:");
$worksheet->write($finished_row, 0, "#Samples finished merging:");
$worksheet->write($analyzed_row, 1, "$analyzed_samples");
$worksheet->write($finished_row, 1, "$finished_lanes");

#print REPORTOUTPUT "\nSamples analyzed: $analyzed_samples\n";
#print REPORTOUTPUT "Samples finished merging: $finished_lanes\n";