#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use File::Glob ':glob';
use File::stat;
use Time::localtime;
use Spreadsheet::WriteExcel;


#############################
##USAGE: perl collect_in_house_samplestatsV2.pl <working_directory> <project> <report_output_file.xls>
#############################



my $dir;

my $input_dir = $ARGV[0];
my $project = $ARGV[1];
my $output_xls = $ARGV[2];
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
my $hsline;
my $dupline;
my $num_files=1;

##Retrieve files containing QC information
my @files = glob "$input_dir/*/*.runtime.log";
my @steplogs = glob "$input_dir/*/*.log";

##Create and open output XLS file
my $workbook  = Spreadsheet::WriteExcel->new($output_xls);
my $worksheet = $workbook->add_worksheet();

##Set color variables
my $green = $workbook->add_format(fg_color => 'green');
my $orange = $workbook->add_format(fg_color => 'orange');
my $red = $workbook->add_format(fg_color => 'red');

#Set italic font
my $italic = $workbook->add_format(italic => 1);

#Initialize coverage graph
my $coverage_chart = $workbook->add_chart( name => 'Coverage chart', type => 'column', embedded => 1);

##Set and print header for report
my $heading = $workbook->add_format(bold => 1);
$worksheet->write(0, 0, "$project", $heading);
$worksheet->write(0, 1, "current_lane_analysis_step (of 26)", $heading);
$worksheet->write(0, 2, "analysis_finished", $heading);
$worksheet->write(0, 4, "purity_filter_reads", $heading);
$worksheet->write(0, 5, "percent_pf_reads_aligned", $heading);
$worksheet->write(0, 6, "percent_pf_reads_aligned_in_pairs", $heading);
$worksheet->write(0, 7, "strand_balance", $heading);
$worksheet->write(0, 9, "median_insert_size", $heading);
$worksheet->write(0, 10, "mean_insert_size", $heading);
$worksheet->write(0, 11, "standard_deviation", $heading);
$worksheet->write(0, 13, "read_pair_duplicates", $heading);
$worksheet->write(0, 14, "read_pair_optical_duplicates", $heading);
$worksheet->write(0, 15, "percent_duplication", $heading);
$worksheet->write(0, 17, "enrichment_kit", $heading);
$worksheet->write(0, 18, "genome_size", $heading);
$worksheet->write(0, 19, "bait_territory", $heading);
$worksheet->write(0, 20, "target_territory", $heading);
$worksheet->write(0, 21, "number_of_pf_reads", $heading);
$worksheet->write(0, 22, "number_of_pf_unique_reads", $heading);
$worksheet->write(0, 23, "percentage_pf_unique_reads_aligned", $heading);
$worksheet->write(0, 24, "on_bait_bases", $heading);
$worksheet->write(0, 25, "near_bait_bases", $heading);
$worksheet->write(0, 26, "off_bait_bases", $heading);
$worksheet->write(0, 27, "on_target_bases", $heading);
$worksheet->write(0, 28, "percentage_off_bait_bases", $heading);
$worksheet->write(0, 29, "mean_bait_coverage", $heading);
$worksheet->write(0, 30, "mean_target_coverage", $heading);
$worksheet->write(0, 31, "percentage_usable_bases_on_target", $heading);
$worksheet->write(0, 32, "percentage_zero_cov_targets", $heading);
$worksheet->write(0, 33, "percentage_target_bases_2x", $heading);
$worksheet->write(0, 34, "percentage_target_bases_10x", $heading);
$worksheet->write(0, 35, "percentage_target_bases_20x", $heading);
$worksheet->write(0, 36, "percentage_target_bases_30x", $heading);
$worksheet->write(0, 38, "number_of_Q20_DP10_SNPs", $heading);
$worksheet->write(0, 39, "percentage_dbSNP", $heading);
$worksheet->write(0, 40, "ti/tv_known", $heading);
$worksheet->write(0, 41, "ti/tv_novel", $heading);
$worksheet->write(0, 42, "percentage_all_comp_het_called_het", $heading);
$worksheet->write(0, 43, "percentage_known_comp_het_called_het", $heading);
$worksheet->write(0, 44, "percentage_nonref_sensitivity", $heading);
$worksheet->write(0, 45, "percentage_nonref_discrepancy", $heading);
$worksheet->write(0, 46, "percentage_overall_concordance", $heading);

##Walk through log file paths
while ($_ = shift @files) {
    #print $_ . "\n";
    my %hash=();
    my $sample;
    my $file;
    my @asuminformation;
    my @isinformation;
    my @hsinformation;
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
    my $enrichment_kit;
    my $genome_size;
    my $bait_ter;
    my $tar_ter;
    my $num_pf_reads;
    my $num_pf_unique_reads;
    my $pct_pf_unique_reads_aligned;
    my $on_bait_bases;
    my $near_bait_bases;
    my $off_bait_bases;
    my $on_target_bases;
    my $pct_off_bait;
    my $mean_bait_coverage;
    my $mean_target_coverage;
    my $pct_usable_bases_on_target;
    my $pct_zero_cov_targets;
    my $pct_target_bases_2x;
    my $pct_target_bases_10x;
    my $pct_target_bases_20x;
    my $pct_target_bases_30x;
    my $conline;
    my $num_snps;
    my $pct_dbsnp;
    my $titv_known;
    my $titv_novel;
    my $all_comp_het_call_het;
    my $known_comp_het_call_het;
    my $nonref_sens;
    my $nonref_disc;
    my $overall_conc;
    my $dedup_devide = 2;
    open (FILE,"<$_") or die "Can't open $_$!";
    my $date_string = ctime(stat($_)->mtime); #Retrieve date of file
    if ($_ =~ m/$input_dir\/(.+)\/(.+L[0-9]{1,1}).+v37.runtime.log/gs){                                         #Retrieve log files
        $sample = $1;
        $file = $2;
        $worksheet->write($num_files, 0, "$file", $italic);
    }elsif($_ =~ m/$input_dir\/(.+)\/(.+).+v37.runtime.log/gs){
        $sample = $1;
        my @fastqclogs = glob "$input_dir/$sample/*.fastqcsummary.log";
        while ($_ = shift @fastqclogs) {
            #print "$_\n";
            if ($_ =~ m/$input_dir\/(.+)\/(.+).HSpe00.+v37.1.fastqcsummary.log/gs){
                $file = $2;
                $worksheet->write($num_files, 0, "$file", $italic);
            }
        }
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
            #print OUTPUT "\t$hash{$key}";                                       #Print all values per key
        }
        my $count_steps = scalar keys %hash;    #Count number of steps in hash
        $worksheet->write($num_files, 1, "$count_steps");
        
        if ($count_steps >= 26){    #If all 26 steps succeeded print timestamp of logfile
            $worksheet->write($num_files, 2, "$date_string");
        }else{  #Else print -
            $date_string = "-";
            $worksheet->write($num_files, 2, "$date_string");
        }
        #Print metrics from picardQC step
        my $open_asum = "$input_dir/$sample/$file.HSpe16.recalibrate.ftl.human_g1k_v37.recal.AlignmentSummaryMetrics";   #Open alignment summary metrics file
        #print "$open_asum\n";
        if (-e $open_asum){
            open (ASUMFILE,"<$open_asum") or die "Can't open $open_asum$!";
            while($asumline = <ASUMFILE>){
                # if line matches data push into one array!
                if ($asumline =~ m/^PAIR\t[0-9]+\t.+/gs){ #Paired-end data
                    chomp $asumline;
                    @asuminformation = split('\t', $asumline);
                    $pf_reads = $asuminformation[2];
                    $pct_pf_reads_aligned = $asuminformation[6]*100;
                    $pct_pf_reads_aligned_in_pairs = $asuminformation[14]*100;
                    $sb = $asuminformation[16]*100;
                    $dedup_devide = 2;
                }elsif ($asumline =~ m/^UNPAIRED\t[0-9]+\t.+/gs){ #Single-end data
                    chomp $asumline;
                    @asuminformation = split('\t', $asumline);
                    $pf_reads = $asuminformation[2];
                    $pct_pf_reads_aligned = $asuminformation[6]*100;
                    $pct_pf_reads_aligned_in_pairs = $asuminformation[14]*100;
                    $sb = $asuminformation[16]*100;
                    $dedup_devide = 1;
                }
            }
        }else{
            $pf_reads = "-";
            $pct_pf_reads_aligned = "-";
            $pct_pf_reads_aligned_in_pairs = "-";
            $sb = "-";
        }
        close (ASUMFILE);
         
        #Open insert size metrics file
        my $open_ismet = "$input_dir/$sample/$file.HSpe16.recalibrate.ftl.human_g1k_v37.recal.CollectInsertSizeMetrics";
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
                }
            }
        }else{
            $median_insert_size = "-";
            $mean_insert_size = "-";
            $sd = "-";
        }
        close (ISFILE);
                
        #Open HS metrics file
        my $open_hsmet = "$input_dir/$sample/$file.HSpe16.recalibrate.ftl.human_g1k_v37.recal.HsMetrics";
        if (-e $open_hsmet){
            open (HSFILE,"<$open_hsmet") or die "Can't open $open_hsmet$!";
            while($hsline = <HSFILE>){
            # if line matches data push into one array!
                if ($hsline =~ m/.+\t[0-9].+\t[0-9].+\t[0-9].+\t.+/gs){
                    chomp $hsline;
                    #print "$hsline\n";
                    @hsinformation = split('\t', $hsline);
                    $enrichment_kit = $hsinformation[0];
                    $genome_size = $hsinformation[1];
                    $bait_ter = $hsinformation[2];
                    $tar_ter = $hsinformation[3];
                    $num_pf_reads = $hsinformation[6];
                    $num_pf_unique_reads = $hsinformation[7];
                    $pct_pf_unique_reads_aligned = $hsinformation[11] * 100;
                    $on_bait_bases = $hsinformation[13];
                    $near_bait_bases = $hsinformation[14];
                    $off_bait_bases = $hsinformation[15];
                    $on_target_bases = $hsinformation[16];
                    $pct_off_bait = $hsinformation[18] * 100;
                    $mean_bait_coverage = $hsinformation[20];
                    $mean_target_coverage = $hsinformation[21];
                    $pct_usable_bases_on_target = $hsinformation[23] * 100;
                    $pct_zero_cov_targets = $hsinformation[25] * 100;
                    $pct_target_bases_2x = $hsinformation[27] * 100;
                    $pct_target_bases_10x = $hsinformation[28] * 100;
                    $pct_target_bases_20x = $hsinformation[29] * 100;
                    $pct_target_bases_30x = $hsinformation[30] * 100;
                }
            }
        }else {
            $enrichment_kit = "-";
            $genome_size = "-";
            $bait_ter = "-";
            $tar_ter = "-";
            $num_pf_reads = "-";
            $num_pf_unique_reads = "-";
            $pct_pf_unique_reads_aligned = "-";
            $on_bait_bases = "-";
            $near_bait_bases = "-";
            $off_bait_bases = "-";
            $on_target_bases = "-";
            $pct_off_bait = "-";
            $mean_bait_coverage = "-";
            $mean_target_coverage = "-";
            $pct_usable_bases_on_target = "-";
            $pct_zero_cov_targets = "-";
            $pct_target_bases_2x = "-";
            $pct_target_bases_10x = "-";
            $pct_target_bases_20x = "-";
            $pct_target_bases_30x = "-";   
        }
        close (HSFILE);
                
        #Open duplicates metrics
        my $open_dupmet = "$input_dir/$sample/$file.HSpe07.mark_duplicates.ftl.human_g1k_v37.dedup.metrics";
        if (-e $open_dupmet){
            open (DUPFILE,"<$open_dupmet") or die "Can't open $open_dupmet$!";
            while($dupline = <DUPFILE>){
            # if line matches data push into one array!
                if ($dupline =~ m/.+\t[0-9].+\t[0-9].+\t.+/gs){
                    chomp $dupline;
                    @dupinformation = split('\t', $dupline);
                    $read_pair_duplicates = $dupinformation[5];
                    $read_pair_optical_duplicates = $dupinformation[6];
                    $pct_duplication = ($dupinformation[7]/$dedup_devide)*100;
                }
            }
        }else{
            $read_pair_duplicates = "-";
            $read_pair_optical_duplicates = "-";
            $pct_duplication = "-";
        }
        close (DUPFILE);
                
        #Open concordance metrics
        my $open_conmet = "$input_dir/$sample/$file.HSpe26.concordance_check.calls_vs_array_concordance.txt";
        if (-e $open_conmet){
            open (CONFILE,"<$open_conmet") or die "Can't open $open_conmet$!";
            while($conline = <CONFILE>){
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
        close (CONFILE);
    
    #Write alignment stats
    $worksheet->write($num_files, 4, "$pf_reads");
    if ($pct_pf_reads_aligned ne "-"){
        if ($pct_pf_reads_aligned >= 90){
            $worksheet->write($num_files, 5, "$pct_pf_reads_aligned");
        }elsif ($pct_pf_reads_aligned >= 80){
            $worksheet->write($num_files, 5, "$pct_pf_reads_aligned", $orange);
        }else{
            $worksheet->write($num_files, 5, "$pct_pf_reads_aligned", $red);
        }
    }else{
        $worksheet->write($num_files, 5, "$pct_pf_reads_aligned");
    }
    if ($pct_pf_reads_aligned_in_pairs ne "-"){
        if ($pct_pf_reads_aligned_in_pairs >= 95){
            $worksheet->write($num_files, 6, "$pct_pf_reads_aligned_in_pairs");
        }elsif ($pct_pf_reads_aligned_in_pairs >= 80){
            $worksheet->write($num_files, 6, "$pct_pf_reads_aligned_in_pairs", $orange);
        }else{
            $worksheet->write($num_files, 6, "$pct_pf_reads_aligned_in_pairs", $red);
        }
    }else{
        $worksheet->write($num_files, 6, "$pct_pf_reads_aligned_in_pairs");
    }
    if ($sb ne "-"){
        if ($sb >= 49.5 && $sb <= 50.5){
            $worksheet->write($num_files, 7, "$sb");
        }elsif ($sb >= 48.5 && $sb <= 51.5){
            $worksheet->write($num_files, 7, "$sb", $orange);
        }else{
            $worksheet->write($num_files, 7, "$sb", $red);
        }
    }else{
        $worksheet->write($num_files, 7, "$sb");
    }
    
    #Write insert size metrics
    $worksheet->write($num_files, 9, "$median_insert_size");
    $worksheet->write($num_files, 10, "$mean_insert_size");
    $worksheet->write($num_files, 11, "$sd");
    
    #Write duplicate stats
    $worksheet->write($num_files, 13, "$read_pair_duplicates");
    $worksheet->write($num_files, 14, "$read_pair_optical_duplicates");
    if ($pct_duplication ne "-"){
        if ($pct_duplication <= 8){
            $worksheet->write($num_files, 15, "$pct_duplication");
        }elsif ($pct_duplication <=15){
            $worksheet->write($num_files, 15, "$pct_duplication", $orange);
        }else{
            $worksheet->write($num_files, 15, "$pct_duplication", $red);
        }
    }else{
        $worksheet->write($num_files, 15, "$pct_duplication");
    }
    
    #Write hybrid selection metrics
    $worksheet->write($num_files, 17, "$enrichment_kit");
    $worksheet->write($num_files, 18, "$genome_size");
    $worksheet->write($num_files, 19, "$bait_ter");
    $worksheet->write($num_files, 20, "$tar_ter");
    if ($num_pf_reads ne "-"){
        if ($num_pf_reads == $pf_reads){
            $worksheet->write($num_files, 21, "$num_pf_reads");
        }else{
            $worksheet->write($num_files, 21, "$num_pf_reads", $red);
        }
    }else{
        $worksheet->write($num_files, 21, "$num_pf_reads");
    }
    $worksheet->write($num_files, 22, "$num_pf_unique_reads");
    $worksheet->write($num_files, 23, "$pct_pf_unique_reads_aligned");
    $worksheet->write($num_files, 24, "$on_bait_bases");
    $worksheet->write($num_files, 25, "$near_bait_bases");
    $worksheet->write($num_files, 26, "$off_bait_bases");
    $worksheet->write($num_files, 27, "$on_target_bases");
    $worksheet->write($num_files, 28, "$pct_off_bait");
    $worksheet->write($num_files, 29, "$mean_bait_coverage");
    $worksheet->write($num_files, 30, "$mean_target_coverage");
    $worksheet->write($num_files, 31, "$pct_usable_bases_on_target");
    $worksheet->write($num_files, 32, "$pct_zero_cov_targets");
    $worksheet->write($num_files, 33, "$pct_target_bases_2x");
    $worksheet->write($num_files, 34, "$pct_target_bases_10x");
    if ($pct_target_bases_20x ne "-"){
        if ($pct_target_bases_20x >= 80){
            $worksheet->write($num_files, 35, "$pct_target_bases_20x");
        }elsif ($pct_target_bases_20x >=75){
            $worksheet->write($num_files, 35, "$pct_target_bases_20x", $orange);
        }else{
            $worksheet->write($num_files, 35, "$pct_target_bases_20x", $red);
        }
    }else{
        $worksheet->write($num_files, 35, "$pct_target_bases_20x");
    }
    $worksheet->write($num_files, 36, "$pct_target_bases_30x");
    
    #Write concordance metrics
    $worksheet->write($num_files, 38, "$num_snps");
    $worksheet->write($num_files, 39, "$pct_dbsnp");
    $worksheet->write($num_files, 40, "$titv_known");
    $worksheet->write($num_files, 41, "$titv_novel");
    $worksheet->write($num_files, 42, "$all_comp_het_call_het");
    $worksheet->write($num_files, 43, "$known_comp_het_call_het");
    $worksheet->write($num_files, 44, "$nonref_sens");
    $worksheet->write($num_files, 45, "$nonref_disc");
    if ($overall_conc ne "-"){
        if ($overall_conc >= 95){
            $worksheet->write($num_files, 46, "$overall_conc");
        }elsif ($overall_conc >= 70){
            $worksheet->write($num_files, 46, "$overall_conc", $orange);
        }else{
            $worksheet->write($num_files, 46, "$overall_conc", $red);
        }
    }else{
        $worksheet->write($num_files, 46, "$overall_conc");
    }

    
    $num_files++;
    undef %hash;
    undef @array;
}
