#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use File::Glob ':glob';

#############################
##USAGE: perl check_output2.pl <working_directory> <output_directory>
#############################

my $dir;

my $input_dir = $ARGV[0];
my $output_dir = $ARGV[1];
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

##Retrieve files containing QC information
my @files = glob "$input_dir/outputs/*.log";
my @steplogs = glob "$input_dir/outputs/FC*";

##Open output directory
open (OUTPUT, ">$output_dir");

##Print output header
$header = "file\tpe0--fastqc\tpe00-bwa-align-pair1\tpe01-bwa-align-pair2\tpe02-bwa-sampe\tpe03-sam-to-bam\tpe04-sam-sort\tpe04a-HsMetrics\tpe04b-picardQC\tpe05-mark-duplicates\tpe06-realign\tpe07-fixmates\tpe08-covariates-before\tpe09-recalibrate\tpe10-sam-sort\tpe11-covariates-after\tpe12-analyze-covariates\n";
print OUTPUT $header;

##Walk through log file paths
while ($_ = shift @files) {
    #print $_ . "\n";
    my %hash=();
    open (FILE,"<$_") or die "Can't open $_$!";
    if ($_ =~ m/.+(FC.+).1.log/gs){                                         #Retrieve log files
        my $file = $1;
        print OUTPUT $file;
    }

    while($logline = <FILE>){                                               #Retrieve all lines starting with Completed and push them into array.
        chomp $logline;                                                     #Bottomline is pushed in array first
        if($logline =~ m/Completed (.+).ftl .+in ([0-9]{1,}) seconds/gs){
            unshift (@array, $logline);                                     
        }
    }
    close (FILE);
    foreach my $element(@array){
        if($element =~ m/Completed (.+).ftl .+in ([0-9]{1,}) seconds/gs){   #Retrieve stepname and elapsed time from lines
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
    print OUTPUT "\n";
    undef %hash;
    undef @array;
}

##Retrieve nodes on which process ran

print OUTPUT "\n\n\n" . $header;

my @filter;
foreach my $element (@steplogs){
    if ($element =~ m/(.+FC.+).[0-9]{1,}.log/gs){
        my $var = $1;
        push(@filter,$var);
        #print "$var\n";
    }
}

#while ($_ = shift @steplogs) {
while ($_ = shift @filter) {
    #print $_ . "\n";
    if ($_ =~ m/.+\/outputs\/(FC.+_L[0-9])/gs){
        $folder = $1;
        #print "$folder\n";
        print OUTPUT $folder;
        my @logfiles = glob "$input_dir/outputs/$folder/logs/*.out";
        while ($_ = shift @logfiles){
            if ($_ =~ m/.+\/(pe[0-9]{1,}.+)out/gs){
                $stepname = $1;
            }
            open (LOGS,"<$_") or die "Can't open $_$!";
            while (my $linlog = <LOGS>){
                if ($linlog =~ m/.+set: (.+).cm.+/gs){
                    $node = $1;
                }
                $loghash{ $stepname } = $node;
            }
        }
        my $key;
        foreach $key (sort (keys(%loghash))) {                                     #Loop through hash sorting by key
            print OUTPUT "\t$loghash{$key}";                                       #Print all values per key
        }
        print OUTPUT "\n";
        undef %loghash;
        undef @logfiles;
    }
}

print OUTPUT "\n\n\n" . $header;


my @logfiles = glob "$input_dir/outputs/*.log";

while ($_ = shift @logfiles) {
    #print $_ . "\n";
    my %starthash=();
    open (LOGFILE,"<$_") or die "Can't open $_$!";
    if ($_ =~ m/.+(FC.+).1.log/gs){                                         #Retrieve log files
        my $filenam = $1;
        print OUTPUT $filenam;
    }
    while($logline = <LOGFILE>){                                            #Retrieve all lines starting with Completed and push them into array.
        chomp $logline;                                                     #Bottomline is pushed in array first
        #Begin pe01-bwa-align-pair2.ftl at Thu Dec 23 10:43:15 CET 2010
        if($logline =~ m/Begin (.+).ftl at.+/gs){
            unshift (@array, $logline);                                     
        }
    }
    close (LOGFILE);
        foreach my $element(@array){
        if ($element =~ m/Begin (.+).ftl at [A-Za-z]{3,3} (.+)/gs){   #Retrieve stepname and elapsed time from lines
            my $step = $1;
            my $starttime = $2;
            if (exists($starthash{$step})){                                      #If step not exists push into hash
                # if the key is found in the hash come here
            }else{
            # come here if the key is not found in the hash
                $starthash{ $step } = $starttime;
                }
            }
        }
    my $key;
    foreach $key (sort (keys(%starthash))) {                                     #Loop through hash sorting by key
        print OUTPUT "\t$starthash{$key}";                                       #Print all values per key
    }
    print OUTPUT "\n";
    undef %starthash;
    undef @array;
}