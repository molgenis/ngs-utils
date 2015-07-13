#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use File::Glob ':glob';

#############################
##USAGE: perl remove_bad_files.pl <working_directory> <run_report>
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

##Retrieve files containing QC information
my @files = glob "$input_dir/outputs/*.log";

##Open output directory
open (OUTPUT, ">$output_dir");

##Walk through log file paths
while ($_ = shift @files) {
    #print $_ . "\n";
    open (FILE,"<$_") or die "Can't open $_$!";
    if ($_ =~ m/.+(FC.+).1.log/gs){                                         #Retrieve log files
        my $file = $1;
        print OUTPUT $file;
    }
    my %hash=();
    while($logline = <FILE>){                                               #Retrieve all lines starting with Completed and push them into array.
        chomp $logline;                                                     #Bottomline is pushed in array first
        if($logline =~ m/Completed (.+).ftl .+in ([0-9]{1,}) seconds/gs){
            unshift (@array, $logline);                                     
        }
    }
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
}
