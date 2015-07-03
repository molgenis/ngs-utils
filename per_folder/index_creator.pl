#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use Cwd;

#############################
##USAGE:
#############################

my $dir;
my @files;
my $file;
my $bwa_dir;
my $samtools_dir;
my $cwd;
my $bwa_index;
my $samtools_index;

##RETRIEVE DIRECTORY IN WHICH FILES ARE LOCATED##
$cwd = getcwd();
opendir DIR, $cwd or die "cannot open directory named $cwd: $!";
@files= grep { $_ ne '.' && $_ ne '..' && $_ ne '.pac' && $_ ne 'index_creator.pl'} readdir DIR;
closedir(DIR);

#@files = ("index_hg19.fa");

##BWA and SAMtools directory##
$bwa_dir = "/data/p255198/tools/bwa-0.5.8a";
$samtools_dir="/data/p25198/tools/samtools-0.1.8";


##Create indices##
foreach $file (@files){
    #print $file,"\n";
    $bwa_index="./bwa index -p $cwd/$file -a is $cwd/$file";
    chdir($bwa_dir);
    #print $bwa_index . "\n";
    system($bwa_index); #Execute command
    
    $samtools_index="./samtools faidx $cwd/$file";
    chdir($samtools_dir);
    #print $samtools_index . "\n";
    system($samtools_index);
    print "######################## Finished creating indices for $file ########################\n\n";
    #/data/home/data/sequences/pipeline/bwa-0.5.8a# ./bwa index -p chr1.fa -a is ../../cluster_test/indices/hg19/chr1.fa
}