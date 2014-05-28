#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use Cwd;

#############################
##USAGE:
#############################

my @files;
my @filtered;
my $file;
my $cwd;
my $command;
my $chr;
my $element;

##RETRIEVE DIRECTORY AND CONTENTS##
my $dir = "/data/fvandijk/resources/hg19/indices";
#$cwd = getcwd("/data/fvandijk/resources/hg19/indices");
opendir DIR, $dir or die "cannot open directory named $dir: $!";
@files= grep { $_ ne '.' && $_ ne '..' && $_ ne '.pac' && $_ ne 'index_creator.pl' && $_ ne '*.fa.*'} readdir DIR;
closedir(DIR);

#filter out everything except .fa files to generate indel intervals on#
######NOTE: hg19.fa and chr1.fa filtered out because indices for these two files were already excisting######
######Remove these files from regular expression if they need to be created#####
foreach $file (@files){
    if ($file =~ m/chr[0-9]{1,2}.fa.+/gs || $file =~ m/index.+/gs || $file =~ m/.+fai/gs || $file =~ m/.+dict/gs || $file =~ m/hg19.fa/gs || $file =~ m/chr1.fa/gs){
        #shift(@files);
        #skip this file
    }else{
        push(@filtered,$file);
    }
}

#####Generate interval file for each chromosome#####
foreach $element (@filtered){
    print $element . "\n";
    $chr = $element;
    $chr =~ s/.fa//mg;
    $command = "java -jar /data/fvandijk/tools/Sting/dist/GenomeAnalysisTK.jar -T RealignerTargetCreator " .
    "-R /data/fvandijk/resources/hg19/indices/" . $element .
    " -D /data/fvandijk/resources/hg19/dbsnp/dbsnp_129_b37_" . $chr . ".rod" .
    " -B:indels,VCF /data/fvandijk/resources/hg19/indels/1kg.pilot_release.merged.indels.sites.hg19." . $chr . ".vcf" .
    " -o /data/fvandijk/resources/hg19/intervals/realign_intervals_hg19_" . $chr . ".intervals";
    system($command);
    #print $command . "\n";
    print "Finished creating intervals for file: " . $element . "\n\n";
}
print "Finished creating intervals for all files!\n";