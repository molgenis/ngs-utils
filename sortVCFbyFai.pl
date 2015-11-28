#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use File::Glob ':glob';


my ($help, $fai_file, $vcf_file, $outputvcf_file);

#### get options
GetOptions(
                "h"                                 => \$help,
                "fastaIndexFile=s"                  => \$fai_file,
                "inputVCF=s"                        => \$vcf_file,
                "outputVCF=s"                       => \$outputvcf_file  
          );
usage() and exit(1) if $help;
# mandatory args
usage() and exit(1) unless $fai_file;
usage() and exit(1) unless $vcf_file;
usage() and exit(1) unless $outputvcf_file;

chomp $fai_file;
chomp $vcf_file;
chomp $outputvcf_file;

print "\nStarted sorting $vcf_file ...\n";

open (OUTPUT, ">", "$outputvcf_file" ) or die $!;

#Load fa.fai file into memory
open(DICT,$fai_file) or die "Can't open $fai_file!\n";
my @contig_order;
my $c=0;
while(<DICT>){
    my ($contig) = $_ =~ /(^.+)\t[0-9]{1,}\t[0-9]{1,}\t[0-9]{1,}\t[0-9]{1,}/;
    $contig_order[$c]=$contig;
    ++$c; 
    #print $contig,"\n";
}
close(DICT);

#Parse VCF file and sorted output
open(VCF,$vcf_file) or die "Can't open $vcf_file!\n";

my %vcf_hash;
my $header;

while(<VCF>){
    if($_=~/^#/){ $header .= $_; } # store header and comment fields
    else{
    chomp($_);

    my @data = split(/\t/,$_);
    #print $_ . "\n";
    my $contig = $data[0];
    my $start = $data[1];
    my $variant = $data[3]."to".$data[4];
    my $line = $_;

    #print $contig,":",$start,"\n";
    
    $vcf_hash{$contig}{$start}{$variant}=$line;
    }

}
close(VCF);

#Print VCF in order of fa.fai file

#Print standard VCF header
print OUTPUT $header;

foreach my $contig (@contig_order) # sort by contig order
    {
    #print $contig,"\n";
    foreach my $start (sort {$a <=> $b} keys %{$vcf_hash{$contig}}) # sort numerically by coordinates
	{
	#print $start,"\n";
	foreach my $variant (keys %{$vcf_hash{$contig}{$start}}) # if overlapping mutation, print each variant
            {
            print OUTPUT $vcf_hash{$contig}{$start}{$variant},"\n";
            }	
	}	
    }

close(OUTPUT);

print "\n\nFinished sorting file and wrote output to $outputvcf_file\n";

sub usage {
        print <<EOF;
#########################################################################################
This script sorts a VCF file based on the sequence of chromosomes/contigs provided in a
fasta index file (*.fa.fai).
#########################################################################################
Usage: ./sortVCFbyFai.pl
\t-fastaIndexFile       The Fasta index file (*.fa.fai) to sort on.
\t-inputVCF             Input VCF file to sort.
\t-outputVCF            Fasta reference sorted VCF file.
#########################################################################################
EOF
 
}
