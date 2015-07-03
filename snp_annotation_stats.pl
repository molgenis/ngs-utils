#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use File::Basename;
use Switch;

######
#Example: perl snp_annotation_stats.pl -vcf_table input.vcf.table -typefile output.type.txt -classfile output.class.txt -impactfile output.impact.txt
######

my $totalsnps = 0;

#Retrieve input and create three output files
my ($vcf_table, $typefile, $classfile, $impactfile, $snptypes, $snpclasses, $snpimpacts, $help);

GetOptions(
                "h"            => \$help,
                "vcf_table=s"  => \$vcf_table,
                "typefile=s"   => \$typefile,
                "classfile=s"  => \$classfile,
                "impactfile=s" => \$impactfile,
                "snptypes=s"   => \$snptypes,
                "snpclasses=s" => \$snpclasses,
                "snpimpacts=s" => \$snpimpacts,
          );  
usage() and exit(0) if $help;
usage() and exit(0) unless $vcf_table;


open (INPUT,"<$vcf_table") or die "Can't open $vcf_table$!";
open (TYPE, ">$typefile") or die "Can't open $typefile$!";
open (CLASS, ">$classfile") or die "Can't open $classfile$!";
open (IMPACT, ">$impactfile") or die "Can't open $impactfile$!";

my %func_type;
my %func_class;
my %func_impact;
my %snp_types;
my %snp_classes;
my %snp_impacts;

#Initialize types, classes and impacts
my @types = split(',', $snptypes);
foreach my $elem (@types){
    $func_type{$elem} = 0;
}
my @classes = split(',', $snpclasses);
foreach my $elem (@classes){
    $func_class{$elem} = 0;
}
my @impacts = split(',', $snpimpacts);
foreach my $elem (@impacts){
    $func_impact{$elem} = 0;
}

#Read header vcf.table
my $head = `head -1 $vcf_table`;
chomp $head;
my @headers = split('\t',$head);

#Retrieve column indices
my %headers_to_indexes = (map { $headers[$_] => $_ } (0 .. $#headers));


#Read vcf.table per line
my $lines;
while ($lines = <INPUT>){
    chomp $lines;
    if ($lines !~ m/^Chromosome\tPosition.+/gs){ #NOT match headerline
        $totalsnps++; #Count total # SNPs
        #print "$lines\n";
        my @info = split('\t', $lines); #Retrieve column indices for columnnames below
        my $type = $info[$headers_to_indexes{ "snpEff_effect" }];
        my $class = $info[$headers_to_indexes{ "snpEff_functional_class" }];
        my $impact = $info[$headers_to_indexes{ "snpEff_impact" }];
        #print "$type\t$class\t$impact\n";
        if (exists($func_type{$type})){                                      #If step not exists push into hash
            # if the key is found in the hash come here
            $func_type{$type}++;
        }else{
            # come here if the key is not found in the hash
            $func_type{$type} = 1;
        }
        if (exists($func_class{$class})){                                      #If step not exists push into hash
            # if the key is found in the hash come here
            $func_class{$class}++;
        }else{
            # come here if the key is not found in the hash
            $func_class{$class} = 1;
        }
        if (exists($func_impact{$impact})){                                      #If step not exists push into hash
            # if the key is found in the hash come here
            $func_impact{$impact}++;
        }else{
            # come here if the key is not found in the hash
            $func_impact{$impact} = 1;
        }
    }
}

#Sort hashes and write output to tab-delimited file
my $typekey;
foreach $typekey (sort (keys(%func_type))) {                                     #Loop through hash sorting by key
    my $percentage = ($func_type{$typekey}/$totalsnps)*100;
    my $typekey1 = uc(substr($typekey, 0,1));
    my $typekey2 = lc(substr($typekey, 1));
    my $merged_typekey = "$typekey1" . "$typekey2";
    $merged_typekey =~ s/Non_synonymous/Non-synonymous/g;
    $merged_typekey =~ s/Splice_site/Splice-site/g;
    $merged_typekey =~ s/_/ /g;
    $merged_typekey =~ s/Utr 3 prime/3'UTR/g;
    $merged_typekey =~ s/Utr 5 prime/5'UTR/g;
    my $perc = sprintf("%.2f", $percentage);
    #print TYPE "$merged_typekey & $func_type{$typekey} & $perc";
    print TYPE "$merged_typekey\t$func_type{$typekey}";
    #print TYPE " \\" . "\\" . "\n";                                       #Print all values per key
    print TYPE "\n";
}

my $classkey;
foreach $classkey (sort (keys(%func_class))) {                                     #Loop through hash sorting by key
    my $percentage = ($func_class{$classkey}/$totalsnps)*100;
    my $classkey1 = uc(substr($classkey, 0,1));
    my $classkey2 = lc(substr($classkey, 1));
    my $perc = sprintf("%.2f", $percentage);
    #print CLASS "$classkey1" . "$classkey2 & $func_class{$classkey} & $perc";
    print CLASS "$classkey1" . "$classkey2\t$func_class{$classkey}";
    print CLASS "\n";                                       #Print all values per key
}

my $impactkey;
foreach $impactkey (sort (keys(%func_impact))) {                                     #Loop through hash sorting by key
    my $percentage = ($func_impact{$impactkey}/$totalsnps)*100;
    my $impactkey1 = uc(substr($impactkey, 0,1));
    my $impactkey2 = lc(substr($impactkey, 1));
    my $perc = sprintf("%.2f", $percentage);
    #print IMPACT "$impactkey1" . "$impactkey2 & $func_impact{$impactkey} & $perc";
    print IMPACT "$impactkey1" . "$impactkey2\t$func_impact{$impactkey}";
    print IMPACT "\n";                                       #Print all values per key
}

print TYPE "Total\t$totalsnps\n";
print CLASS "Total\t$totalsnps\n";
print IMPACT "Total\t$totalsnps\n";


########SUBS#######
sub usage {
        print <<EOF;
#########################################################################################
This script takes an input.vcf.table (generated with snpEff v2.0.5 ONLY) and extracts
the functional information per SNP. Per type, class and impact a LATEX ready output
file is generated.
#########################################################################################
Usage: ./snp_annotation_stats.pl
\t-vcf_table    input.vcf.table
\t-typefile     output file containing counts per type
\t-classfile    output file containing counts per class
\t-impactfile   output file containing counts per impact
\t-snptypes     list (comma separated) of snptypes to count and output
\t-snpclasses   list (comma separated) of snpclasses to count and output
\t-snpimpacts   list (comma separated) of snpimpacts to count and output
#########################################################################################
EOF
 
}