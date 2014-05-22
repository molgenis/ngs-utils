#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use Getopt::Long;


my ($help, $vcffile, $snpfile, $phasing, $output);

#### get options
GetOptions(
                "h"              => \$help,
                "vcffile=s"      => \$vcffile,
                "snpfile=s"      => \$snpfile,
                "phasing=s"      => \$phasing,
                "outputfile=s"   => \$output
          );
usage() and exit(1) if $help;
# mandatory args
usage() and exit(1) unless $vcffile;
usage() and exit(1) unless $snpfile;
usage() and exit(1) unless $phasing;
usage() and exit(1) unless $output;

my $symbol;
if ($phasing eq "phased"){
    $symbol = "|";
}elsif ($phasing eq "unphased"){
    $symbol = "/";
}else{
    print "ERROR: Incorrect value specified for \"phasing\"!\n";
    exit(1);
}

print "Reading VCF file...\n";

open (VCF, "<", $vcffile) or die $!;
open (OUTPUT, ">", $output) or die $!;
my $lines;
my @vcf;
while ($lines = <VCF>){
    chomp $lines;
    if ($lines =~ m/^#.+/gs){
        #Write header away
        print OUTPUT "$lines\n";
    }else{
        push(@vcf, "$lines\n");
    }
}

print "Reading SNP file...\n";

open (LIST, "<", $snpfile) or die $!;
my @list = <LIST>;
close (LIST);

my $idx_last_item = $#list;

print "Starting allele swapping...\n";
for (my $i=0; $i<=$idx_last_item; $i++){
    my $snpline = $list[$i];
    chomp $snpline;
    my $vcfline = $vcf[$i];
    chomp $vcfline;
    
    my @snparray = split ("\\s",$snpline);
    my $listRsID = $snparray[0];
    my $listRefAll = $snparray[1];
    my $listAltAll = $snparray[2];
    
    my @vcfarray = split ("\t", $vcfline);
    my $vcfline_idx_last_item = $#vcfarray;
    my $vcfRsID = $vcfarray[2];
    my $vcfRefAll = $vcfarray[3];
    my $vcfAltAll = $vcfarray[4];
    
    if ($listRsID eq $vcfRsID){
        if ($listRefAll eq $vcfRefAll && $listAltAll eq $vcfAltAll){
            #print normal line to outputfile
            print OUTPUT $vcfarray[0] . "\t" . $vcfarray[1] . "\t" . $vcfarray[2] . "\t" . $vcfarray[3] . "\t" . $vcfarray[4] . "\t" . $vcfarray[5] . "\t" . $vcfarray[6] . "\t" . $vcfarray[7] . "\t" . $vcfarray[8];
            for (my $j = 9; $j<= $vcfline_idx_last_item; $j++){
                my @indvInfo = split (";",$vcfarray[$j]);
                my $to_swap = $indvInfo[0];
                $to_swap =~ s/\|/$symbol/g;
                $to_swap =~ s/\//$symbol/g;
                my $out_genotypes = $to_swap;
                print OUTPUT "\t$out_genotypes";
                my $indvInfo_idx_last_item = $#indvInfo;
                for (my $l = 1; $l<= $indvInfo_idx_last_item; $l++){
                    print OUTPUT ";" . $indvInfo[$l];
                }
            }
            print OUTPUT "\n";
        }elsif ($listRefAll eq $vcfAltAll || $listAltAll eq $vcfRefAll){#Check if swapping of alleles is necessary
            print OUTPUT $vcfarray[0] . "\t" . $vcfarray[1] . "\t$vcfRsID\t$listRefAll\t$listAltAll\t" . $vcfarray[5] . "\t" . $vcfarray[6] . "\t" . $vcfarray[7] . "\t" . $vcfarray[8];
            for (my $k = 9; $k<= $vcfline_idx_last_item; $k++){
                my @indvInfo = split (";",$vcfarray[$k]);
                my $to_swap = $indvInfo[0];
                my $out_genotypes;
                if ($to_swap eq "0/0" || $to_swap eq "0|0"){
                    $out_genotypes = "1" . $symbol . "1";
                }elsif ($to_swap eq "1/1" || $to_swap eq "1|1"){
                    $out_genotypes = "0" . $symbol . "0";
                }elsif ($to_swap eq "0/1" || $to_swap eq "0|1"){
                    $to_swap =~ s/\|/$symbol/g;
                    $to_swap =~ s/\//$symbol/g;
                    $out_genotypes = "$to_swap";
                }elsif ($to_swap eq "./."){
                    $out_genotypes = "$to_swap";
                }else{
                    print "ERROR: Undefined genotype for SNP $vcfRsID!\n";
                    exit (1);
                }
                print OUTPUT "\t$out_genotypes";
                my $indvInfo_idx_last_item = $#indvInfo;
                for (my $l = 1; $l<= $indvInfo_idx_last_item; $l++){
                    print OUTPUT ";" . $indvInfo[$l];
                }
            }
            print OUTPUT "\n";
        }else{
            print "ERROR: Alleles of snp $vcfRsID in snpfile do not match alleles in VCF file!\n";
            exit (1);
        }
    }else{
        print "ERROR: rsIDs in snpfile not sorted in the same order as SNPs in VCF file!\n";
        exit (1);
    }
}

print "\nFinished swapping alleles!\n";


sub usage {
        print <<EOF;
#########################################################################################
swap_alleles_vcfV1.pl (version 1)

This script swaps alleles in a VCF file based on an input list of SNPs. The phase of the
genotypes can be changed as well.
#########################################################################################
Usage: ./swap_alleles_vcfV1.pl
\t-vcffile              Input VCF file to swap alleles in
\t-snpfile              Input snp file to swap alleles with
\t                      Format should be (rsID ReferenceAllele AlternativeAllele)
\t-phasing              Phasing of the output genotypes
\t                      This can be "phased" or "unphased"
\t-outputfile           Output file.
#########################################################################################
EOF

}