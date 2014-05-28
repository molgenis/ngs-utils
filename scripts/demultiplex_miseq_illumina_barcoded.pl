#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use File::Glob ':glob';
use File::stat;
use Getopt::Long;
use File::Basename;
use Switch;
use Text::CSV;
use Scalar::Util;

#perl ../demultiplex.pl -end1 test_1.fq -end2 test_2.fq -index index.fq -csvfile SampleSheetConverted.csv -output ~/miseq/testoutput


my ($help, $end1, $end2, $indexfile, $csvfile, $outputfolder);

#### get options
GetOptions(
                "h"         => \$help,
                "end1=s"    => \$end1,
                "end2=s"    => \$end2,
                "index=s"   => \$indexfile,
                "csvfile=s" => \$csvfile,
                "output=s"  => \$outputfolder  
          );
usage() and exit(1) if $help;
# mandatory args
usage() and exit(1) unless $end1;
usage() and exit(1) unless $end2;
usage() and exit(1) unless $indexfile;
usage() and exit(1) unless $csvfile;
usage() and exit(1) unless $outputfolder;

#Check if output dir is empty
print "Checking if output directory is empty..\n";
opendir(DIR, $outputfolder) or die "$!";
readdir DIR; # reads .
readdir DIR; # reads ..
if (readdir DIR) { #now there should be a file or folder
    print "Output directory is not empty!\n";
    print "Please manually clean the directory or choose a different output folder!\n";
    exit (1);
}
else {
   print "Output directory is empty!\n";
   print "Starting analysis now..\n";
}
close DIR;

print "Started!\n";

#Open files
open (END1, "<", $end1) or die("Unable to open file $end1");
open (END2, "<", $end2) or die("Unable to open file $end2");
open (INDEX, "<", $indexfile) or die("Unable to open file $indexfile");
 
#Read files into array
my @end1file = <END1>;
my @end2file = <END2>;
my @indexfile = <INDEX>;

#Close files
close(END1);
close(END2);
close(INDEX);

my $end1lastindex = $#end1file;

#read index file and retrieve indx from barcode column
#Read CSV file
my $csv = Text::CSV->new();
open (CSV, "<", $csvfile) or die $!;

my @headers;
while (<CSV>) {
    next if ($. != 1);
    if ($csv->parse($_)) {
        @headers = $csv->fields();
    } else {
        my $err = $csv->error_input;
        print "Failed to parse line: $err";
    }
}
close(CSV);

####Retrieve column indices
my %headers_to_indexes = (map { $headers[$_] => $_ } (0 .. $#headers));

my $sampleID;
my $samplebarcode;
my %sample_barcode_map;
open (CSV, "<", $csvfile) or die $!;
while (<CSV>) {
    next if ($. == 1);
    if ($csv->parse($_)) {
        my @columns = $csv->fields();
        $sampleID = $columns[$headers_to_indexes{"Sample_ID"}];
        $samplebarcode = $columns[$headers_to_indexes{"index"}];
        if (exists $sample_barcode_map{$sampleID}) {
            #Do nothing
        }else{
            #Add to hash
            $sample_barcode_map{ $sampleID } = $samplebarcode;
        }
    }
}
close(CSV);

#Print sample barcode combinations
print "Found the following sample barcode combinations:\n";
foreach my $key (sort keys %sample_barcode_map){
    print "$key\t$sample_barcode_map{$key}\n";
}

my $tot = ($end1lastindex+1)/4;
print "\n\nTotal number of reads found: $tot\n\n";

#Cut file extensions
my $prefix_out_end1 = $end1;
my $prefix_out_end2 = $end2;
$prefix_out_end1 =~ s/.fq//gs;
$prefix_out_end2 =~ s/.fq//gs;

#Actual code
for (my $i=0; $i <= $end1lastindex; $i = $i + 4){
    #print "$i\n";
    #print $end1file[$i] . $end1file[$i+1] . $end1file[$i+3] . "\n"
    #Retrieve information from read IDs
    my $end1idline = $end1file[$i];
    my $end1seq = $end1file[$i+1];
    my $end1qual = $end1file[$i+3];
    my $end1id;
    my $end1idx;
    my $end2idline = $end2file[$i];
    my $end2seq = $end2file[$i+1];
    my $end2qual = $end2file[$i+3];
    my $end2id;
    my $end2idx;
    my $indexidline = $indexfile[$i];
    my $indexid;
    my $indexidx;
    if ($end1idline =~ m/(.+) .+:.+:.+:(.+)/gs){
        $end1id = $1;
        $end1idx = $2;
        chomp $end1idx;
        #print "$end1id\t$end1idx\n";
    }
    if ($end2idline =~ m/(.+) .+:.+:.+:(.+)/gs){
        $end2id = $1;
        $end2idx = $2;
        chomp $end2idx;
    }
    if ($indexidline =~ m/(.+) .+:.+:.+:(.+)/gs){
        $indexid = $1;
        $indexidx = $2;
        chomp $indexidx;
    }
    #Check if all info matches between end1, end2 and index file
    if ($end1id eq $end2id && $end1id eq $indexid && $end2id eq $indexid && $end1idx eq $end2idx && $end1idx eq $indexidx && $end2idx eq $indexidx){
        #Retrieve barcode corresponding to id
        my $brcd = $sample_barcode_map{ "Sample$end1idx" };
        #print "$brcd\n";
        #Write output away        
        open (OUT1, ">>", "$outputfolder/$prefix_out_end1\_$brcd.fq") or die("Unable to open file $outputfolder/$prefix_out_end1\_$brcd.fq");
        open (OUT2, ">>", "$outputfolder/$prefix_out_end2\_$brcd.fq") or die("Unable to open file $outputfolder/$prefix_out_end2\_$brcd.fq");
        print OUT1 $end1idline . $end1seq . "+\n" . $end1qual;
        print OUT2 $end2idline . $end2seq . "+\n" . $end2qual;
        close(OUT1);
        close(OUT2);
        
    }else{
        print "Reads with ID $end1id don't match sequence ID or barcode ID\nExiting now!\n";
        exit(1);
    }
}

print "\nGathering statistics..\n";

#Count statistics
my @files1 = glob "$outputfolder/$prefix_out_end1*";
my @files2 = glob "$outputfolder/$prefix_out_end2*";

my $files1index = $#files1;
my $files2index = $#files2;

if ($files1index != $files2index){
    print "Number of 1st en 2nd end output files is not equal!\n";
    print "Please check if files are missing!\nExiting now!\n";
    exit(1);
}

open (LOG, ">", "$outputfolder/demultiplex_log.log") or die("Unable to open file $outputfolder/demultiplex_log.log");

#Print statistics per barcode
print LOG "\n\nFILE\tTOTAL_NUMBER_OF_READ_PAIRS\tPERCENTAGE_OF_TOTAL\n";
for (my $k=0; $k <= $files1index; $k++){
    my $first = `wc -l $files1[$k] | awk '{print \$1}'`;
    my $second = `wc -l $files2[$k] | awk '{print \$1}'`;
    chomp $first;
    chomp $second;
    #print "$first\t$second\n";
    if ($first == $second){
        #calculate statistics
        my $reads = ($first/4);
        my $perc = ($reads/$tot)*100;
        print LOG $files1[$k] . "\t$reads\t$perc\n";
    }else{
        print "Detected different number of lines in " . $files1[$k] . " and " . $files2[$k] . "\n";
        print "Exiting now!";
        exit(1);
    }
}

close(LOG);

print "\n\nFinished!\n";


sub usage {
        print <<EOF;
#########################################################################################
This script demultiplexes paired-end illumina barcoded ngs data generated by the MiSeq.
#########################################################################################
Usage: ./demultiplex_miseq_illumina_barcoded.pl
\t-end1         First end of pair to demultiplex
\t-end2         Second end of pair to demutliplex
\t-index        Third file containing the barcodes (I file)
\t-csvfile      CSV file containing the sample IDs and used barcodes
\t-output       Folder to write demultiplexed output to
#########################################################################################
EOF
 
}