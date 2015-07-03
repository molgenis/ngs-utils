#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use File::Glob ':glob';
use File::Basename;

my ($help, $studyId, $studyDir, $pseudoFile, $testDir, $refStudyDir, $report);

#Get options
GetOptions(
                "h"                                     => \$help,
                "studyId=s"                             => \$studyId,
                "studyDir=s"                            => \$studyDir,
                "pseudoFile=s"                          => \$pseudoFile,
                "testDir=s"                             => \$testDir,
                "refStudyDir=s"                         => \$refStudyDir,
                "report=s"                              => \$report,
          );
usage() and exit(1) if $help;
#Mandatory args
usage() and exit(1) unless $studyId;
usage() and exit(1) unless $studyDir;
usage() and exit(1) unless $pseudoFile;
usage() and exit(1) unless $testDir;
usage() and exit(1) unless $refStudyDir;
usage() and exit(1) unless $report;

#Chomp inputs
chomp $studyId;
chomp $studyDir;
chomp $pseudoFile;
chomp $testDir;
chomp $refStudyDir;
chomp $report;

print "\nStarting tests, please be patient..\n\n";

###Check if directories exist
if ( !-d "$studyDir"){
        print "$studyDir does not exist, please check your input.\n";
        exit (1);
}
if ( -d "$testDir") { #If testDir exists, empty previous one. Else create new one
        print "$testDir already exists, emptying it now.\n";
        `rm $testDir/*`;
}else {
        print "$testDir does not yet exist, creating it now.\n";
        `mkdir $testDir`;
}

#Check if files exist
if ( !-e "$pseudoFile"){
        print "$pseudoFile does not exist, please check your input.\n";
        exit (1);
}

#Create sorted file containing original pseudoFile IDs
`awk '{print \$2}' $pseudoFile | tr -d \'\015\' | sort > $testDir/original.ids`;


#Glob directories for files to use
my @studyDoseFiles = glob "$studyDir/*.dose";
my @pedFiles = glob "$studyDir/*.ped";
my @mapFiles = glob "$studyDir/*.map";
my @refDoseFiles = glob "$refStudyDir/*.dose";
my @refMapFiles = glob "$refStudyDir/*.map";

print "\nCreating files needed for tests in $testDir\n\n";

###Extract info from generated study files
#Create *.ids and *.snps files from study dosage files to use in test
print "Creating study dosefiles..\n";
for my $doseFile (@studyDoseFiles){
        print "$doseFile\n";
        my $doseFileName = basename $doseFile;
        #Alter Beagle imputed dosage file
        #`head -n 1 $doseFile | sed \'s/SNP\\tA1\\tA2\\t1\\t//g\' | sed \'s/\\t1\\t/\\n/g\' | tr -d \'\015\' > $testDir/$doseFileName.ids`;
        #Alter Minimac imputed dosage file
        `head -n 1 $doseFile | sed \'s/SNP\\tA1\\tA2\\t1\\s//g\' | sed \'s/\\t1\\s/\\n/g\' | tr -d \'\015\' > $testDir/$doseFileName.ids`;
        `awk '{ print \$1 }' $doseFile | sed 1d > $testDir/$doseFileName.snp`;
}
print "Creating study IDfiles..\n";
#Create *.ids file from PED file
for my $pedFile (@pedFiles){
        print "$pedFile\n";
        my $pedFileName = basename $pedFile;
        `awk '{print \$2}' $pedFile > $testDir/$pedFileName.ids`;
}
print "Creating study snpfiles..\n";
#Create *.snp file from MAP file
for my $mapFile (@mapFiles){
        print "$mapFile\n";
        my $mapFileName = basename $mapFile;
        `awk '{print \$2}' $mapFile > $testDir/$mapFileName.snp`;
}

print "Creating snpfiles from reference dosefiles..\n";
###Extract info from original reference data
#Create *.snp file from reference dose file
for my $refDoseFile (@refDoseFiles){
        print "$refDoseFile\n";
        my $refDoseFileName = basename $refDoseFile;
        `awk '{print \$1}' $refDoseFile | sed 1d > $testDir/$refDoseFileName.snp`;
}

print "Creating snpfiles from reference mapfiles..\n";
#Create *.snp file from reference map file
for my $refMapFile (@refMapFiles){
        print "$refMapFile\n";
        my $refMapFileName = basename $refMapFile;
        `awk '{print \$2}' $refMapFile > $testDir/$refMapFileName.snp`;
}


###Start test by looping through 22 chromosomes
for (my $i=1; $i < 23; $i++){
        #Create filenames to compare
        my $compareRefDoseSnp = "$testDir/output.chr$i.dose.snp";
        my $compareRefMapSnp = "$testDir/output.chr$i.map.snp";
        my $compareDoseFileSnp = "$testDir/$studyId\_chr$i.dose.snp";
        my $compareDoseFileIds = "$testDir/$studyId\_chr$i.dose.ids";
        my $comparePedFileIds = "$testDir/$studyId\_chr$i.ped.ids";
        my $compareMapFileSnp = "$testDir/$studyId\_chr$i.map.snp";
        print "Starting test for chromosome $i..\n";
        #Check if files exist
        if ( !-e "$compareRefDoseSnp"){
                print "$compareRefDoseSnp does not exist, exiting now!\n";
                exit (1);
        }
        if ( !-e "$compareRefMapSnp"){
                print "$compareRefMapSnp does not exist, exiting now!\n";
                exit (1);
        }
        if ( !-e "$compareDoseFileSnp"){
                print "$compareDoseFileSnp does not exist, exiting now!\n";
                exit (1);
        }
        if ( !-e "$compareDoseFileIds"){
                print "$compareDoseFileIds does not exist, exiting now!\n";
               exit (1);
        }
        if ( !-e "$comparePedFileIds"){
                print "$comparePedFileIds does not exist, exiting now!\n";
                exit (1);
        }
        if ( !-e "$compareMapFileSnp"){
                print "$compareMapFileSnp does not exist, exiting now!\n";
                exit (1);
        }
        ###Start tests
        #Test zero diff lines
        my $in1;
        my $in2;

        $in1 = $compareRefDoseSnp;
        $in2 = $compareRefMapSnp;
        my $cmd1 = "diff --side-by-side --suppress-common-lines $in1 $in2 | wc -l";
        my $diff = `$cmd1`;
        print "TESTING: $cmd1\nRESULT:$diff\n";
        if ($diff != 0){ #Difference between both files, exit.
                print "\nERROR:\nTest on chromosome $i failed.\nFailing test: $cmd1\n\n";
                exit (1);
        }
        $in1 = $compareRefMapSnp;
        $in2 = $compareDoseFileSnp;
        $cmd1 = "diff --side-by-side --suppress-common-lines $in1 $in2 | wc -l";
        $diff = `$cmd1`;
        print "TESTING: $cmd1\nRESULT:$diff\n";
        if ($diff != 0){ #Difference between both files, exit.
                print "\nERROR:\nTest on chromosome $i failed.\nFailing test: $cmd1\n\n";
                exit (1);
        }
        $in1 = $compareDoseFileSnp;
        $in2 = $compareMapFileSnp;
        $cmd1 = "diff --side-by-side --suppress-common-lines $in1 $in2 | wc -l";
        $diff = `$cmd1`;
        print "TESTING: $cmd1\nRESULT:$diff\n";
        if ($diff != 0){ #Difference between both files, exit.
                print "\nERROR:\nTest on chromosome $i failed.\nFailing test: $cmd1\n\n";
                exit (1);
        }
        $in1 = $compareDoseFileIds;
        $in2 = "$testDir/original.ids";
        my $cmd2 = "sort $in1 | diff - --side-by-side --suppress-common-lines $in2 | wc -l";
        my $diff2 = `$cmd2`;
        print "TESTING: $cmd2\nRESULT:$diff2\n";
        if ($diff2 != 0){ #Difference between both files, exit.
                print "\nERROR:\nTest on chromosome $i failed.\nFailing test: $cmd2\n\n";
                exit (1);
        }
        $in1 = $comparePedFileIds;
        $in2 = "$testDir/original.ids";
        $cmd2 = "sort $in1 | diff - --side-by-side --suppress-common-lines $in2 | wc -l";
        $diff2 = `$cmd2`;
        print "TESTING: $cmd2\nRESULT:$diff2\n";
        if ($diff2 != 0){ #Difference between both files, exit.
                print "\nERROR:\nTest on chromosome $i failed.\nFailing test: $cmd2\n\n";
                exit (1);
        }
        print "Done comparing chr$i\n\n\n";
}

print "\n+---------- ALL TESTS PASSED ----------+\n";

#Count individuals in test and print
my $indv = `cat $testDir/original.ids | wc -w`;

print "\nNumber of participants: $indv\n";

print "\n\nWriting report to $report\n\n";

###Create output report
#Open output report
open (REPORT, ">", $report ) or die $!;

print REPORT "-----------------------\n";
print REPORT "STARTING REPORT OUTPUT\n";
print REPORT "-----------------------\n";

print REPORT "\n\nInput parameters:\n";
print REPORT "studyId: $studyId\n";
print REPORT "studyDir: $studyDir\n";
print REPORT "pseudoFile: $pseudoFile\n";
print REPORT "testDir: $testDir\n";
print REPORT "refStudyDir: $refStudyDir\n";
print REPORT "report: $report\n";
print REPORT "\n\n";
print REPORT "Number of participants: $indv\n";
print REPORT "\n\nAll tests PASSED\n";
print REPORT "\n\n\n\n";
print REPORT "+----------+\n";
print REPORT "| TESTS OK |\n";
print REPORT "+----------+\n";
print REPORT "\n\n\n";
print REPORT "-----------------------\n";
print REPORT "ENDING REPORT OUTPUT\n";
print REPORT "-----------------------\n";

#Close output report filehandler
close (REPORT);

###Describe usage of script
sub usage {
        print <<EOF;
###############################################################################
This script performs several validation tests on lifelines data.
###############################################################################
Usage: ./validate.pl
\t-studyId                   Study ID.
\t-studyDir                  Directory which contains study output.
\t-pseudoFile                File containing sample mappings (pseudo file).
\t-testDir                   Directory to which tests are written.
\t-refStudyDir               Reference study directory used for imputation.
\t-report                    Output test report.
###############################################################################
EOF

}