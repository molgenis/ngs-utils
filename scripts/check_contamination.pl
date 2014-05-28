#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;

####Usage
##check_contamination.pl <sample> <*.self.merged> <*.best.merged>
####

my $sample = $ARGV[0];
my $selfmerged = $ARGV[1];
chomp $selfmerged;
my $bestmerged = $ARGV[2];
chomp $bestmerged;
my $selflines;
my $bestlines;
my $node;
my @array;

open (SELF,"<$selfmerged") or die "Can't open $selfmerged$!"; #Open file
open (BEST,"<$bestmerged") or die "Can't open $bestmerged$!"; #Open file
#my $first = <INPUT>;

my $path = "\n\n\nPath to files:\n$selfmerged\n$bestmerged\n";

my $documentation = "\n#Each sample or lane can be checked in this way. When [SELF-IBD] < 1 AND [%MIX] > 0 AND [REF-A2%] > 0.01, meaning 1% or more of non-reference bases are observed in reference sites, we recommend to examine the data more carefully for the possibility of contamination.
#We recommend to check each lane for the possibility of sample swaps. When [SELF-IBD] << 1 AND [%MIX] ~ 0, then it is possible that the sample is swapped with another sample. When [BEST-IBD] ~ 1, [BEST_SM] might be actually the swapped sample. Otherwise, the swapped sample may not exist in the genotype data you have compared against.
#When genotype data is not available but allele-frequency-based estimates of [%MIX] >= 0.03 and [BESTMIXLLK-] is large (greater than 100), then it is possible that the sample is contaminated with other sample. We recommend to use per-sample data rather than per-lane data for checking this for low coverage data, because the inference will be more confident when there are large number of bases with depth 2 or higher.
#If [EXHET] << 1.00, it indicates that the sequence reads have excessive homozygosity, and this might probably be due to the small library size. In this case, we recommend to check the duplication rate of the BAM file, for example, using samtools flagstat, and see whether the library size is too small to call heterozygous sites appropriately.";

while ($selflines=<SELF>){ #Read concatenated VerifyBamID SELF file per line
    chomp $selflines;
    my @columns = split('\t', $selflines);
        #my %column_element_index;
        #@column_element_index{ @columns } = (0 .. $#columns);
        #my $selfibd = $column_element_index{ 'SELFIBD' }; #-
        #my $mix = $column_element_index{ '%MIX' };         #-
        #my $ref = $column_element_index{ 'REF-A2%' };      #-All needed for conditional check. insert statement later
        #my $bestibd = $column_element_index{ 'BEST-IBD' };
        #my $bestsm = $column_element_index{ 'BEST-SM' };
        #my $bestmixllk = $column_element_index{ 'BESTMIXLLK-' };
        #my $exhet = $column_element_index{ 'EXHET' };
        #print "$lines\n";
    if ($selflines =~ m/$sample.+/gs){ #Check for sample lines and retrieve values from colums
        my $selfibd = $columns[3];
        my $mix = $columns[24];
        my $ref = $columns[15];
        my $exhet = $columns[23];
        my $bestmixllk_;
        
        $selfibd =~ s/N\/A/1/g; #Substitute N/A to 1
        $ref =~ s/N\/A/1/g;
        
        #if ($selfibd != 1 && $ref != 1){ #Check for possible contamination
        if ($selfibd < 1 && $mix > 0 && $ref > 0.01){ #Check for possible contamination
            my $contamination = "Possible contamination in sample: $sample\nThe following results were found:\nSELF-IBD: $selfibd\n%MIX: $mix\nREF-A2%: $ref\n";
            if ($selfibd < 1 && $mix < 0.1){ #Check for sample swaps
                while ($bestlines=<BEST>){ #Retrieve sample lines and values from columns from VerifyBamID BEST file
                    my @bestcolumns = split('\t',$bestlines);
                    if ($selflines =~ m/$sample.+/gs){
                        my $bestibd = $bestcolumns[3];
                        my $best_sm = $bestcolumns[2];
                        $bestmixllk_ = $bestcolumns[26];
                        my $sample_swap = ();
                        #print "BESTIBD: $bestibd\n";
                        if ($bestibd >= 0.9 && $bestibd <= 1.1){ #Check for possible sample swap and give best option for swap
                            $sample_swap = "\nPossible sample swap detected!\nThis swap might have occured between sample $sample and $best_sm\n";
                        }else{ #Suggest to check manually for sample swap (swaps in "grey zone")
                            $sample_swap = "\nNo sample swap detected, please manually check for this\n";
                        }
                        my $message_to_send = $contamination . $sample_swap . $path . $documentation; #Compose message from detected contamination, swaps etc.
                        #print "$message_to_send";
                        sendEmail("freerk.van.dijk\@gmail.com", "verifyBAMid\@go.nl", "Possible contamination in sample $sample", "$message_to_send"); #Send e-mail to fvandijk
                        sendEmail("laurent.francioli\@gmail.com", "verifyBAMid\@go.nl", "Possible contamination in sample $sample", "$message_to_send"); #Send e-mail to lfrancioli
                        exit;
                    }
                }
            } ##INTEGRATE CORRECT MAILADRESSES FOR IN-HOUSE SAMPLE USE
        }elsif ($mix >= 0.03 && $bestmixllk_ > 100){ #Check for contamination with another sample
            my $without_genotype_data_message = "#######################\nOnly take action if no genotype data was used during the verifyBAMid analysis!\n#######################\n\nSample $sample has possibly been contaminated with another sample\nThe following results were found:\n%MIX: $mix\nBESTMIXLLK-: $bestmixllk_\n";
            #sendEmail("freerk.van.dijk\@gmail.com", "verifyBAMid\@inhouse.nl", "Possible contamination in sample $sample", "$without_genotype_data_message");
            #sendEmail("m.dijkstra.work\@gmail.com", "verifyBAMid\@inhouse.nl", "Possible contamination in sample $sample", "$without_genotype_data_message");
        }elsif ($exhet < 0.9){
            my $exces_het_message = "Sequence reads from sample $sample have excessive homozygosity\nThe following results were found:\nEXHET: $exhet\n";
            #sendEmail("freerk.van.dijk\@gmail.com", "verifyBAMid\@inhouse.nl", "Possible contamination in sample $sample", "$exces_het_message");
            #sendEmail("m.dijkstra.work\@gmail.com", "verifyBAMid\@inhouse.nl", "Possible contamination in sample $sample", "$exces_het_message");
        }
    }
}

#Sub used to send e-mail
sub sendEmail
{
    my ($to, $from, $subject, $message) = @_;
    my $sendmail = '/usr/lib/sendmail';
    open(MAIL, "|$sendmail -oi -t");
        print MAIL "From: $from\n";
	print MAIL "To: $to\n";
	print MAIL "Subject: $subject\n\n";
	print MAIL "$message\n";
    close(MAIL);
}