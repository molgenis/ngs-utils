#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use File::Basename;

####
# unfold_lanes_worksheet.pl
# WBK 20120829
####
# Script to unfold worksheet entries where lane is "1,2" or whatever lane numbers.
# Each such entry is duplicated, one version gets 1 and the other 2 for lane, like you would manually.
# input: GAF list or worksheet subset, input_worksheet.csv
# output: output_worksheet.csv or, if not given, input_worksheet_unfolded.csv 
# help: usage and exit.
####

#### get options
my ($help, $input_worksheet);
my $output_worksheet = 0;
GetOptions(
	'h|help'=> \$help,
        'i=s'   => \$input_worksheet,
        'o=s'   => \$output_worksheet
);

#### parse options
usage() and exit(0) if $help;
# mandatory argument: input
usage() and exit(0) unless $input_worksheet;
# optional argument: output
my ($input_worksheet_basename, $input_worksheet_dirs, $input_worksheet_suffix) = fileparse($input_worksheet, qr|\.[^.]*|);
unless ($output_worksheet) {$output_worksheet = $input_worksheet_dirs . $input_worksheet_basename . "_unfolded" . $input_worksheet_suffix};

#### read and write worksheets
open (INPUT,"<$input_worksheet") or die "Can't open $input_worksheet [$!]";
open (OUTPUT, ">$output_worksheet") or die "Can't open $output_worksheet [$!]";
while(<INPUT>) {
	if (m/,internalSampleID,externalSampleID,/) { # header
		print OUTPUT $_;
		next;
	} elsif (m/,\"([1-8]{1}),([1-8]{1})\",/) { # double lane entry
		my ($first, $second) = ($_, $_);
		$first =~ s/,\"([1-8]{1}),([1-8]{1})\",/,$1,/;
		$second =~ s/,\"([1-8]{1}),([1-8]{1})\",/,$2,/;
		print OUTPUT $first . $second;
	} else { # other entries
		print OUTPUT $_;
	}
}
close INPUT;
close OUTPUT;

#### write log message to STDOUT
print "Done, unfolded worksheet written to " . $output_worksheet . "\n";

#### sub: usage
# prints info message and usage
sub usage {
        print <<EOF;
#######################################################################################################
# Script to unfold worksheet entries where lane is "1,2" or whatever lane numbers.
# Each such entry is duplicated, one version gets 1 and the other 2 for lane, like you would manually.
# input: GAF list or worksheet subset, input_worksheet.csv
# output: output_worksheet.csv or, if not given, defaults to input_worksheet_unfolded.csv 
#######################################################################################################
Usage: ./unfold_lanes_worksheet.pl -i input_worksheet.csv -o output_worksheet.csv
\t-i\t\tInput GAF list or worksheet, mandatory
\t-o\t\tOutput worksheet, optional
\t-h, --help\tThis message
#######################################################################################################
EOF
} 
