#!/usr/bin/env perl
#
# This script takes a sample sheet in GAF *.csv format and creates one in Illumina *.csv format
# suitable for conversion of BCL to FastQ files using Illumina's bcl2fastq tool.
#
#
# =====================================================
# $Id: CreateIlluminaSampleSheet.pl 1054 2013-04-22 17:07:07Z pneerincx $
# $URL: http://www.bbmriwiki.nl/svn/ngs_scripts/trunk/scripts_data_archiving/CreateIlluminaSampleSheet.pl $
# $LastChangedDate: 2013-04-22 19:07:07 +0200 (Mon, 22 Apr 2013) $
# $LastChangedRevision: 1054 $
# $LastChangedBy: pneerincx $
# =====================================================
# 

use strict;
use warnings;
use Getopt::Std;
use File::Basename;
use Text::CSV;
use Log::Log4perl qw(:easy);

my %log_levels = (
	'ALL'   => $ALL,
	'TRACE' => $TRACE,
	'DEBUG' => $DEBUG,
	'INFO'  => $INFO,
	'WARN'  => $WARN,
	'ERROR' => $ERROR,
	'FATAL' => $FATAL,
	'OFF'   => $OFF,
);

my @gaf_column_names = ('internalSampleID','externalSampleID','sequencingStartDate','sequencer','run','flowcell','lane','barcode','barcodeType','seqType');
my %gaf_barcode_types = ('None' => {
							'IlluminaDemultiplexing'	=> 0},
						 'GAF' => {
							'IlluminaDemultiplexing'	=> 0},
						 'RPI' => {
							'IlluminaDemultiplexing'	=> 1},
						 'NEX' => {
                                                        'IlluminaDemultiplexing'        => 1},
						 'LEX' => {
                                                        'IlluminaDemultiplexing'        => 1},
						 'AGI' => {
							'IlluminaDemultiplexing'	=> 1},
						 'MON' => {
							'IlluminaDemultiplexing'	=> 1},
						 'AG8' => {
                                                        'IlluminaDemultiplexing'        => 1},
						);
my $log_section_break = '##################################################################################################################################';

#
# Get options.
#
my %opts;
Getopt::Std::getopts('i:o:r:l:', \%opts);

my $input		= $opts{'i'};
my $output		= $opts{'o'};
my $run_to_find	= $opts{'r'};
my $log_level	= $opts{'l'};

#
# Configure logging.
#
# Provides default if user did not specify log level:
$log_level = (defined($log_level) ? $log_level : 'INFO');

# Reset log level to default if user specified illegal log level.
$log_level = (
	defined($log_levels{$log_level})
	? $log_levels{$log_level}
	: $log_levels{'INFO'});

#Log::Log4perl->init('log4perl.properties');
Log::Log4perl->easy_init(
	{
		level  => $log_level,
		file   => "STDOUT",
		layout => '%d L:%L %p> %m%n'
	},
);
my $logger = Log::Log4perl::get_logger();

#
# Check user input.
#
unless (defined($input) && defined($output) && defined($run_to_find)) {
	_Usage();
	exit(1);
}
if ($input =~ /^$/ || $output =~ /^$/) {
	_Usage();
	exit(1);
}
if ($input eq $output) {
	_Usage();
	$logger->fatal('Output file is the same as the input file. Select a different file for the output.');
	exit(1);
}
unless (-f $input && -r $input) {
	$logger->fatal('Cannot read/access file: ' . $input);
	exit(1);
}
unless ($run_to_find =~ m/^[0-9]{4}$/) {
	_Usage();	
	$logger->fatal('Run number in unsupported format: must be a for digit integer like for example 0488.');
	exit(1);
}

#
# Check GAF list a.k.a the input samples sheet.
#
$logger->info('Checking input sample sheet in GAF *.csv format.');

my $input_fh;
eval {
	open($input_fh, '<', $input);
};
if ($@) {
	$logger->fatal('Cannot open CSV file: ' . $@);
	exit(1);
}

#
# Retrieve header from CSV file
#
my $csv = Text::CSV->new();
my @column_labels;
while (<$input_fh>) {
    next unless ($. == 1); # header only
    if ($csv->parse($_)) {
        @column_labels = $csv->fields();
    } else {
		$logger->fatal('Cannot parse CSV header line: ' . $csv->error_input);
		exit(1);
    }
    last;
}
close($input_fh);

#
# Check header.
#
# * Retrieve column indices.
# * Create hash with header name as key and 0-based column index as value.
# * Check if the required columns were present in the header line.
#
my %header = (map { $column_labels[$_] => $_ } (0 .. $#column_labels));
foreach my $column_name (@gaf_column_names) {
	unless (defined($header{$column_name})) {
		$logger->fatal('Missing column ' . $column_name . '. Please check if the input sample sheet in GAF *.csv format has a proper header.');
		exit(1);
	}
}

$logger->info($log_section_break);
$logger->info('Found the following samples:');
my $format = '%16s  %25.25s  %19s  %9s  %4s  %10s  %4s  %7s  %11s  %7s';
my $log_header = sprintf($format, @gaf_column_names);
$logger->info($log_header);

#
# Retrieve values from relevant columns using column indices.
# 
# Create experiment hash for each "measurement" a.k.a sample + lane + optional barcode combination.:
#   keys:   lanes
#   values: seqType, barcodeType, barcode
#
my %experiments; # keys: lanes, values: sequencetypes

$csv = Text::CSV->new();
eval {
	open($input_fh, '<', $input);
};
if ($@) {
	$logger->fatal('Cannot open CSV file: ' . $@);
	exit(1);
}

my @sequence_types;
# extract these columns:
my ($internalSampleID, $externalSampleID, $sequencingStartDate, $sequencer, $run, $flowcell, $lane, $barcode, $barcodeType, $seqType);
my $experiment_count = 0;

while (<$input_fh>) {
	next if ($. == 1); # skip header
	if ($csv->parse($_)) {
		my @fields = $csv->fields();
		$internalSampleID		= $fields[$header{'internalSampleID'}];
		$externalSampleID		= $fields[$header{'externalSampleID'}];
		$sequencingStartDate	= $fields[$header{'sequencingStartDate'}];
		$sequencer				= $fields[$header{'sequencer'}];
		$run					= $fields[$header{'run'}];
		$flowcell				= $fields[$header{'flowcell'}];
		$lane					= $fields[$header{'lane'}];
		$barcode				= $fields[$header{'barcode'}];
		$barcodeType			= $fields[$header{'barcodeType'}];
		$seqType				= $fields[$header{'seqType'}];
		if (defined($run) && $run eq $run_to_find) { # check if input arg $run_to_find is the same as the run num from the input GAF list.
			$experiment_count++;
			my $truncated_externalSampleID;
			if (length($externalSampleID) > 25) {
				$truncated_externalSampleID = substr($externalSampleID, 0, 22) . '...';
			} else {
				$truncated_externalSampleID = $externalSampleID;
			}
			my $log_message = sprintf($format, $internalSampleID, $truncated_externalSampleID, $sequencingStartDate, $sequencer, $run, $flowcell, $lane, $barcode, $barcodeType, $seqType);
			$logger->info($log_message);
			# Perform checks on input and create hashes of lane/barcode to analyse,
			&_ProcessExperiment($internalSampleID, $externalSampleID, $sequencingStartDate, $sequencer, $run, $flowcell, $lane, $barcode, $barcodeType, $seqType);
		} 
	} else {
		$logger->fatal('Cannot parse CSV line ' . $. . ': ' . $csv->error_input);
		#exit(1);
	}
}

close ($input_fh);

#
# Check if we found any experiments that matched the requested run number.
#
if ($experiment_count == 0) {
	$logger->fatal('Cannot find any (optionally barcoded) samples matching the requested run number ' . $run_to_find);
	exit(1);
}

#
# Check if the sequence type is consistent across the whole run.
#
my $first_sequence_type = $sequence_types[0]; 
foreach my $sequence_type (@sequence_types) {
	unless ($sequence_type eq $first_sequence_type) {
		$logger->fatal('Inconsistent sequence types detected for the (optionally barcoded) samples matching the requested run number ' . $run_to_find . ': ' . 
						$first_sequence_type . ' != ' . $sequence_type);
		exit(1);
	}
}

$logger->info($log_section_break);

&_WriteIlluminaSampleSheet($output, \%experiments);

$logger->info('Finished!');
exit(0);

################################################
########SUBS#######SUBS########SUBS#############
################################################

sub _ProcessExperiment {
	
	my ($internalSampleID, $externalSampleID, $sequencingStartDate, $sequencer, $run, $flowcell, $lanes, $barcode, $barcodeType, $seqType) = @_;
	
	push (@sequence_types, $seqType);
	
	#
	# Check if obligatory fields are not empty.
	#
	foreach my $required_value (@_) {
		if ($required_value eq '') {
			$logger->fatal('Cannot process this GAF sample sheet as one or more of the obligatory fields are empty!');
			exit(1);
		}
	}
	
	#
	# Check if barcodeType is one for which we know how to handle it. 
	#
	unless (exists($gaf_barcode_types{$barcodeType})) {
		$logger->fatal('Detected unknown barcodeType ' . $barcodeType . ' for internalSampleID ' . $internalSampleID . '.');
		$logger->fatal('Known barcodeType values are ' . join(' || ', keys(%gaf_barcode_types)) . '.');
		exit(1);
	}
	
	my @lanes = split(',',$lanes); # Separate multiple lane values like "1,2".
	
	# Need to create:
	# $experiments{$lane}{$barcodeType}{$barcode}{'flowcell'}
	# $experiments{$lane}{$barcodeType}{$barcode}{'internalSampleID'}
	
	foreach my $lane (@lanes) {
		
		if ($lane < 1 || $lane > 8) {
			$logger->fatal('Lane number must be between 1 and 8.');
			exit(1);
		}
		
		if (exists($experiments{$lane})) {
					
			foreach my $previous_barcodeType (keys(%{$experiments{$lane}})) {	
				#
				# Check if barcodeType matches: all samples multiplexed in the same lane must use the same barcoding system.
				#
				if ($previous_barcodeType ne $barcodeType) {
					$logger->fatal('All samples multiplexed in the same lane must use the same barcode type, but multiple barcode types detected for lane ' . $lane);
					exit(1);
				}			
				#
				# Check if barcode not already exists: they must be unique!  
				#
				foreach my $previous_barcode (keys(%{$experiments{$lane}{$previous_barcodeType}})) { 
					if ($previous_barcode eq $barcode) {
						$logger->fatal('Samples multiplexed in the same lane must not share the same barcode, but barcode ' . 
										$barcode . ' is not unique in lane ' . $lane);
						exit(1);
					}
				}
			}
			
			#
			# All OK: Append values for new sample in existing lane. 
			#
			$experiments{$lane}{$barcodeType}{$barcode}{'flowcell'}			= $flowcell;
			$experiments{$lane}{$barcodeType}{$barcode}{'internalSampleID'}	= $internalSampleID;
						
		} else {
			
			#
			# Append values for new sample in new lane.
			#
			$experiments{$lane}{$barcodeType}{$barcode}{'flowcell'}			= $flowcell;
			$experiments{$lane}{$barcodeType}{$barcode}{'internalSampleID'}	= $internalSampleID;		
			
		}
	}
}

# 
# Write new Illumina sample sheet file.
#
sub _WriteIlluminaSampleSheet {
	
	my ($output, $experiments) = @_;
	
	$logger->info('Writing Illumina sample sheet: ' . $output);
	
	# New Illumina example:
	#FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject
	#FC626BWAAXX,4,lane4,Unknown,,'DefaultSample',N,,,FC626BWAAXX
	#
	# Old GAF example:
	#FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator
	#D1P31ACXX,1,D33012_6pM,hg19,CAGATC,1356,N,R1,GAF
	#D1P31ACXX,1,D36656_6pM,hg19,CCGTCC,1357,N,R1,GAF
	
	my $output_fh;
	eval {
		open($output_fh, '>', $output);
	};
	if ($@) {
		$logger->fatal('Cannot create CSV file: ' . $@);
		exit(1);
	}
	
	#
	# Write header.
	#
	print($output_fh 'FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject' . "\n");
	
	#
	# Write records.
	#
	foreach my $lane (sort(keys(%{$experiments}))) {
		foreach my $barcodeType (sort(keys(%{${$experiments}{$lane}}))) {
			
			if ($gaf_barcode_types{$barcodeType}{'IlluminaDemultiplexing'}) {
			
				foreach my $barcode (sort(keys(%{${$experiments}{$lane}{$barcodeType}}))) {
					my $record 	 = ${$experiments}{$lane}{$barcodeType}{$barcode}{'flowcell'} . ',';
					$record		.= $lane . ',';
					# We don't use the sample IDs here: these are used to create FastQ file names
					# Therefore use the lane numbers to create predictable FastQ filenames.
					$record		.= 'lane' . $lane . ',';
					# We don't use the Illumina software for alignment of reads vs. a reference genome.
					# Hence the referene genome is irrelevant; use "Unknown".
					$record		.= 'Unknown' . ',';
					$record		.= $barcode . ',';
					# Use the GAF internalSampleID as "Description".
					$record 	.= 'GAF_sample_ID:' . ${$experiments}{$lane}{$barcodeType}{$barcode}{'internalSampleID'} . ',';
					$record 	.= 'N' . ',';
					$record 	.= 'R1' . ',';
					$record 	.= 'GAF' . ',';
					# We don't use the project names here: these are used to create subfolders for the FastQ files.
					# Therefore use the flowcells to create predictable subdirs.
					$record 	.= ${experiments}{$lane}{$barcodeType}{$barcode}{'flowcell'} . "\n";
					print($output_fh $record);
				}
				
			} else {
				
				my $record;
				my $flowcell;
				my $description = 'GAF_sample_ID:';
				
				foreach my $barcode (sort(keys(%{${$experiments}{$lane}{$barcodeType}}))) {
					$description .= ${$experiments}{$lane}{$barcodeType}{$barcode}{'internalSampleID'} . '-';
					$flowcell     = ${$experiments}{$lane}{$barcodeType}{$barcode}{'flowcell'};
				}
				$description =~ s/-$//;
				
				$record 	 = $flowcell . ',';
				$record		.= $lane . ',';
				# We don't use the sample IDs here: these are used to create FastQ file names
				# Therefore use the lane numbers to create predictable FastQ filenames.
				$record		.= 'lane' . $lane . ',';
				# We don't use the Illumina software for alignment of reads vs. a reference genome.
				# Hence the referene genome is irrelevant; use "Unknown".
				$record		.= 'Unknown' . ',';
				# Only barcodeType = RPI will be demultiplexed by the Illumina software.
				# Samples without barcode or with GAF barcodes must be listed with an empty index value.
				$record		.= ',';
				# Use the GAF internalSampleID as "Description".
				$record 	.= $description . ',';
				$record 	.= 'N' . ',';
				$record 	.= 'R1' . ',';
				$record 	.= 'GAF' . ',';
				# We don't use the project names here: these are used to create subfolders for the FastQ files.
				# Therefore use the flowcells to create predictable subdirs.
				$record 	.= $flowcell . "\n";
					
				print($output_fh $record);
				
			}
		}
	}
	
	close($output_fh);
	
}

sub _Usage {
	
	print <<EOF;
#########################################################################################################
# This script takes a sample sheet in GAF *.csv format and creates one for a specific run in            #
# Illumina *.csv format suitable for conversion of BCL to FastQ files using Illumina's bcl2fastq tool.  #
#########################################################################################################
Usage: CreateIlluminaSmpleSheet.pl [options]
Options:
   -i [file]   Input file in GAF *.csv format.
   -o [file]   Output file in Illumina *.csv format.
   -r [0-9]{4} Run number for which the Illumina sample sheet should be created.
               Must be a four digit number; Use left padding with zeros if the number is less then 1000.
   -l [LEVEL]  Log4perl log level. One of: ALL, TRACE, DEBUG, INFO (default), WARN, ERROR, FATAL or OFF.
#########################################################################################################
EOF

}
