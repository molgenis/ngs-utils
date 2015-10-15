#!/usr/bin/env perl

#
# This script takes as input:
# * a sample sheet in GAF *.csv format and 
# * an input dir with Genome Studio files - potentially containing multiple samples.
# to create Genome Studio files per sample.
#
# =====================================================
# $Id: CreateFinalReportPerSample.pl 1055 2013-04-22 17:53:50Z pneerincx $
# $URL: http://www.bbmriwiki.nl/svn/ngs_scripts/trunk/scripts_data_archiving/CreateFinalReportPerSample.pl $
# $LastChangedDate: 2013-04-22 19:53:50 +0200 (Mon, 22 Apr 2013) $
# $LastChangedRevision: 1055 $
# $LastChangedBy: pneerincx $
# =====================================================
# 


use strict;
use warnings;
use Getopt::Std;
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

my @gaf_column_names = ('internalSampleID','externalSampleID','sequencingStartDate','sequencer','run','flowcell','lane','barcode','barcodeType','seqType','arrayFile','arrayID');
my %gaf_barcode_types = ('None' => {
							'IlluminaDemultiplexing'	=> 0},
						 'GAF' => {
							'IlluminaDemultiplexing'	=> 0},
						 'RPI' => {
							'IlluminaDemultiplexing'	=> 1},
						 'AGI' => {
							'IlluminaDemultiplexing'	=> 1},
						 'MON' => {
							'IlluminaDemultiplexing'	=> 1},
						);
my $log_section_break = '##################################################################################################################################';

#
# Get options.
#
my %opts;
Getopt::Std::getopts('i:o:s:r:l:', \%opts);

my $inputdir	= $opts{'i'};
my $outputdir	= $opts{'o'};
my $run_to_find	= $opts{'r'};
my $file		= $opts{'s'};
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

# mandatory args
_Usage() and exit(1) unless $run_to_find;
_Usage() and exit(1) unless $file;
_Usage() and exit(1) unless $inputdir;
_Usage() and exit(1) unless $outputdir;

chomp $run_to_find;
chomp $file;
chomp $inputdir;
chomp $outputdir;

$logger->info('Checking ' . $file);

#
# Retrieve column headers from CSV file.
#
my $csv = Text::CSV->new();
open (CSV, "<", $file) or die $!;
my @headers;

while (<CSV>) {
	next if ($. != 1);
	if ($csv->parse($_)) {
		@headers = $csv->fields();
	} else {
		my $err = $csv->error_input;
		$logger->fatal('Failed to parse line: ' . $err);
		exit(1);
	}
}

close CSV;

#
# Retrieve column indices.
#
my %headers_to_indexes = (map { $headers[$_] => $_ } (0 .. $#headers));

$logger->info($log_section_break);
$logger->info('Found the following samples:');
my $format = '%16s  %25.25s  %19s  %9s  %4s  %10s  %4s  %7s  %11s  %7s  %25.25s  %2s';
my $log_header = sprintf($format, @gaf_column_names);
$logger->info($log_header);

#
# Retrieve values.
#
my %mappings;
my ($internalSampleID, $externalSampleID, $sequencingStartDate, $sequencer, $run, $flowcell, $lane, $barcode, $barcodeType, $seqType, $arrayFile, $arrayID);
$csv = Text::CSV->new();
open (CSV, "<", $file) or die $!;
while (<CSV>) {
	next if ($. == 1);
	if ($csv->parse($_)) {
		my @columns = $csv->fields();
		$internalSampleID = $columns[$headers_to_indexes{"internalSampleID"}];
		$externalSampleID = $columns[$headers_to_indexes{"externalSampleID"}];
		$sequencingStartDate = $columns[$headers_to_indexes{"sequencingStartDate"}];
		$sequencer = $columns[$headers_to_indexes{"sequencer"}];
		$run = $columns[$headers_to_indexes{"run"}];
		$flowcell = $columns[$headers_to_indexes{"flowcell"}];
		$lane = $columns[$headers_to_indexes{"lane"}];
		$barcode = $columns[$headers_to_indexes{"barcode"}];
		$barcodeType = $columns[$headers_to_indexes{"barcodeType"}];
		$seqType = $columns[$headers_to_indexes{"seqType"}];
		$arrayFile = $columns[$headers_to_indexes{"arrayFile"}];
		# Retrieve the arrayfilename from the complete path to the file (on a windows SAMBA share).
		if ($arrayFile =~ m/([^\\\/]+)$/) {
			$arrayFile = $1;
		}
		$arrayID = $columns[$headers_to_indexes{"arrayID"}];
		if ($run eq "$run_to_find") {
			my $truncated_externalSampleID = _TruncateLongValue($externalSampleID, 25);
			my $truncated_arrayFile = _TruncateLongValue($arrayFile, 25);
			my $log_message = sprintf($format, $internalSampleID, $truncated_externalSampleID, $sequencingStartDate, $sequencer, $run, $flowcell, $lane, $barcode, $barcodeType, $seqType, $truncated_arrayFile, $arrayID);
			$logger->info($log_message);
			
			# Perform checks on input and create hashes of sample => arrayFile/arrayID combo.
			if ($externalSampleID eq '' || $arrayFile eq '' || $arrayID eq ''){
				# No files can be created for this sample.
				$logger->warn('Cannot create *FinalReport.txt: arrayFile and/or arrayID missing for sample ' . $internalSampleID . ' (' . $externalSampleID . ')');
			}else{
				my $concat = $arrayFile . '::' . $arrayID;
				$mappings {$externalSampleID} = $concat;
			}
		}
	}
}

#
# Check for non unique sample <-> arrayFile_arrayID combinations.
#
my %reverse;
while (my ($key, $value) = each %mappings) {
	push @{$reverse{$value}}, $key;
}

while (my ($key, $value) = each %reverse) {
	next unless @$value > 1;
	$logger->fatal('Found ' . scalar(@$value) . ' samples with the same arrayFile and arrayID!');
	$logger->fatal("\t" . 'Samples ' . @$value . ' have the same arrayFile::arrayID ' . $key . ' combination!');
	exit(1);
}

#
# Generate per sample finalreports.
#
$logger->info('Generating sample finalreports...');

foreach my $key (keys %mappings){
	$logger->debug('KEY: ' . $key . ' ' . $mappings{$key});
	my @array = split("::", $mappings{$key});
	my $file = $array[0];
	my $id = $array[1];
	
	#Check if finalreport exists
	my $to_check = $inputdir . '/' . $file . '_FinalReport.txt';
	if (-e $to_check && -r $to_check) {
		#File exists, do nothing;
		$logger->trace('Found FinalReport: ' . $to_check);
	} else {
		$logger->fatal('FinalReport ' . $to_check .' does not exist or is not readable.');
		exit(1);
	}

	# MD5 the parent FinalReport.
	my $checksum_file = $to_check . '.md5';
	if (-e $checksum_file) {
		# Checksum exists; don't overwrite!
	} else {
		my $result = `md5sum $to_check > $checksum_file`;
		if ($result ne '' || $? > 0) {
			$logger->fatal('Creating md5 checksum for ' . $to_check . ' failed: '. $result);
			exit(1);
		} else {
			# md5summing was a success.
			`chmod 0640 $checksum_file`;
			$logger->trace('Created md5 checksum file: ' . $checksum_file . '.');
		}
	}
	
	# Retrieve tenth line from file (contains header)
	my $header = `head -10 $inputdir/$file\_FinalReport.txt | tail -1`;
	$header =~ s/\R//g; #\R is the general linebreak!

	# Retrieve column indices.
	my @headers = split("\t", $header);

	my %head_to_index = (map { $headers[$_] => $_ } (0 .. $#headers));

	my $snpname = $head_to_index{"SNP Name"};
	my $sampleindex = $head_to_index{"Sample Index"};
	my $all1_fwd = $head_to_index{"Allele1 - Forward"};
	my $all2_fwd = $head_to_index{"Allele2 - Forward"};
	my $gcscore = $head_to_index{"GC Score"};
	my $chr = $head_to_index{"Chr"};
	my $pos = $head_to_index{"Position"};
	
	# Create per sample finalreport.
	my $sample_final_report_file_name = $key . '_FinalReport.txt';
	my $sample_final_report_file_path = $outputdir . '/' . $sample_final_report_file_name;
	`touch $sample_final_report_file_path`;
	`chmod 0640 $sample_final_report_file_path`;
	`head -10 $inputdir/$file\_FinalReport.txt > $sample_final_report_file_path`;
	`perl -pi -e 's/Sample Index/Sample Name/' $sample_final_report_file_path`;

	my $retrieveinfo = "awk '\$" . ($sampleindex+1) . " == \"" . $id . "\" \{ if \(\$" . ($all1_fwd+1) . " \!\= \"\-\" \|\| \$" . ($all2_fwd+1) . " \!\= \"\-\"\) print \$" . ($snpname+1) . ",\"" . $key . "\",\$" . ($all1_fwd+1) . ",\$" . ($all2_fwd+1) . ",\$" . ($gcscore+1) . ",\$" . ($chr+1) . ",\$" . ($pos+1) . "\}\;' OFS=\"\\t\" " . "$inputdir\/$file\_FinalReport.txt" . " >> " . "$sample_final_report_file_path";

	# Execute above defined awk statement.
	`$retrieveinfo`;

	# MD5 the new FinalReport.
	chdir($outputdir);
	$checksum_file = $sample_final_report_file_path . '.md5';
	my $result = `md5sum $sample_final_report_file_name > $checksum_file`;
	if ($result ne '' || $? > 0) {
		$logger->fatal('Creating md5 checksum for ' . $sample_final_report_file_path . ' failed: '. $result);
		exit(1);
	} else {
		# md5summing was a success.
		`chmod 0640 $checksum_file`;
		$logger->trace('Created md5 checksum file: ' . $checksum_file . '.');
	}

	$logger->info('Created ' . $sample_final_report_file_path);
}

$logger->info('Finished!');
exit(0);

################################################
########SUBS#######SUBS########SUBS#############
################################################

sub _TruncateLongValue {

	my ($value, $max_length) = @_;
	my $truncated_value;
	
	if (length($value) > $max_length) {
		$truncated_value = substr($value, 0, $max_length - 3) . '...';
	} else {
		$truncated_value = $value;
	}
	
	return($truncated_value);
	
}

sub _Usage {
	
	print <<EOF;
#########################################################################################################
# This script uses a final report generated with GenomeStudio to create a finalreport per sample.       #
#########################################################################################################
Usage: ./CreateFinalReportPerSample.pl
   -i dir       Input directory containing *FinalReport.txt files with multiple samples.
   -o dir       Output directory to save per sample *FinalReport.txts.
   -s file.csv  Sample sheet in *.csv format.
   -r [0-9]{4}  Run number from the sample sheet for which to create per sample *FinalReport.txts.
                Must be a four digit number; Use left padding with zeros if the number is less then 1000.
   -l [LEVEL]   Log4perl log level. One of: ALL, TRACE, DEBUG, INFO (default), WARN, ERROR, FATAL or OFF.
#########################################################################################################
EOF
 
}
