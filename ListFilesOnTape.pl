#!/usr/bin/perl -w
#
# This script parse Tivoli dsmc output to list files available on backup tapes.
#

#
# initialise environment
#
use strict;
use Getopt::Long;
use Log::Log4perl qw(:easy);
use Vcf;
#use Data::Dumper;

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
my %match_logic_types = (
	'literal'   => 1,
	'regex' 	=> 1,
);

#
# Get options.
#
my $input;
my $output;
my $files;
my $log_level = 'INFO';			# The default log level.

Getopt::Long::GetOptions (
	'i=s'	=> \$input,
	'o=s'	=> \$output,
	'll=s'		=> \$log_level);

#
# Configure logging.
#
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
unless (defined($input) && defined($output)) {
	_Usage();
	exit;
}
if ($input =~ /^$/ || $output =~ /^$/) {
	_Usage();
	exit;
}
if ($input eq $output) {
	$logger->fatal('Output file is the same as the input file. Select a different file for the output.');
	exit;
}

#
# Check input file.
#
unless (-e $input && -f $input && -r $input) {
	$logger->fatal('Cannot read from input file ' . $input . ': ' . $!);
	exit;
}

$files=_ParseFileList($input, $output);

$logger->info('Processed ' . $files . ' files.');
$logger->info('Finished!');

#
##
### Internal subs.
##
#

sub _ParseFileList {
	
	my ($input_file_path, $output_file_path) = @_;
	
	$logger->info('Parsing ' . $input_file_path . '...');
	
	my $input_file_fh;
	my $output_file_fh;
	my $files_parsed = 0;
	
	eval {
		open($input_file_fh, "<$input_file_path");
	};
	if ($@) {
		$logger->fatal('Cannot read from file: ' . $@);
		exit;
	}
	eval {
		open($output_file_fh, ">$output_file_path");
	};
	if ($@) {
		$logger->fatal('Cannot write to file: ' . $@);
		exit;
	}
	
	local $/ = 'Compressed:  NO		Encryption Type:        None';
	LINE: while (my $record = <$input_file_fh>) {
		
		#$record =~ s/[\r\n]+//g;
		if ($record =~ m/(\/target\/gpfs2\/gcc\/[^\n]+)\n.*Modified:\s([0-9]{2})\/([0-9]{2})\/([0-9]{2})\s+([0-9]{2}:[0-9]{2}):[0-9]{2}/) {
			my $file_on_tape = $1 . "\t" . $4 . '/' . $2 . '/'. $3 . '   ' . $5;
			$logger->debug('Found file ' . $file_on_tape);
			$files_parsed++;
			print $output_file_fh $file_on_tape . "\n";
		}
	}
	
	close($input_file_fh);
	close($output_file_fh);
	$logger->info('Finished parsing.');
	
	return ($files_parsed);
	
}

sub _Usage {

	print STDERR "\n" .
	  'ListFilesOnTape.pl:' . "\n" .
	  '   . ' . "\n\n" .
	  'Usage:' . "\n" .
	  '   ListFilesOnTape.pl [options].' . "\n\n" .
	  'Available options are:' . "\n" .
	  '   -i    [file]   Input file.' . "\n" .
	  '   -o    [file]   Output file.' . "\n" .
	  '   -ll   [LEVEL]  Log4perl log level. One of: ALL, TRACE, DEBUG, INFO (default), WARN, ERROR, FATAL or OFF.' . "\n" .
	  "\n";
	exit;

}
