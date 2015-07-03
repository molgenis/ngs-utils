#!/usr/bin/env perl -w

#
# Script to create a BED file with genomic regions for HGNC GeneSymbols
# using annotation as provided by Ensembl and quertied remotely via the 
# Ensembl API. 
#
# Make sure Ensembl Perl API and BioPerl are
#  * installed and
#  * relevant paths are added to your $PERL5LIB environment variable if not installed in a default location.
#

#
# initialise environment
#
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use Log::Log4perl qw(:easy);
use DBI;
use DBD::mysql;
use Bio::EnsEMBL::Registry;

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

my %targets = (
	'all' => 1, # any kind of transcript.
	'CDS' => 2, # protein Coding DNA Sequence; hence excluding UTR.
	'CT'  => 3, # protein Coding Transcript; hence including UTR.
);

#
# Get options.
#
my $input;
my $output;
my $log_level = 'INFO';					# Default log level.


Getopt::Long::GetOptions (
	'i=s'		=> \$input,
	'o=s'		=> \$output,
	'l=s'		=> \$log_level,
);

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
##
### Main
##
#

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
unless (-e $input && -r $input) {
	$logger->fatal($input . ': File does not exist or is not readable.');
	exit;
}

my ($resolved_biotypes_counter) = _Parse($input, $output);

$logger->info('===================================================================');
$logger->info('Found biotype for ' . $resolved_biotypes_counter . ' transcripts.');
$logger->info('===================================================================');
$logger->info('Finished!');

#
##
### Subs
##
#

#
# Usage
#
sub _Usage {
	
	print STDERR "\n"
	  . 'GeneSymbol2BED.pl:' . "\n\n"
	  . '   Script to create a BED file with genomic regions for HGNC GeneSymbols' . "\n"
	  . '   using annotation as provided by Ensembl and queried remotely via the Ensembl API.' . "\n\n"
	  . 'Usage:' . "\n\n"
	  . '   GeneSymbol2BED.pl options' . "\n\n"
	  . 'Available options are:' . "\n\n"
	  . '   -i [file]    Input file in tab delimited format. Must contain at least 2 columns: ' . "\n"
	  . '   -o [file]    Output file in tab delimited BED format.' . "\n"
	  . '   -l [LEVEL]   Log4perl log level. One of: ALL, TRACE, DEBUG, INFO (default), WARN, ERROR, FATAL or OFF.' . "\n"
	  . "\n";
	exit;
	
}

#
# Extract info from Ensembl features and return that as tab delimited string.
#

sub _DumpForDebugging {
	my $data_to_dump = @_;
	use Data::Dumper;
	print Data::Dumper->Dump($data_to_dump);
	exit 1;
}

sub _GetFeatureDetails {
	
	my $feature		= shift;
	my $stable_id	= $feature->stable_id();
	my $seq_name	= $feature->slice->seq_region_name();
	my $start		= $feature->start();
	my $end			= $feature->end();
	my $strand		= $feature->strand();
	
	return($stable_id, $seq_name, $start, $end, $strand);
	
}

sub _Parse {
	
	my ($file_path_in, $file_path_out) = @_;
	
	$logger->info('Parsing Ensembl transcript IDs from ' . $file_path_in . '...');
	$logger->info('and storing their biotypes in ' . $file_path_out . '...');
	
	my $file_path_in_fh;
	my $file_path_out_fh;
	my $line_counter = 0;
	my $resolved_biotypes_counter = 0;
	
	#
	# Create filehandles.
	#
	eval {
		open($file_path_in_fh, "<$file_path_in");
	};
	if ($@) {
		$logger->fatal('Cannot read input file: ' . $@);
		exit;
	}
	eval {
		open($file_path_out_fh, ">$file_path_out");
	};
	if ($@) {
		$logger->fatal('Cannot write output file: ' . $@);
		exit;
	}
	
	#
	# Connect to public Ensembl database.
	#
	my $registry = 'Bio::EnsEMBL::Registry';
	$registry->load_registry_from_db(
		-host => 'ensembldb.ensembl.org',    # alternatively 'useastdb.ensembl.org'
		-user => 'anonymous'
	);
	
	my $transcript_adaptor = $registry->get_adaptor('Human', 'Core', 'Transcript');
	
	LINE: while (my $line = <$file_path_in_fh>) {
		
		$line_counter++;
		$logger->debug('Parsing record/line ' . $line_counter . '.');
		my $transcript_id;
		
		#
		# Skip empty lines and chomp newline, return and formfeed chars.
		#
		$line =~ s/[\r\n\f]+//g;
		next LINE if ($line =~ m/^$/);
		if ($line =~ m/^(ENST[0-9]+)$/) {
			$transcript_id = $1;
		} else {
			$logger->warn($transcript_id . ' is not a valid Ensembl transcript ID.');
			next LINE;
		}
		
		#
		# Get the transcript(s) from Ensembl core database.
		#
		my $transcript = $transcript_adaptor->fetch_by_stable_id($transcript_id);
		
		#
		# Retrieve transcript info.
		#
		my $id = $transcript->stable_id();
		if ($id eq $transcript_id) {
			$logger->debug('Found Ensembl ID ' . $transcript_id . '.');
		} else {
			$logger->fatal('Ensembl ID ' . $transcript_id . ' does not match found ' . $id . '.');
			exit 1;
		}
		
		my $biotype = $transcript->biotype();
		my $status = $transcript->status();
		#my ($tr_id, $tr_seq, $tr_start, $tr_end, $tr_strand) = _GetFeatureDetails($transcript);
		#$logger->trace("\t" . join("\t", $tr_id, $tr_seq, $tr_start, $tr_end, $tr_strand));
		
		#
		# Store biotype.
		#
		print $file_path_out_fh join("\t", $transcript_id, $status , $biotype) . "\n";
		$resolved_biotypes_counter++;
		
	}
	
	close($file_path_in_fh);
	close($file_path_out_fh);
	
	return($resolved_biotypes_counter);
	
}
