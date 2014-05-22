#!/usr/bin/perl -w

#
# Extract a subset of rows from a *.csv format and save them to a new *.csv file.
# Subset is determined by a query term and the query can be restricted to a specific column.
#

#
# initialise environment
#
use strict;
use Getopt::Long;
use File::Basename;
use Text::CSV;
use Log::Log4perl qw(:easy);

my $ps = '/'; # Path Separator.
my %log_levels = (
	'ALL'	=> $ALL,
	'TRACE'	=> $TRACE,
	'DEBUG'	=> $DEBUG,
	'INFO'	=> $INFO,
	'WARN'	=> $WARN,
	'ERROR'	=> $ERROR,
	'FATAL'	=> $FATAL,
	'OFF'	=> $OFF,
);

my $gaflist;
my $gaflist_subset;
my $query_term;
my $column;
my $log_level;

Getopt::Long::GetOptions(
	'i=s'	=> \$gaflist,			# Input file
	'o=s'	=> \$gaflist_subset,	# Output file
	'q=s'	=> \$query_term,
	'c=s'	=> \$column,			# Restrict queryterm to a column in the GAF list
	'll:s'	=> \$log_level
);

#
# Configure logging.
#
# Provides default if user did not specify log level:
$log_level = (defined($log_level) ? $log_level : 'INFO');
# Reset log level to default if user specified illegal log level.
$log_level = (defined($log_levels{$log_level}) ? $log_levels{$log_level} : $log_levels{'INFO'});
#Log::Log4perl->init('log4perl.properties');
Log::Log4perl->easy_init(
	{ level    => $log_level,
	  file     => "STDOUT",
	  layout   => '%d L:%L %p> %m%n' },
);
my $logger = Log::Log4perl::get_logger();

#
# Check user input.
#
foreach my $input ($gaflist, $gaflist_subset, $query_term, $column) {
	if (defined($input) && length($input) <= 0) {
		$input = undef();
	}
}
#
unless (defined($gaflist) && defined($gaflist_subset) && defined($query_term) && defined($column)) {
	_Usage();
	exit;
}
# Make sure we don't overwrite the input.
foreach my $output ($gaflist_subset) {
	foreach my $input ($gaflist) {
		if (defined($input) && defined($output) && $output eq $input) {
		    $logger->fatal('Output file ' . $output . ' is the same as input file ' . $input . '.');
		    $logger->fatal("\t" . 'Please choose a different file for the output.');
		    exit;
		}
	}
}
# Make sure input files exist.
foreach my $input ($gaflist) {
	unless (-e $input && -f $input && -r $input) {
	    $logger->fatal('Cannot read from input file ' . $input . ': ' . $!);
	    exit;
	}
}

my $headers = _GetHeaders($gaflist);

_GetSubset($gaflist, $headers, $query_term, $column, $gaflist_subset);

$logger->info('Finished!');

#
##
### Subs
##
#

sub _GetHeaders {
	
	my ($gaflist) = @_;
	my @col_headers;
	
	#
	# Parse header.
	#
	my $csv = Text::CSV->new({ binary => 1, eol => $/ });
	open(my $read_csv_fh, '<', $gaflist) or die $!;
	while (my $row = $csv->getline($read_csv_fh)) {
		@col_headers = @{$row};
		#foreach my $col (@col_headers) {
		#	$logger->trace('Found column label/name ' . $col);
		#}
		last if (scalar(@col_headers) > 0);
	}
	close($read_csv_fh);
	
	#
	# Retrieve column indices.
	#
	my %header_indices = (map{ $col_headers[$_] => $_ } (0 .. $#col_headers));
	
	#foreach my $key (keys(%header_indices)) {
	#	$logger->trace('Key: ' . $key . '   Col: ' . $header_indices{$key});
	#}
		
	return(\%header_indices);
	
}

sub _GetSubset {
	
	my ($gaflist, $headers, $query_term, $column, $gaflist_subset) = @_;
	
	unless (defined(${$headers}{$column})) {
		$logger->fatal('Cannot find column with label/name ' . $column);	
	}

	my $csv = Text::CSV->new({ binary => 1, eol => $/ });
	my $csv_subset = Text::CSV->new({ binary => 1, eol => $/ });
	open(my $read_csv_fh, '<', $gaflist) or die $!;
	open(my $write_csv_fh, '>', $gaflist_subset) or die $!;
	
	#
	# Write header to output file.
	#
	my $headerrow = $csv->getline($read_csv_fh);
	$csv_subset->print($write_csv_fh, $headerrow);
	
	#
	# Filter the rows.
	#
	my @values;
	while (my $row = $csv->getline($read_csv_fh)) {
		my @values = @{$row};
		if (defined($values[${$headers}{$column}])) {
			my $value = $values[${$headers}{$column}];
			if ($value eq $query_term) {
				#
				# Store this row.
				#
				$csv_subset->print($write_csv_fh, $row);
			}
		}
	}
	close($read_csv_fh);
	close($write_csv_fh);
	
}

sub _Usage {
        print <<EOF;
###################################################################################################
Extract a subset of rows from a *.csv file and save them to a new *.csv file.
(p.k.a. extract_samples_from_GAF_list.pl)
###################################################################################################
Usage: ./filterCSV.pl
\t --i  [file]    Input file:  *.csv file .
\t --o  [file]    Output file: subset of rows of the input *.csv file.
\t --q  [string]  Query term.
\t --c  [string]  Column name/label from the header line to restrict the search to.
\t --ll [LEVEL]   Log4perl log level.
\t                One of: ALL, TRACE, DEBUG, INFO (default), WARN, ERROR, FATAL or OFF.
###################################################################################################
EOF
 
}