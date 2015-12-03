#!/usr/bin/env perl

#
# ReplaceFastaHeaders.pl
#
# Replaces sequence headers of FASTA files (in various customisable ways).
#

#
# Initialize evironment
#
use strict;
use Getopt::Std;
use Log::Log4perl qw(:easy);
use File::Basename;

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

my $me_myself_and_i = basename($0);

#
# Get options.
#
my %opts;
Getopt::Std::getopts('f:o:l:e:r:i:j:p:q:', \%opts);
my $input				= $opts{'f'};
my $output				= $opts{'o'};
my $log_level			= $opts{'l'};
my $extension			= $opts{'e'};
my $replacements		= $opts{'r'};
my $input_index			= $opts{'i'};
my $replacement_index	= $opts{'j'};
my $input_prefix		= $opts{'p'};
my $replacement_prefix	= $opts{'q'};

#
# Configure logging.
#
# Provides default if user did not specify log level:
$log_level = (defined($log_level) ? $log_level : 'WARN');
# Reset log level to default if user specified illegal log level.
$log_level = (defined($log_levels{$log_level}) ? $log_levels{$log_level} : $log_levels{'WARN'});
#Log::Log4perl->init('log4perl.properties');
Log::Log4perl->easy_init(
	#{ level    => $log_level,
	#  file     => ">>$me_myself_and_i.log",
	#  layout   => '%F{1}-%L-%M: %m%n' },
	{ level    => $log_level,
	  file     => "STDERR",
	  layout   => '%d L:%L %p> %m%n' },
);
my $logger = Log::Log4perl::get_logger();

#
# Start the conversion process.
#
$logger->info("Starting...");

#
# Check user input.
#
# Provides default if user did not specify log level:
$log_level = (defined($log_level) ? $log_level : 'WARN');
# Reset log level to default if user specified illegal log level.
$log_level = (defined($log_levels{$log_level}) ? $log_levels{$log_level} : $log_levels{'WARN'});

# Provide default if user did not specify fasta filename extension:
$extension = (defined($extension) ? $extension : 'fa');

if ($input =~ /^$/ || $replacements =~ /^$/ || $output =~ /^$/) {
	# Indir and outdir cannot be empty.
	_Usage();
}
if ($input eq $output) {
	$logger->fatal("Output dir/file is the same as the input dir/file. Please choose a different one.");
	exit;
}

for my $path ($input, $replacements) {
	unless (-e $path && -r $path) {
		$logger->fatal("Input $input does not exist or is not readable: $!");
		exit;
	}
}

my $find_input_header_key_by_type;
my $find_input_header_key_by_value;
if (defined($input_index)) {
	if (defined($input_prefix)) {
		$logger->fatal("Both -i and -p were specified, but these options are mutually exclusive.");
		exit;
	} else {
		$find_input_header_key_by_type  = 'index';
		$find_input_header_key_by_value = $input_index;
	}
} elsif (defined($input_prefix)) {
	$find_input_header_key_by_type  = 'prefix';
	$find_input_header_key_by_value = $input_prefix;
} else {
	$logger->fatal("Not -i nor -p was specified; must supply one or the other.");
	exit;
}

my $find_replacement_header_key_by_type;
my $find_replacement_header_key_by_value;
if (defined($replacement_index)) {
	if (defined($replacement_prefix)) {
		$logger->fatal("Both -j and -q were specified, but these options are mutually exclusive.");
		exit;
	} else {
		$find_replacement_header_key_by_type  = 'index';
		$find_replacement_header_key_by_value = $replacement_index;
	}
} elsif (defined($replacement_prefix)) {
	$find_replacement_header_key_by_type  = 'prefix';
	$find_replacement_header_key_by_value = $replacement_prefix;
} else {
	$logger->fatal("Not -j nor -q was specified; must supply one or the other.");
	exit;
}

#
# Create hash with replacement headers.
#
my $replacement_headers = _ParseReplacements($replacements, $find_replacement_header_key_by_type, $find_replacement_header_key_by_value);

#
# Check if input is a single file or a directory. 
#
if (-f $input) {
	
	#
	# We've got an input file.
	#
	my $file;
	if ($input =~ m/(.+\/)([^\/]+)$/) {
		$file = $2;
	} else {
		$file = $input;
	}
	
	$logger->info('Parsing ' . $file . "...\n");
	
	_ReplaceFastaHeaders($input, $output, $replacement_headers, $find_replacement_header_key_by_type, $find_replacement_header_key_by_value);
	
	$logger->info('Converted ' . $file);
	
} else {

	#
	# We've got an input directory.
	# Assume the output is also a directory.
	# Append trailing path separators if they was missing.
	#
	my $indir;
	my $outdir;
	unless ($input =~ m/\/$/) {
		$input  = $input .+ '/';
	}
	unless ($output =~ m/\/$/) {
		$output = $output .+ '/';
	}
	#
	# Make sure the input dir is a directory. 
	#
	unless (-d $input) {
	    $logger->fatal("Input $input is not a file nor a directory: $!");
	    exit;
	} else {
		$indir	= $input;
		$outdir	= $output;
	}
	
	#
	# Get all FASTA files from the input dir.
	#
	my $files = _GetFiles($indir, $outdir, $extension);
	
	#
	# Create the output directory if did not exist yet.
	#
	if (-e $outdir && -d $outdir) {
		unless (-w $outdir) {
			$logger->fatal("Cannot write to output directory $outdir. Check for permission errors, read-only file systems, etc.");
			exit;	
		}
	} else {
		$logger->info("Creating output directory $outdir...");
		eval{mkdir($outdir);};
		if ($@) {
			$logger->fatal("Cannot create output directory $outdir: $@");
			exit;	
		}
	}
	
	#
	# Convert FASTA files.
	#
	foreach my $file (@{$files}) {
		
		$logger->info('Parsing ' . $file . "...\n");
		
		my $pathfrom = $indir .+ $file;
		my $pathto   = $outdir .+ $file;
		
		_ReplaceFastaHeaders($pathfrom, $pathto, $replacement_headers, $find_replacement_header_key_by_type, $find_replacement_header_key_by_value);
		
		$logger->info('Converted ' . $file);
	
	}
}

$logger->info('Finished!');

#
##
### Internal subs.
##
#

sub _GetFiles {
	
	my ($indir, $outdir, $extension) = @_;
	my @files;
	
	#
	# Get the relative path to the outdir.
	# Use this to remove it from the list of files/folders that need to be processed
	# in case it's a subfolder of the input directory.
	#
	$outdir =~ m/\/([^\/]+)\/$/;
	my $outdir_rel = $1;
	
	#
	# Get and parse all files from the input dir.
	#
	eval{
		opendir (INDIR, $indir);
		@files = grep { /.+\.$extension/i and not /^\..*/ and not /$outdir_rel/} readdir INDIR;
		closedir INDIR;
	};
	if ($@) {
		$logger->fatal("Cannot read files from input directory $indir: $@");
		exit;
	}

	return(\@files);
}

sub _ParseReplacements {
	
	$logger->debug('_ParseReplacements sub');
	
	my ($replacements, $find_replacement_header_key_by_type, $find_replacement_header_key_by_value) = @_;
	my %replacement_headers;
	
	open(READ,"<$replacements") or die "\tcan't open input file $replacements: $!";
	
	while (my $line = <READ>) {
		
		$line =~ s/[\n\r\f]+//; # Remove line end characters.
		next if $line =~ /^$/;  # Skip empty lines.
		my @values = split(/\s/, $line);
		my $key;
		
		if ($find_replacement_header_key_by_type eq 'index') {
			$key = $values[$find_replacement_header_key_by_value -1];
			if (defined($key)) { 
				$replacement_headers{$key} = $line;
				$logger->debug('Parsed replacement key|header: ' . $key . '|' . $line);
			} else {
				$logger->warn('Failed to parse line: ' . $line);
			}
		} elsif ($find_replacement_header_key_by_type eq 'prefix') {
			my $found_key_by_prefix = 0;
			foreach my $value (@values) {
				if ($value =~ m/^$find_replacement_header_key_by_value(.+)/) {
					$key = $1;
					$replacement_headers{$key} = $line;
					$logger->debug('Parsed replacement key|header: ' . $key . '|' . $line);
					$found_key_by_prefix = 1;
				}
			}
			unless ($found_key_by_prefix == 1) {
				$logger->warn('Failed to parse line: ' . $line);
			}
		} else {
			$logger->fatal('The find_replacement_header_key_by_type variable must contain either index or prefix, but I got: ' . $find_replacement_header_key_by_type);
			exit;
		}
		
	}
	
	return (\%replacement_headers);
	
}

sub _ReplaceFastaHeaders {
		
	$logger->debug('_ReplaceFastaHeaders sub');
	
	my ($pathfrom, $pathto, $replacement_headers, $find_replacement_header_key_by_type, $find_replacement_header_key_by_value) = @_;

	my $header_count = 0;
	
	open(READ,"<$pathfrom") or die "\tcan't open input file $pathfrom: $!";
	open(SAVE,">$pathto") or die "\tcan't open output file $pathto: $!";
	while (my $line = <READ>) {
		
		$line =~ s/[\n\r\f]+//; # Remove line end characters.
		next if $line =~ /^$/;  # Skip empty lines.
		my $new_line;
			
		if ($line =~ /^>/) {
			
			#
			# It's a header line.
			#				
			$header_count++;
			my @values = split(/\s/, $line);
			my $key;
			
			if ($find_input_header_key_by_type eq 'index') {
				$key = $values[$find_input_header_key_by_value -1];
			} elsif ($find_input_header_key_by_type eq 'prefix') {
				my $found_key_by_prefix = 0;
				foreach my $value (@values) {
					if ($value =~ m/^$find_input_header_key_by_value(.+)/) {
						$key = $1;
						$found_key_by_prefix = 1;
					}
				}
				unless ($found_key_by_prefix == 1) {
					$logger->warn('Failed to parse line: ' . $line);
				}
			} else {
				$logger->fatal('The find_input_header_key_by_type variable must contain either index or prefix, but I got: ' . $find_input_header_key_by_type);
				exit;
			}
			
			if (defined($key) && defined(${$replacement_headers}{$key})) { 
				$logger->debug('Parsed input key|header: ' . $key . '|' . $line);
				$new_line = ${$replacement_headers}{$key};
				$logger->debug('Found replacement key|header: ' . $key . '|' . $new_line);
			} else {
				$logger->warn('Failed to find replacement header for line: ' . $line);
				$new_line = $line;
			}
			
		} else {
			
			#
			# It must be a sequence line.
			#
			$new_line = $line;
			
		}
		
		# Save (modified) line.
		print SAVE $new_line . "\n" or die "\tcan't save output to file $pathto: $!";
		
	}
	
	close(READ);
	close(SAVE);
	
}

sub _Usage {

	print "\n";
	print "$me_myself_and_i - Replaces sequence headers of FASTA files.\n";
	print "\n";
	print "Usage:\n";
	print "\n";
	print "   $me_myself_and_i options\n";
	print "\n";
	print "Available options are:\n";
	print "\n";
    print "   -f [dir/file]    Input can be a single FASTA file or a directory containing FASTA files.\n";
    print "   -e [ext]         File name extension for the FASTA files in case the input is a directory. (default = fa)\n";
    print "   -r [file]        Input file with replacement headers.\n";
    print "   -o [dir/file]    Output file or directory where the result(s) will be saved.\n";
    print "   -i [index]       Position of the identifier in the space seperated header of the FastA file, which will be used to find the replacment header.\n";
    print "   -j [index]       Position of the identifier in the space seperated header of the file with replacement headers, which will be used to find the replacment header.\n";
    print "   -p [prefix]      Prefix to find the identifier in the space seperated header of the FastA file, which will be used to find the replacment header.\n";
    print "   -q [prefix]      Prefix to find the identifier in the space seperated header of the file with replacement headers, which will be used to find the replacment header.\n";
    print "   -l [LEVEL]       Log4perl log level. One of: ALL, TRACE, DEBUG, INFO (default), WARN, ERROR, FATAL or OFF.\n";
    print "\n";
    print "Examples:\n";
    print "\n";
    print "    Given this example header in the FastA file:\n";
    print "        >chr1  AC:CM000663.2  gi:568336023  LN:248956422  rl:Chromosome  M5:6aef897c3d6ff0c78aff06ac189178dd  AS:GRCh38\n";
    print "    and this example header in the file with replacement headers :\n";
    print "        >NC_000001.11(1)  AC:CM000663.2  LN:248956422  rl:Chromosome  M5:6aef897c3d6ff0c78aff06ac189178dd  AS:GRCh38\n";
    print "    the following will replace the header based on prefixes:\n";
    print "            $me_myself_and_i -i FastA.fa -r new_headers.txt -p 'AC:' -q 'AC:'\n";
    print "    the following will replace the header based on prefixes:\n";
    print "            $me_myself_and_i -i FastA.fa -r new_headers.txt -i 2 -j 2\n";
	print "\n";
	exit;

}
