#!/usr/bin/perl -w
#
# This script fills the GenoType Count (GTC) sub-fields of the INFO fields in a VCF file.
# based on a list of accession numbers / IDs
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
my $filter_samples;				# Filter on Sample IDs!
my $filter_FILTER;				# Filter on VCF FILTER column values.
my %vcf_filters;
my $drop = 0;					# The default is to keep the sample IDs supplied in the filter list.
my $match_logic = 'literal';
my $strip_info = 0;				# Whether to strip all INFO subfield except for AN, AC and the ones from $info_sub_fields_to_keep.
								# The default is to keep all INFO subfields.
my $keep_info;					# Keys of INFO subfields to keep in addition to AN and AC.
my $strip_samples = 0;			# The default is to keep the info per sample.
								# When enabled we remove sample details and only keep summary stats in the INFO field.
my $log_level = 'INFO';			# The default log level.

Getopt::Long::GetOptions (
	'vcfi=s'	=> \$input,
	'vcfo=s'	=> \$output,
	'fs=s'		=> \$filter_samples,
	'fv=s'		=> \$filter_FILTER,
	'd!'		=> \$drop,
	'ml=s'		=> \$match_logic,
	'si!'		=> \$strip_info,
	'ki=s'		=> \$keep_info,
	'ss!'		=> \$strip_samples,
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

	#{ level    => $log_level,
	#  file     => ">>vcf-fill-gtc.log",
	#  layout   => '%F{1}-%L-%M: %m%n' },
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
# Check input VCF file.
#
unless (-e $input && -f $input && -r $input) {
	$logger->fatal('Cannot read from input file ' . $input . ': ' . $!);
	exit;
}

#
# Check filter settings.
# 
#  * Determine whether we need to filter for certain samples.
#  * Determine whetehr we should use literal or regex matching for the names/IDs of the samples.
#
if (defined($filter_samples)) {
	unless (-e $filter_samples && -f $filter_samples && -r $filter_samples) {
		$logger->fatal('Cannot read from filter file ' . $filter_samples . ': ' . $!);
		exit;
	}
	$match_logic = (defined($match_logic) ? $match_logic : 'literal');
	unless (defined($match_logic_types{$match_logic})) {
		$logger->fatal('Unkown logic ' . $match_logic . ' specified for sample filtering.');
		exit;
	}
	if ($match_logic eq 'regex') {
		$logger->fatal('Matching logic ' . $match_logic . ' for sample filtering not implemented (yet).');
		exit;
	}
}
if (defined($filter_FILTER)) {
	my @filters = split(',', $filter_FILTER);
	foreach my $filter (@filters) {
		$vcf_filters{$filter} = 1;
	}
}

#
# Create a list of sample IDs to search the VCF file for.
#
my $samples = [];
my $samples_to_filter;
if (defined($filter_samples)) {
	$samples = _CreateLookupList($filter_samples, $match_logic);
	$samples_to_filter = scalar(@{$samples});
	$logger->info('Number of samples to filter for: ' . $samples_to_filter);
}

#
# Create a lis of INFO sub field keys to keep when stripping INFO fields due to -si or -fs.
#
my %info_sub_fields_to_keep;
if (defined($keep_info) && $keep_info ne '') {
	%info_sub_fields_to_keep = map({ $_ => 1 } split(/,/, $keep_info));
}

#
# Create file handle for the new VCF file
#
my $output_fh;
eval {
	open($output_fh, ">$output");
};
if ($@) {
	$logger->fatal('Cannot write VCF file: ' . $@);
	exit;
}

#
# Parse VCF records
#
my ($stats) = _ParseVcf($input, $output_fh, $samples, $drop, $strip_samples, \%vcf_filters, $strip_info, \%info_sub_fields_to_keep);

#
# Clean up.
#
close($output_fh);

$logger->info('New VCF file contains ' . ${$stats}{'variants_stored'} . '/' . ${$stats}{'variants_parsed'} . ' variants from ' 
										. ${$stats}{'samples_stored'} . '/' . ${$stats}{'samples_parsed'} . ' samples.');
$logger->info('Finished!');

#
##
### Internal subs.
##
#

sub _CreateLookupList {
	
	my ($file_path, $match_logic) = @_;
	
	$logger->info('Parsing ' . $file_path . '...');
	
	my @wanted;
	my $file_path_fh;
	
	eval {
		open($file_path_fh, "<$file_path");
	};
	if ($@) {
		$logger->fatal('Cannot read filter file: ' . $@);
		exit;
	}
	
	LINE: while (my $line = <$file_path_fh>) {
	
		$line =~ s/[\r\n]+//g;
		$logger->debug('Found ID ' . $line);
		push(@wanted, $line);
	
	}
	
	close($file_path_fh);
	$logger->info('Created ID lookup list.');
	
	return (\@wanted);
	
}

sub _ParseVcf {
	
	my ($path_from, $path_to_fh, $samples, $drop, $strip_samples, $vcf_filters, $strip_info, $info_sub_fields_to_keep) = @_;
	my %stats;
	
	my $vcf = Vcf->new(file=>$path_from, version => '4.1');
	$vcf->parse_header();
	$logger->trace('Parsed header.');
	
	#use Data::Dumper;
	#$logger->fatal('######### Header: ##########');
	#$logger->fatal(Data::Dumper->Dump([$vcf]));
	#$logger->fatal('############################');
	
	#
	##
	### Create new header.
	##
	#
	
	my %sample_id_values;	
	if ($strip_samples) {
		
		#
		# Create a list of sample IDs by parsing the sample column headers.
		#
		# When sample details must be removed from the VCF file,
		# we use the list of sample IDs created here to check later on
		# whether any sample ID related values are present in any of 
		# the remaining meta-data header lines.
		#
		my (@sample_ids) = $vcf->get_samples();
		foreach my $sample_id (@sample_ids) {
			$logger->trace('Found sample ID: ' . $sample_id . '.');
			my @values = split("[\._\-]", $sample_id);
			foreach my $value (@values) {
				$sample_id_values{$value} = 'drop';
			}
		}
		$logger->trace('Finished parsing sample IDs.');
		
		#
		# When keeping only the aggregate data and removing the sample specific data:
		# Remove the FORMAT field's meta-data header line.
		# Note:  the FORMAT field's column label will be removed later on when formatting the new header.
		#
		$vcf->remove_header_line(key=>'FORMAT');
		$logger->trace('Removed FORMAT header line.');
		
		#
		# Remove any header line starting with the word "source" and containing "vcf".
		#
		while (my ($hline_key,$hline_hashes) = each %{$$vcf{header}}) {
			if ($hline_key =~ m/source/i) {
				$logger->trace('Header key: '. $hline_key);
				foreach my $hline_hash (@{$hline_hashes}) {
					foreach my $key (keys(%{$hline_hash})) {
						if (${$hline_hash}{$key} =~ m/vcf/i) {
							$vcf->remove_header_line(key=>$hline_key);
							$logger->info('Removed ' . $hline_key . ' header line, because it contains sample details.');
						}
					}
				}
			}
		}
		
		#
		# Remove INFO:NS and INFO:SF subfields as those also provide details of the "source".
		#
		my $header_lines = $vcf->get_header_line(key=>'INFO');
		foreach my $header_line_id (keys(%{${$header_lines}[0]})) {
			$logger->trace('INFO header line ID: ' . $header_line_id);
			if ($header_line_id eq 'NS' || $header_line_id eq 'SF') {
				$vcf->remove_header_line(key=>'INFO', ID=>$header_line_id);
				$logger->info('Removed INFO:' . $header_line_id . ' header lines.');
			}
		}
		
	}
	
	#
	# Determine which samples to process.
	#
	#  When removing specific samples from a VCF: force ${strip_info}; 
	#    We force dropping all INFO sub fields except INFO:AC and INFO:AN 
	#    as those two are automagically updated when samples are removed/added, 
	#    but the others are not.
	#
	my @samples_parsed = $vcf->get_samples();
	if (scalar(@{$samples}) > 0 && $drop == 0) {
		$logger->trace('Keeping samples from filter list.');
		$vcf->set_samples(include=>$samples);
		$strip_info = 1;
	} elsif (scalar(@{$samples}) && $drop == 1 ) {
		$logger->trace('Dropping samples from filter list.');
		$vcf->set_samples(exclude=>$samples);
		$strip_info = 1;
	} else {
		$logger->trace('Keeping all samples from input VCF.');
		$vcf->set_samples(exclude=>[]);
	}
	#my @samples_included = $vcf->get_samples();
	my @samples_included;
	for (my $column_offset = 9; $column_offset < scalar(@{$$vcf{columns}}); $column_offset++) {
		if ($$vcf{samples_to_parse}[$column_offset] == 1) {
			push (@samples_included, $$vcf{columns}[$column_offset]);
		}
	}
	
	#
	# Count how many samples are in the original input VCF and how many are left.
	#
	$stats{'samples_parsed'} = scalar(@samples_parsed);
	$stats{'samples_stored'} = scalar(@samples_included);
	
	#
	# When removing specific samples from a VCF:
	#  * By default we drop all INFO sub fields except INFO:AC and INFO:AN 
	#    as those two are automagically updated when samples are removed/added, 
	#    but the others are not.
	#  * You can manually override this behanviour by supplying a list of INFO 
	#    field keys that should be retained in addition to INFO:AC and INFO:AN.
	#
	if ($strip_info > 0) {
		my $header_lines = $vcf->get_header_line(key=>'INFO');
		foreach my $header_line_id (keys(%{${$header_lines}[0]})) {
			$logger->trace('INFO header line ID: ' . $header_line_id);
			unless ($header_line_id eq 'AN' || $header_line_id eq 'AC' || defined(${$info_sub_fields_to_keep}{$header_line_id})) {
				$vcf->remove_header_line(key=>'INFO', ID=>$header_line_id);
				$logger->info('Removed INFO header lines for sub field key ' . $header_line_id  . ', because it\'s value is not automagically updated and you are removing or adding samples.');
			}
		}
	} else {
		$logger->info('Keeping (remaining) INFO header lines.');
	}
	
	#
	# Add header line for the INFO:GTC field
	#
	$vcf->add_header_line({key=>'INFO',ID=>'GTC',Number=>'G',Type=>'Integer',
		Description=>'GenoType Counts listed in the same order as the ALT alleles: 0/0,0/1,1/1,0/2,1/2,2/2,0/3,1/3,2/3,3/3,etc. ' . 
					'Phasing is ignored; hence 0/1, 1/0, 0|1 and 1|0 are all counted as 0/1. ' . 
					'Incomplete gentotypes (./., ./0, ./1, ./2, etc.) are completely discarded for calculating GTC.'});
	
	my $header_new = $vcf->format_header(\@samples_included);
	if ($strip_samples) {
		#
		# When keeping only the aggregate data and removing the sample specific data:
		# Keep column labels only for the first 8 field.
		#
		$header_new =~ s/(#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO)\t.*/$1/;
		$logger->info('Removed FORMAT and sample column labels from column header line.');
	}
	$logger->debug('VCF header new: ' . $header_new);
	
	#
	# Sanity check when individual level details should be removed: 
	# find remants of the sample IDs in the new header and issue a warning if found.
	#
	if ($strip_samples) {
		foreach my $sample_id_value (keys(%sample_id_values)) {
			if (index($header_new, $sample_id_value) != -1) {
				$logger->warn('Found potential remnant of a sample ID value in the new VCF header: ' . $sample_id_value . '.');
			} 
		}
	}
	
	#
	# Write new header to new VCF file.
	#
	print($output_fh $header_new);
	
	#
	##
	### Create new records.
	##
	#
	
	$logger->info('Parsing VCF records...');
	$stats{'variants_parsed'} = 0;
	$stats{'variants_stored'} = 0;
	
	#
	# Parse variant records.
	#
	VCF_RECORD:while (my $record=$vcf->next_data_hash()) {
		
		my $seq_id   = ${$record}{'CHROM'};
		my $pos      = ${$record}{'POS'};
		#my $var_id  = ${$record}{'ID'};
		my $ref      = ${$record}{'REF'};
		my $alts     = ${$record}{'ALT'};
		#my $qual    = ${$record}{'QUAL'};
		my $filters  = ${$record}{'FILTER'};
		#my $infos   = ${$record}{'INFO'};
		my $variant  = $seq_id . ':g.' . $pos . $ref . '>';
		my %gt_counts;
		my @gt_order;
		$stats{variants_parsed}++;
		
		#
		# Check VCF filters
		#
		if (scalar(keys(%{$vcf_filters}) > 0)) {
			$logger->debug('Filtering on VCF FILTER fields...');
			my $keep_record = 0;
			#my @applied_filters = split(';', $filters);
			CHECK_FILTERS:foreach my $applied_filter (@{$filters}) {
				$logger->debug('Checking if we should keep applied VCF FILTER ' . $applied_filter . '...');
				if (defined(${$vcf_filters}{$applied_filter})) {
					$keep_record = 1;
					$logger->debug('YES; keep record based on applied VCF FILTER ' . $applied_filter . '.');
					last CHECK_FILTERS;
				}
			}
			next VCF_RECORD unless($keep_record);
			$logger->debug('Keeping variant ' . $variant . '.');
		}
		
		if ($strip_samples) {
			#
			# Remove INFO:NS and INFO:SF subfields as those also provide details of the "source".
			#
			foreach my $key (keys(%{${$record}{'INFO'}})) {
				if ($key eq 'NS' || $key eq 'SF') {
					delete ${$record}{'INFO'}{$key};
					$logger->debug('Dropped INFO:' . $key . 'field.');
				}
			}
		}
		
		if ($strip_info > 0) {
			#
			# Drop all INFO sub fields except AC and AN as those two are automagically updated
			# when samples are removed/added, but the others are not. 
			#
			foreach my $key (keys(%{${$record}{'INFO'}})) {
				unless ($key eq 'AC' || $key eq 'AN' || defined(${$info_sub_fields_to_keep}{$key})) {
					delete ${$record}{'INFO'}{$key};
					$logger->debug('Dropped INFO:' . $key . 'field.');
				}
			}
		}
		
		#
		# Initialize hash and array for genotype counts.
		#
		# Length of hash and array depends on total number of alt alleles seen.
		# For example: 0/0,0/1,1/1,0/2,1/2,2/2,0/3,1/3,2/3,3/3,etc.
		# Hash is used to count genotypes, whereas the array is used to control 
		# the order of the genotypes (== keys of the hash).
		#
		my $this_allele_id = 0;
		my @added_allele_ids = (0);
		$gt_counts{'0/0'} = 0;
		push(@gt_order, '0/0');
		for (1..scalar(@{$alts})) {
			$this_allele_id++;
			foreach my $added_allele_id (@added_allele_ids) {
				my $gt = $added_allele_id . '/' . $this_allele_id;
				$gt_counts{$gt} = 0;
				push(@gt_order, $gt);
			}
			my $gt = $this_allele_id . '/' . $this_allele_id;
			$gt_counts{$gt} = 0;
			push(@gt_order, $gt);
			push(@added_allele_ids, $this_allele_id);
		}
		
		#
		# Parse samples for this variant.
		#
		SAMPLE:for my $sample (keys(%{${$record}{'gtypes'}})) {
			my $gt = ${$record}{'gtypes'}{$sample}{'GT'};	
			$logger->debug("\t$sample: GT=$gt");
			# Ditch phasing info.
			$gt =~ s'\|'/';
			if ($gt =~ m'\.') {
				$logger->trace("\t" . $sample . ': Skipping GT with one or more uncalled alleles (' . $gt . ').');
				next SAMPLE;
			} elsif ($gt =~ m/^(\d+)\/(\d+)$/) {
				my $allele1 = $1;
				my $allele2 = $2;
				# Re-order alleles to make sure the second allele's ID number >= first allele's ID number. 
				if ($allele2 < $allele1) {
					$gt = $allele2 . '/' . $allele1;
				}
				$logger->trace("\t$sample: GT=$gt (unphased and reordered).");
				if (defined($gt_counts{$gt})) {
					$logger->trace("\tGTC for $gt was:        $gt_counts{$gt}.");
				} else {
					$logger->fatal("\t" . 'GTC for ' . $gt . ' was:        undefined.');
					$logger->fatal("\t" . 'Cannot increment count for unexpected genotype with undefined count.');
					$logger->fatal("\t" . 'Sample ' . $sample . ': GT in unsupported format (' . $gt . ').');
					$logger->fatal("\t" . 'Defined GTCs for this variant (' . $seq_id . ':' . $pos . '):');
					foreach my $hkey (keys(%gt_counts)) {
						$logger->fatal("\t\t" . 'GTC ' . $hkey . ' = ' . $gt_counts{$hkey} . '.');
					}
					$logger->fatal('Processing ABORTED and hence output incomplete!');
					exit(1);
				}
				$gt_counts{$gt}++;
				$logger->debug("\tGTC for $gt updated to: $gt_counts{$gt}.");	
			} else {
				$logger->fatal("\t" . 'Sample ' . $sample . ': GT in unsupported format (' . $gt . ').');
				$logger->fatal('Processing ABORTED and hence output incomplete!');
				exit(1);
			}
		}
		
		#
		# Create INFO:GTC field with genotype counts in the right order. 
		#
		my $info_gtc;
		foreach my $gt (@gt_order) {
			$info_gtc .= $gt_counts{$gt} . ',';
		}	
		$info_gtc =~ s/,$//;
		$logger->debug("\tINFO:GTC = $info_gtc.");
		
		#
		# Append INFO:GTC field and format VCF record.
		#
		my $record_new;
		${$record}{'INFO'}{'GTC'} = $info_gtc;
		$record_new = $vcf->format_line($record, \@samples_included);
		
		#
		# If we only want to keep the summary data:
		# Keep only the first 8 columns == Remove FORMAT and sample columns.
		#
		if ($strip_samples) {
			$record_new =~ m/(([^\t]+\t){7}[^\t]+)/;
			$record_new = $1 . "\n";
		}
		
		#
		# Write new record to new VCF file.
		#
		print($path_to_fh $record_new);
		$logger->debug('New VCF record: ' . $record_new);
		$stats{'variants_stored'}++;
	}
	
	$vcf->close();
	return (\%stats);
	
}

sub _Usage {

	print STDERR "\n" .
	  'vcf-fill-gtc.pl:' . "\n" .
	  '   This script calculates GenoType Counts for samples in a VCF file and adds them to a GTC sub-field of the INFO field.' . "\n" .
	  '   The VCF can be filtered on the fly for certain samples using their IDs. ' . "\n\n" .
	  'Usage:' . "\n" .
	  '   vcf-fill-gtc.pl [options].' . "\n\n" .
	  'Available options are:' . "\n" .
	  '   -vcfi [file]   Input file in VCF format.' . "\n" .
	  '   -vcfo [file]   Output file in VCF format.' . "\n" .
	  '   -fs   [file]   Filter Samples: file containing sample IDs that shoud be filtered for.' . "\n" .
	  '                  (One ID per line).' . "\n" .
	  '                  Note that this will affect INFO subfields: see the -si and -ki options below.' . "\n" .
	  '   -d             When specied the sample IDs listed in the filter file will be dropped from the resulting VCF.' . "\n" .
	  '                  (Default is to keep the samples listed in the filter file).' . "\n" .
#	  '   -ml [string]   Match Logic that defines how the sample IDs from the filter file will be.' . "\n" .
#	  '                  matched to those in the VCF file. Supported logic types:' . "\n" .
#	  '                     literal     for exact matching.' . "\n" .
#	  '                     regex       for fuzzy matching using simple regular expressions (in Perl regex syntax).' . "\n" .
	  '   -si            Strip INFO subfields except AN, AC and other keys specified with -ki: removes others INFO subfields from the VCF.' . "\n" .
	  '                  This will happen automatically when removing samples with -fs and cannot be disabled in that case,' . "\n" .
	  '                  because only the AN and AC are updated automatically when adding or removing samples.' . "\n" .
	  '                  When samples were added or removed in a pre-processing step and values in the INFO subfields ' . "\n" . 
	  '                  (except for AN, AC and the ones specified with -ki) are no longer correct, these fields can be removed with -si.' . "\n" . 
	  '   -ki            Keep specific INFO subfields.' . "\n" . 
	  '                  Comma separated list of INFO subfield keys that in addition to AN and AC will be preserved when -fs or -si are used.' . "\n" . 
	  '   -ss            Strip Samples: removes detailed sample info from the VCF.' . "\n" .
	  '                  This retains only the first 8 "standard" VCF columns and removes the FORMAT column as well as all following sample columns.' . "\n" .
	  '   -fv   [string] Filter on VCF FILTERs: Comma separated list of filters.' . "\n" .
	  '                  Hence this will filter on the FILTER (7th) column of the VCF file.' . "\n" .
	  '                  Variants that match at least one of the filters will be preserved in the output.' . "\n" .
	  '   -ll   [LEVEL]  Log4perl log level. One of: ALL, TRACE, DEBUG, INFO (default), WARN, ERROR, FATAL or OFF.' . "\n" .
	  "\n";
	exit;

}
