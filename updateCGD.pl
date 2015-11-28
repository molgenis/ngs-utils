#!/usr/bin/perl -w

#
# initialise environment
#
use strict;
use Getopt::Long;
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

my %formats = (
	'cgd'	=> 1,
);

#
# Get options.
#
my $input;
my $format = 'cgd'; 
my $output;
my $hgnc_custom_download;
my $hgnc_symbol_checker;
my $log_level = 'INFO';			# The default log level.

Getopt::Long::GetOptions (
	'i=s'		=> \$input,
	'f=s'		=> \$format,
	'o=s'		=> \$output,
	'cd=s'		=> \$hgnc_custom_download,
	'sc=s'		=> \$hgnc_symbol_checker,
	'll=s'		=> \$log_level,
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
$format = (defined($formats{$format}) ? $format : 'cgd');

#
# Check input files.
#
for my $file ($input, $hgnc_custom_download, $hgnc_symbol_checker) { 
	unless (-e $file && -f $file && -r $file) {
		$logger->fatal('Cannot read from input file ' . $file . ': ' . $!);
		exit;
	}
}

#
# Create lookup tables (hashes) of relevant Gene Symbols and related IDs, synonyms, etc.
#
my ($hgnc_ids, $hgnc_symbols) = _ParseHGNC($hgnc_custom_download, $hgnc_symbol_checker);

#
# Parse input file with Gene Symbols and 
# update the symbols to (newer) official HGNC Gene symbols
# where necessary.
#
if ($format eq 'cgd') {
	_UpdateCGD($input, $hgnc_ids, $hgnc_symbols, $output);
} else {
	$logger->fatal('Unsupported input file format ' . $format. '.');
	exit 1;
}
$logger->info('Finished!');

#
##
### Internal subs.
##
#

#
# Parse HGNC symbol info.
#
sub _ParseHGNC {

	my ($cd_file_path_in, $sc_file_path_in) = @_;
	
	#
	# Create lookup table for HGNC IDs -> EntrezGene IDs.
	#
	$logger->info('Parsing HGNC custom download output ' . $cd_file_path_in . '...');
	
	my $file_path_in_fh;
	my %hgnc2entrezgene_ids;
	my $line_counter = 0;
	
	eval {
		open($file_path_in_fh, "<$cd_file_path_in");
	};
	if ($@) {
		$logger->fatal('Cannot read input file: ' . $@);
		exit;
	}
	
	my $col_offset_hgnc_id;
	my $col_offset_hgnc_gs;
	my $col_offset_entrezgene_id_by_hgnc;
	my $col_offset_entrezgene_id_by_ncbi;
	
	CDLINE: while (my $line = <$file_path_in_fh>) {
		
		$line_counter++;
		$logger->debug('Parsing record/line ' . $line_counter . '.');
		
		$line =~ s/[\r\n\f]+//g;
		# Skip blank and header lines.
		next if ($line eq '');
		
		if ($line_counter == 1) {
			#
			# Parse header line.
			#
			my @header_values = split("\t", $line);
			for my $offset (0 .. $#header_values) {
				
				$logger->trace('Column ' . $header_values[$offset] . ' has offset ' . $offset . '.');
				
				if ($header_values[$offset] eq 'HGNC ID') {
					$col_offset_hgnc_id = $offset;
					$logger->debug('Found HGNC IDs               in column with offset ' . $col_offset_hgnc_id . '.');
				} elsif ($header_values[$offset] eq 'Approved Symbol') {
					$col_offset_hgnc_gs = $offset;
					$logger->debug('Found HGNC Gene Symbols      in column with offset ' . $col_offset_hgnc_gs . '.');
				} elsif ($header_values[$offset] eq 'Entrez Gene ID') {
					$col_offset_entrezgene_id_by_hgnc = $offset;
					$logger->debug('Found EntrezGene IDs by HGNC in column with offset ' . $col_offset_entrezgene_id_by_hgnc . '.');
				} elsif ($header_values[$offset] eq 'Entrez Gene ID(supplied by NCBI)') {
					$col_offset_entrezgene_id_by_ncbi = $offset;
					$logger->debug('Found EntrezGene IDs by NCBI in column with offset ' . $col_offset_entrezgene_id_by_ncbi . '.');
				}
			}
			
		} else {
			
			unless (defined($col_offset_hgnc_id)) {
				$logger->fatal('Could not find column named "HGNC ID" in header line.');
				exit 1;
			}
			unless (defined($col_offset_hgnc_gs)) {
				$logger->fatal('Could not find column named "Approved Symbol" in header line.');
				exit 1;
			}
			unless (defined($col_offset_entrezgene_id_by_hgnc)) {
				$logger->fatal('Could not find column named "Entrez Gene ID" in header line.');
				exit 1;
			}
			unless (defined($col_offset_entrezgene_id_by_ncbi)) {
				$logger->fatal('Could not find column named "Entrez Gene ID" in header line.');
				exit 1;
			}
			
			my @values = split("\t", $line);
			
			my $hgnc_id					= $values[$col_offset_hgnc_id];
			my $hgnc_gs					= $values[$col_offset_hgnc_gs];
			my $entrezgene_id_by_hgnc	= $values[$col_offset_entrezgene_id_by_hgnc];
			my $entrezgene_id_by_ncbi	= $values[$col_offset_entrezgene_id_by_ncbi];
			
			if (defined($entrezgene_id_by_hgnc) && $entrezgene_id_by_hgnc ne '') {
			
				$hgnc2entrezgene_ids{$hgnc_id}{'GeneSymbol'}			= $hgnc_gs;
				$hgnc2entrezgene_ids{$hgnc_id}{'EntrezGeneID'}			= $entrezgene_id_by_hgnc;
				$hgnc2entrezgene_ids{$hgnc_id}{'EntrezGeneID_source'}	= 'HGNC';
				$logger->debug('Found HGNC ID ' . $hgnc_id . ' with Gene Symbol ' . $hgnc_gs . ' and EntrezGene ID ' . $entrezgene_id_by_hgnc . '.');
				
			} elsif (defined($entrezgene_id_by_ncbi)) {
			
				$hgnc2entrezgene_ids{$hgnc_id}{'GeneSymbol'}   = $hgnc_gs;
				$hgnc2entrezgene_ids{$hgnc_id}{'EntrezGeneID'} = $entrezgene_id_by_ncbi;
				$hgnc2entrezgene_ids{$hgnc_id}{'EntrezGeneID_source'}	= 'NCBI';
				$logger->debug('Found HGNC ID ' . $hgnc_id . ' with Gene Symbol ' . $hgnc_gs . ' and EntrezGene ID ' . $entrezgene_id_by_ncbi . '.');
			
			} else {
				unless ($hgnc_gs =~ m/withdrawn/) {
					$logger->warn('HGNC ID ' . $hgnc_id . ' with Gene Symbol ' . $hgnc_gs . ' lacks an EntrezGene ID.');
				}
			}
		}
	}
	
	close($file_path_in_fh);
	
	$logger->info('Finished parsing HGNC custom download output ' . $cd_file_path_in . '.');
	
	$logger->info('Parsing HGNC symbol checker output ' . $sc_file_path_in . '...');
	
	my %hgnc_symbols;
	my %hgnc_ids;
	$line_counter = 0;
	my $translated_symbols_counter = 0;
	my $missing_symbols_counter = 0;
	
	
	eval {
		open($file_path_in_fh, "<$sc_file_path_in");
	};
	if ($@) {
		$logger->fatal('Cannot read input file: ' . $@);
		exit;
	}
	
	#
	# File format example lines:
	# 
	# Input	Match Type	Approved Symbol	Approved Name	HGNC ID	Location
	# A2M	matches	A2M	alpha-2-macroglobulin	HGNC:7	12p13.31
	#
	# Match Types:
	#
	# * Approved symbol - the input is the same as an approved symbol (note this search is case insensitive).
	# * Synonym - the input has been added by an HGNC editor as a synonym of the approved symbol.
	# * Previous symbol- the input was previously an approved gene symbol, but the gene has since been updated with the approved symbol shown.
	# * Entry withdrawn - the input was an approved symbol for a gene that has since been withdrawn.
	# * Unmatched - the input does not match or map to any HGNC approved gene symbol. 
	#               Note that there must be an exact character match between the input and matched term.
	#               For example, OCT7 will match as a synonym of POU3F2 but OCT-7 will not.
	#
	
	SCLINE: while (my $line = <$file_path_in_fh>) {
		
		$line_counter++;
		$logger->debug('Parsing record/line ' . $line_counter . '.');
		
		$line =~ s/[\r\n\f]+//g;
		# Skip blank and header lines.
		next if ($line eq '' | $line =~ m/^Input\t/);
		
		my @values = split("\t", $line);
		
		my $input_symbol			= $values[0];
		my $match_type				= $values[1];
		my $approved_symbol			= $values[2];
		my $approved_name			= $values[3];
		my $hgnc_id					= $values[4];
		my $location				= $values[5];
		my $entrezgene_id			= 'missing';
		my $entrezgene_id_source	= 'missing';
		
		if ($match_type eq 'Unmatched') {
			$missing_symbols_counter++;
			next SCLINE;
		} elsif ($match_type ne 'Approved symbol') { 
			#
			# Sometimes the input symbol is an alias, which only differs in uppercase/lowercase spelling, 
			# which may be lost in the output from the HGNC symbol checker...
			# making the alias and offical symbol exactly the same: skip those! 
			#
			if ($input_symbol eq $approved_symbol) {
				$logger->warn('Skipping alias ' . $input_symbol . ' for ' . $approved_symbol . ' (HGNC ID ' . $hgnc_id . ') on line ' . $line_counter . '.');
				next SCLINE;
			}
		}
		
		if (defined($hgnc_ids{$hgnc_id})) {
			unless ($hgnc_ids{$hgnc_id} eq $approved_symbol) {
				$logger->fatal('HGNC ID ' . $hgnc_id . ' already present/parsed. Duplicate detected on line ' . $line_counter . ':');
				$logger->fatal($line);
				exit 1;
			}
		} else {
			$hgnc_ids{$hgnc_id} = $approved_symbol;
			$translated_symbols_counter++;
		}
		
		#
		# Create truncated location for checking location of the gene with less resolution.
		#
		my $truncated_location;
		
		unless (defined($location)) {
			unless ($match_type =~ m/withdrawn/) {
				$logger->fatal('Cannot parse location for HGNC ID ' . $hgnc_id . ' on line ' . $line_counter . ':');
				$logger->fatal($line);
				exit 1;
			} else {
				$location = '';
			}
		}
		
		if ($location =~ m/(((([0-9]{1,2})|[XY])([pq][0-9]*)?)|M)/i) {
			$truncated_location = $1;
		} else {
			unless ($match_type =~ m/withdrawn/) {
				$logger->fatal('Cannot parse location for HGNC ID ' . $hgnc_id . ' on line ' . $line_counter . ':');
				$logger->fatal($line);
				exit 1;
			}
		}
		
		#
		# Get EntrezGene ID for this symbol.
		#
		if (defined($hgnc2entrezgene_ids{$hgnc_id}{'GeneSymbol'}) && defined($hgnc2entrezgene_ids{$hgnc_id}{'EntrezGeneID'})) {
			
			#
			# Sanity Check
			#
			if ($approved_symbol eq $hgnc2entrezgene_ids{$hgnc_id}{'GeneSymbol'}) {
				$entrezgene_id			= $hgnc2entrezgene_ids{$hgnc_id}{'EntrezGeneID'};
				$entrezgene_id_source	= $hgnc2entrezgene_ids{$hgnc_id}{'EntrezGeneID_source'};
			} else {
				$logger->fatal('Conflicting Approved Symbols for HGNC ID ' . $hgnc_id . ' on line ' . $line_counter . ': ' . $approved_symbol . ' != ' . $hgnc2entrezgene_ids{$hgnc_id}{'GeneSymbol'});
				$logger->fatal($line);
				exit 1;
			}
		}
		
		push(@{$hgnc_symbols{$input_symbol}}, {
			'match_type'			=> $match_type, 
			'symbol' 				=> $approved_symbol, 
			'hgnc_id'				=> $hgnc_id, 
			'entrezgene_id'			=> $entrezgene_id, 
			'entrezgene_id_source'	=> $entrezgene_id_source, 
			'truncloc'				=> $truncated_location, 
			'location'				=> $location}
		);
	
	}
	
	close($file_path_in_fh);
	
	$logger->info('Parsed input.');
	
	return (\%hgnc_ids, \%hgnc_symbols);
	
}

sub _UpdateCGD {
	
	my ($file_path_in, $hgnc_ids, $hgnc_symbols, $file_path_out) = @_;

	$logger->info('Parsing CGD ' . $file_path_in . '...');

	my $file_path_in_fh;
	my $file_path_out_fh;
	my $line_counter = 0;
	my $symbol_current = 0;			# Up-to-date gene symbols.
	my $symbol_updated_for_CS = 0;	# More or less the same symbol, but there was a Case Sensitivity issue.
	my $symbol_updated = 0;			# Outdated aliases, synonyms, homonyms, etc. that were updated to current HGNC symbols. 
	my $symbol_problematic = 0;		# Lost symbols, missing symbols, ambiguous symbols, ....
	my $entrezgene_id_curated_by_hgnc = 0;
	my $entrezgene_id_curated_by_ncbi = 0;
	
	eval {
		open($file_path_in_fh, "<$file_path_in");
	};
	if ($@) {
		$logger->fatal('Cannot read input file: ' . $@);
		exit;
	}
	eval {
		open($file_path_out_fh,	">$file_path_out");
	};
	if ($@) {
		$logger->fatal('Cannot write updated file: ' . $@);
		exit;
	}
	
	#
	# CGD original file format header:
	# #GENE	ENTREZ GENE ID	CONDITION	INHERITANCE	AGE GROUP	ALLELIC CONDITIONS	MANIFESTATION CATEGORIES	INTERVENTION CATEGORIES	COMMENTS	INTERVENTION/RATIONALE	REFERENCES
	#
	
	#
	# Write new header.
	#
	print $file_path_out_fh "#HGNC_Symbol	HGNC_ID	GENE	ENTREZ GENE ID	CONDITION	INHERITANCE	AGE GROUP	ALLELIC CONDITIONS	MANIFESTATION CATEGORIES	INTERVENTION CATEGORIES	COMMENTS	INTERVENTION/RATIONALE	REFERENCES\n";
	
	LINE: while (my $line = <$file_path_in_fh>) {

		$line_counter++;
		$logger->debug('Parsing record/line ' . $line_counter . '.');
		
		# Skip header and meta-data lines.
		next if ($line =~ m/^#/);
		$line =~ s/[\r\n\f]+//g;
		my $new_line;
		
		my @values = split("\t", $line);
		
		my $cgd_gene_symbol		= $values[0];
		my $cgd_entrezgene_id	= $values[1];
		
		if (scalar(@{${$hgnc_symbols}{$cgd_gene_symbol}}) == 1) {
			
			my $approved_symbol			= ${$hgnc_symbols}{$cgd_gene_symbol}[0]{'symbol'};
			my $hgnc_id					= ${$hgnc_symbols}{$cgd_gene_symbol}[0]{'hgnc_id'};
			my $entrezgene_id			= ${$hgnc_symbols}{$cgd_gene_symbol}[0]{'entrezgene_id'};
			my $entrezgene_id_source	= ${$hgnc_symbols}{$cgd_gene_symbol}[0]{'entrezgene_id_source'};
			my $match_type				= ${$hgnc_symbols}{$cgd_gene_symbol}[0]{'match_type'};
			
			unless ($cgd_entrezgene_id eq $entrezgene_id) {
				$logger->fatal('Single hit mismatch: CGD Gene Symbol only present once in HGNC data, but linked to different EntrezGene ID!');
				$logger->fatal('          CGD entry: ' . $line);
				$logger->fatal('   HGNC Gene Symbol: ' . $approved_symbol . ' // HGNC ID: ' . $hgnc_id . ' // EntrezGene ID: ' . $entrezgene_id . ' curated by ' . $entrezgene_id_source . '.');
				exit 1;
			}
			
			$new_line = $approved_symbol . "\t" . $hgnc_id . "\t" . $line;
			$logger->debug('Single hit: ' . $new_line);
			$logger->trace('    ----->: match_type: ' . $match_type);
			
			if ($match_type eq 'Approved symbol') {
				if ($cgd_gene_symbol ne $approved_symbol && $cgd_gene_symbol =~ m/$approved_symbol/i) {
					#
					# Updated for case sensitivity issue.
					#
					$symbol_updated_for_CS++;
				} else {
					#
					# Nothing to update: symbol is current.
					#
					$symbol_current++;
				}
			} else {
				$symbol_updated++;
				$logger->info('SH: ' . $cgd_gene_symbol . ' updated to HGNC Gene Symbol: ' . $approved_symbol . ' // HGNC ID: ' . $hgnc_id . ' // EntrezGene ID: ' . $entrezgene_id . ' curated by ' . $entrezgene_id_source . '.');
			}
			if ($entrezgene_id_source eq 'HGNC') {
				$entrezgene_id_curated_by_hgnc++;
			} elsif ($entrezgene_id_source eq 'NCBI') {
				$entrezgene_id_curated_by_ncbi++;
			} else {
				$logger->fatal('Unknown curator for EntrezGene ID in:');
				$logger->fatal('HGNC Gene Symbol: ' . $approved_symbol . ' // HGNC ID: ' . $hgnc_id . ' // EntrezGene ID: ' . $entrezgene_id . ' curated by ' . $entrezgene_id_source . '.');
				exit 1;
			}
			
		} elsif (scalar(@{${$hgnc_symbols}{$cgd_gene_symbol}}) > 1) {
			
			my $found_match=0;
			
			foreach my $entry (@{${$hgnc_symbols}{$cgd_gene_symbol}}) {
			
				my $approved_symbol			= ${$entry}{'symbol'};
				my $hgnc_id					= ${$entry}{'hgnc_id'};
				my $entrezgene_id			= ${$entry}{'entrezgene_id'};
				my $entrezgene_id_source	= ${$entry}{'entrezgene_id_source'};
				my $match_type				= ${$entry}{'match_type'};
				my $truncated_location		= ${$entry}{'truncloc'};
				my $location				= ${$entry}{'location'};
				
				if ($cgd_entrezgene_id == $entrezgene_id) {
					
					$logger->trace(' Multi hit match: CGD Gene Symbol present in HGNC data and linked to same EntrezGene ID!');
					$logger->trace('       CGD entry: ' . $line);
					$logger->trace('HGNC Gene Symbol: ' . $approved_symbol . ' // HGNC ID: ' . $hgnc_id . ' // EntrezGene ID: ' . $entrezgene_id . ' curated by ' . $entrezgene_id_source . '.');
					
					$found_match++;
					
					unless ($found_match = 1) {
						$logger->fatal('CGD entry alreay linked to HGNC! Redundant data present...');
						$logger->fatal('       CGD entry: ' . $line);
						$logger->fatal('HGNC Gene Symbol: ' . $approved_symbol . ' // HGNC ID: ' . $hgnc_id . ' // EntrezGene ID: ' . $entrezgene_id . ' curated by ' . $entrezgene_id_source . '.');
						exit 1;
					}
					
					if ($match_type eq 'Approved symbol') {
						if ($cgd_gene_symbol ne $approved_symbol && $cgd_gene_symbol =~ m/$approved_symbol/i) {
							#
							# Updated for case sensitivity issue.
							#
							$symbol_updated_for_CS++;
						} else {
							#
							# Nothing to update: symbol is current.
							#
							$symbol_current++;
						}
					} else {
						$symbol_updated++;
						$logger->info('MH: ' . $cgd_gene_symbol . ' updated to HGNC Gene Symbol: ' . $approved_symbol . ' // HGNC ID: ' . $hgnc_id . ' // EntrezGene ID: ' . $entrezgene_id . ' curated by ' . $entrezgene_id_source . '.');
					}
					if ($entrezgene_id_source eq 'HGNC') {
						$entrezgene_id_curated_by_hgnc++;
					} elsif ($entrezgene_id_source eq 'NCBI') {
						$entrezgene_id_curated_by_ncbi++;
					} else {
						$logger->fatal('Unknown curator for EntrezGene ID in:');
						$logger->fatal('HGNC Gene Symbol: ' . $approved_symbol . ' // HGNC ID: ' . $hgnc_id . ' // EntrezGene ID: ' . $entrezgene_id . ' curated by ' . $entrezgene_id_source . '.');
						exit 1;
					}
					
					$new_line = $approved_symbol . "\t" . $hgnc_id . "\t" . $line;
					$logger->debug('Multi hit: ' . $new_line);
					
				} else {
					$logger->trace('Multi hit mismatch: CGD Gene Symbol present in HGNC data, but linked to different EntrezGene ID!');
					$logger->trace('         CGD entry: ' . $line);
					$logger->trace('  HGNC Gene Symbol: ' . $approved_symbol . ' // HGNC ID: ' . $hgnc_id . ' // EntrezGene ID: ' . $entrezgene_id . ' curated by ' . $entrezgene_id_source . '.');
				}
			}
			
			unless ($found_match = 1) {
				
				$new_line = 'N/A' . "\t" . 'N/A' . "\t" . $line;
				$logger->warn('Gene Symbol from CGD "' . $cgd_gene_symbol . '" ambiguously present in HGNC, but not a single entry shares the same EntrezGene ID...');
				$symbol_problematic++;
			
			}
			
		} else {
			
			$new_line = 'N/A' . "\t" . 'N/A' . "\t" . $line;
			$logger->warn('Gene Symbol from CGD "' . $cgd_gene_symbol . '" does not exist in HGNC.');
			$symbol_problematic++;
			
		}
		
		print $file_path_out_fh $new_line . "\n";
		
	}
	
	close($file_path_in_fh);
	close($file_path_out_fh);
	
	$logger->info('Parsed CGD file.');
	$logger->info('==============');
	$logger->info('Problematic gene symbols:       ' . $symbol_problematic . '.');
	$logger->info('==============');
	$logger->info('Resolved gene symbols:');
	$logger->info('  Up-to-date HGNC gene symbols:................... ' . $symbol_current . '.');
	$logger->info('  Gene symbols with fixed case sensitivity issue:  ' . $symbol_updated_for_CS . '.');
	$logger->info('  Aliases / synomyms updated to HGNC gene symbols: ' . $symbol_updated . '.');
	$logger->info('  ==============');
	$logger->info('  HGNC IDs resolved via EntrezGene IDs curated by HGNC: ' . $entrezgene_id_curated_by_hgnc . '.');
	$logger->info('  HGNC IDs resolved via EntrezGene IDs curated by NCBI: ' . $entrezgene_id_curated_by_ncbi . '.');
	$logger->info('==============');
	
}

#
# Usage
#
sub _Usage {

	print STDERR "\n"
	  . 'updateCGD.pl:' . "\n\n"
	  . '   TODO.' . "\n\n"
	  . 'Usage:' . "\n\n"
	  . '   updateCGD.pl options' . "\n\n"
	  . 'Available options are:' . "\n\n"
	  . '   -i [file]    Input file in tab delimited format. Specify exact format details with -f' . "\n"
	  . '   -f [format]  Input file format. One of:' . "\n"
	  . '                  cgd - CGD Excel download exported to tab delimited *.txt. (default)' . "\n"
	  . '   -cd [file]   HGNC tab delimited "Custom Download" file' . "\n"
	  . '                  See: http://www.genenames.org/cgi-bin/download' . "\n"
	  . '                  Must contain at least the columns "HGNC ID", "Approved Symbol" and "Entrez Gene ID"' . "\n"
	  . '   -sc [file]   HGNC tab delimited file from online "Symbol Checker" tool.' . "\n"
	  . '                  See: http://www.genenames.org/cgi-bin/symbol_checker' . "\n"
	  . '   -o [file]    Output file in tab delimited format.' . "\n"
	  . '                Prepends two columns for updated HGNC gene symbols and HGNC gene IDs.' . "\n"
	  . '   -l [LEVEL]   Log4perl log level. One of: ALL, TRACE, DEBUG, INFO (default), WARN, ERROR, FATAL or OFF.' . "\n"
	  . "\n";
	exit;

}
