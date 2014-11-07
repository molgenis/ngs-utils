#!/usr/bin/env perl

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
my $hgnc_gene_column_label = 'GENE';	# Default column label.
my $hgnc_id_column_label   = 'HGNC ID';	# Default column label.
my $output;
my $target = 'CT';						# Default target.
my $log_level = 'INFO';					# Default log level.


Getopt::Long::GetOptions (
	'i=s'		=> \$input,
	'g=s'		=> \$hgnc_gene_column_label,
	'h=s'		=> \$hgnc_id_column_label,
	'o=s'		=> \$output,
	't=s'		=> \$target,
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
# Reset target to default if user specified illegal target.
unless (defined($targets{$target})) {
	$logger->warn('Unknown target specified (' . $target . '). Resetting target to default: CT.');
	$target = 'CT';
}

my ($counter_found, $counter_missing, $counter_discarded) = _Parse($target, $input, $hgnc_gene_column_label, $hgnc_id_column_label, $output);

$logger->info('===================================================================');
$logger->info('Genes Symbols found and represented in BED file: ' . $counter_found . '.');
$logger->info('Genes Symbols discarded and lacking in BED file: ' . $counter_discarded . '. (Discarded because sequence not part of target ' . $target . ').');
$logger->info('Genes Symbols missing and lacking in BED file:   ' . $counter_missing . '.');
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
	  . '                  HGNC Gene Symbols.' . "\n"
	  . '                  HGNC IDs.' . "\n"
	  . '   -g           Column label in header line for the HGNC Gene Symbols. Default = GENE.' . "\n"
	  . '   -h           Column label in header line for the HGNC IDs. Default = HGNC ID.' . "\n"
	  . '   -o [file]    Output file in tab delimited BED format.' . "\n"
	  . '   -t [target]  Target regions of the genes listed in the input to cover in the output BED file. One of:' . "\n"
	  . '                  all     All complete exons.' . "\n"
	  . '                  CDS     Only the protein coding parts of the exons.' . "\n"
	  . '                          Hence excluding UTR.' . "\n"
	  . '                  CT      All complete exons, but only of the protein Coding Transcripts.' . "\n"
	  . '                          Hence including UTR.' . "\n"
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
	
	my ($target, $file_path_in,  $col_label_hgnc_gs, $col_label_hgnc_id, $file_path_out) = @_;
	
	$logger->info('Parsing HGNC Gene Symbols from ' . $file_path_in . '...');
	$logger->info('and storing (genomic) regions for targets "' . $target . '" to ' . $file_path_out . '...');
	
	my $file_path_in_fh;
	my $file_path_out_fh;
	my $line_counter = 0;
	my $genes_found = 0;
	my $genes_missing = 0;
	my $genes_found_but_discarded = 0;
	
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
	
	my $gene_adaptor = $registry->get_adaptor('Human', 'Core', 'Gene');
	my $col_offset_hgnc_id;
	my $col_offset_hgnc_gs;
	
	LINE: while (my $line = <$file_path_in_fh>) {
		
		$line_counter++;
		$logger->debug('Parsing record/line ' . $line_counter . '.');
		
		#
		# Skip empty lines and chomp newline, return and formfeed chars.
		#
		$line =~ s/[\r\n\f]+//g;
		next LINE if ($line =~ m/^$/);
		
		#
		# Parse header line.
		#
		if ($line_counter == 1) {
			
			$line =~ s/^#+//g;
			my @header_values = split("\t", $line);
			for my $offset (0 .. $#header_values) {
				
				$logger->trace('Column ' . $header_values[$offset] . ' has offset ' . $offset . '.');
				
				if ($header_values[$offset] eq "$col_label_hgnc_id") {
					$col_offset_hgnc_id = $offset;
					$logger->debug('Found HGNC IDs     in column with offset ' . $col_offset_hgnc_id . '.');
				} elsif ($header_values[$offset] eq "$col_label_hgnc_gs") {
					$col_offset_hgnc_gs = $offset;
					$logger->debug('Found Gene Symbols in column with offset ' . $col_offset_hgnc_gs . '.');
				}
			}
			
			unless (defined($col_offset_hgnc_id)) {
				$logger->fatal('Could not find column named "' . $col_label_hgnc_id . '" in header line.');
				exit 1;
			}
			unless (defined($col_offset_hgnc_gs)) {
				$logger->fatal('Could not find column named "' . $col_label_hgnc_gs . '" in header line.');
				exit 1;
			}
			
			next LINE;
			
		}
		
		my @column_values	= split("\t", $line);
		my $hgnc_symbol		= $column_values[$col_offset_hgnc_gs];
		my $hgnc_id			= $column_values[$col_offset_hgnc_id];
		#
		# Remove optional HGNC ID prefix
		#
		$hgnc_id =~ s/^HGNC://i;
		
		#
		# Get the gene from Ensembl core database.
		#
		
		#
		# The name is known as the gene's display_xref, which is a DBEntry object.
		# The following call $gene->display_xref->display_id should display the name like on the web page.
		#
		# my $all_genes=$gene_adaptor->fetch_all_by_external_name($key,'HGNC');
		# -> if $gene->slice->is_reference will skip all genes which are on patches and the like.
		#
		# my $hgnc_xrefs = $gene->get_all_DBLinks('HGNC');
		# foreach my $xref (@$hgnc_xrefs) {
		# 	print $xref->display_id . "\n";
		# }
		#
		# my $name = 'CLK2';
		# my ($gene) = grep { $_->display_xref()->display_id() eq $name } @{$gene_adaptor->fetch_all_by_external_name($name, 'HGNC')};
		#
		#
		
		#
		# 
		# If there was a gene object with an official HGNC Gene Symbol as displayname 
		# it will be in most cases the preferred/first gene returned,
		# but fetch_by_display_label may fetch somthing else than the gene with the official HGNC Gene Symbol
		# if another gene object with higher "priority" happed to have a display label value that is not the HGNC Gene symbol we were looking for, 
		# but happened to be the same.
		# This is problematic for example for gene "ACE", which links to multiple ENSG0000xxx IDs: ENSG00000264813 and ENSG00000159640.
		# One of which has an display label "ACE" based on an HGNC Gene Symbol and the other based on a UniProt gn field.
		#
		#my $gene = $gene_adaptor->fetch_by_display_label($hgnc_symbol);
		#
		# This will fetch only entries based on HGNC Gene Symbol.
		# The extra check using "$_->display_xref()->display_id() eq $hgnc_symbol" is necessary to filter for synomyms/aliases 
		# that only differ in case sensitivity compared to an approved gene symbol as the fetch_all_by_external_name method 
		# returns case insensitive results. 
		#
		my @genes = grep($_->display_xref()->display_id() eq $hgnc_symbol, @{$gene_adaptor->fetch_all_by_external_name($hgnc_symbol, 'HGNC')});
		
		#
		# Count the total amount of genes located anywhere:
		#
		my $ensembl_genes_found = scalar(@genes);
		#
		# Count the amount of genes located on large chromosomal reference sequences.
		# Hence counter excludes the ones located on extrachromosomal contigs, alternate assemblies, patches, etc.:
		#
		my $ensembl_genes_found_located_on_ref_chr = 0;
		#
		# Boolean flag to flag wheter any sequence region was added to the *.BED file for this HGNC Gene Symbol.
		#
		my $gene_represented_in_bed = 0;
		
		if ($ensembl_genes_found == 0) {
			$logger->error('Cannot find gene ' . $hgnc_symbol . ' in Ensembl core database!');
			$genes_missing++;
			$logger->trace('Genes missing incremented to ' . $genes_missing);
			next LINE;
		} elsif ($ensembl_genes_found == 1) {
			$logger->trace('Found one Ensembl gene for ' . $hgnc_symbol . '.');
		} elsif ($ensembl_genes_found > 1) {
			$logger->debug('Multiple Ensembl genes found for HGNC Gene Symbol ' . $hgnc_symbol . '.');
		} else {
			$logger->fatal('Cannot parse Ensembl search result for HGNC Gene Symbol ' . $hgnc_symbol . '.');
			exit 1;
		}
		
		GENE:foreach my $gene (@genes) {
			
			my $gene_id = $gene->stable_id();
			$logger->debug('Found Ensembl ID ' . $gene_id . ' for gene ' . $hgnc_symbol . '.');
			
			if ($gene->slice->is_reference) {
				$logger->trace("\t" . 'Ensembl ID ' . $gene_id . ' for gene ' . $hgnc_symbol . ' is located on a chromosomal reference sequence.');
				$ensembl_genes_found_located_on_ref_chr++;
			} else {
				$logger->debug("\t" . 'Skipping Ensembl ID ' . $gene_id . ' for gene ' . $hgnc_symbol . ', because it is located on an extrachromosomal contig, alternate assembly, patch or the like: ' . $gene->slice->seq_region_name() . '.');
				next GENE;
			}
			
			my $exon_count = 0;
			
			#
			# Retrieve transcript info.
			#
			foreach my $transcript (@{$gene->get_all_Transcripts()}) {
				
				my ($tr_id, $tr_seq, $tr_start, $tr_end, $tr_strand) = _GetFeatureDetails($transcript);
				$logger->trace("\t" . join("\t", $tr_id, $tr_seq, $tr_start, $tr_end, $tr_strand));
				
				#
				# Retrieve relevant exons based on requested target for this transcript.
				#
				my $exons = [];
				if ($target eq 'all') {
					$exons = $transcript->get_all_Exons();
				} elsif ($target eq 'CDS') {
					$exons = $transcript->get_all_translateable_Exons();
				} elsif ($target eq 'CT') {
					my $biotype = $transcript->biotype();
					if ($biotype eq 'protein_coding') {
						$exons = $transcript->get_all_Exons();
					} else {
						$logger->warn("\t" . 'Skipping transcript ' . $tr_id . ', because biotype is not protein_coding, but ' . $biotype . '.');
					}
				} else {
					$logger->fatal('Unknown target type specified: ' . $target);
					exit 1;
				}
				
				#
				# Parse exon info and store regions in BED file.
				#
				foreach my $exon (@{$exons}) {
					
					my ($ex_id, $ex_seq, $ex_start, $ex_end, $ex_strand) = _GetFeatureDetails($exon);
					$logger->trace("\t\t" . join("\t", $ex_id, $ex_seq, $ex_start, $ex_end, $ex_strand));
					
					#
					# BED starts are zero-based and BED ends are one-based.
					# Hence we need to substract 1 from the exon start to cover the complete region in the BED file.
					#
					
					my $bed_seq			= $ex_seq;
					my $bed_start		= $ex_start - 1;
					my $bed_end			= $ex_end;
					my $bed_annotation	= join('|', 'HGNC_Symbol=' . $hgnc_symbol, 'HGNC_ID=' . $hgnc_id, $ex_id, $tr_id, $gene_id);
					
					#
					# Store region in TAB delimited BED file.
					#
					print $file_path_out_fh join("\t", $bed_seq, $bed_start, $bed_end, $bed_annotation) . "\n";
					
					$exon_count++;
					$gene_represented_in_bed = 1;
				}
			}
			
			#
			# Check if any exons for any of this Ensembl gene object's transcripts were represented in our BED file.
			#
			if ($exon_count > 0) {
				if ($ensembl_genes_found > 1) {
					$logger->debug("\t" . 'Found target (' . $target . ') for gene ' . $hgnc_symbol . ' (' . $gene_id . ') and appended target regions to the output BED file.');
				}
			}
		}
		
		#
		# Check if:
		#  * any of the Ensembl gene objects for this HGNC Gene Symbol is (partially) represented in our BED file.
		#  * not more than one Ensembl gene object on a chromosomal reference sequences was used for this HGNC Gene Symbol.
		#
		if ($gene_represented_in_bed == 1) {
			$genes_found++;
			$logger->trace('Genes found incremented to ' . $genes_found);
			if ($ensembl_genes_found_located_on_ref_chr > 1) {
				$logger->warn('Gene ' . $hgnc_symbol . ' represented by more than one Ensembl gene in output BED file: check for annotation errors!');
			}
		} else {
			$genes_found_but_discarded++;
			$logger->warn('No target (' . $target . ') found for gene ' . $hgnc_symbol . '. Therefore gene not represented in output BED file.');
		}
	}
	
	close($file_path_in_fh);
	close($file_path_out_fh);
	
	return($genes_found, $genes_missing, $genes_found_but_discarded);
	
}
