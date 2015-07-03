#!/usr/bin/perl -w

#
# initialise environment
#
use strict;
use Getopt::Std;
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

#
# Get options.
#
my %opts;
Getopt::Std::getopts('i:v:f:o:l:', \%opts);

my $input						= $opts{'i'};
my $vcf							= $opts{'v'};
my $fasta						= $opts{'f'};
my $output						= $opts{'o'};
my $log_level					= $opts{'l'};

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

	#{ level    => $log_level,
	#  file     => ">>FlattenTDF.log",
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
# Check input files.
#
for my $file ($input, $vcf, $fasta) { 
	unless (-e $file && -f $file && -r $file) {
		$logger->fatal('Cannot read from input file ' . $file . ': ' . $!);
		exit;
	}
}

#
# Parse input file with mutations.
#
my $mutations = _Parse($input);
my $refseq = _LoadRef($fasta);
_Check($vcf, $mutations, $refseq, $output);
$logger->info('Finished!');

#
##
### Internal subs.
##
#

sub _Parse {

	my ($file_path_in, $file_path_out) = @_;

	$logger->info('Parsing ' . $file_path_in . '...');

	my $file_path_in_fh;
	my $file_path_out_fh;
	my $line_counter = 0;
	
	eval {
		open($file_path_in_fh, "<$file_path_in");
	};
	if ($@) {
		$logger->fatal('Cannot read input file: ' . $@);
		exit;
	}
	
	my %hgmd_records;
	
	LINE: while (my $line = <$file_path_in_fh>) {

		$line_counter++;
		$logger->debug('Parsing record/line ' . $line_counter . '.');
		
		$line =~ s/[\r\n\f]+//g;
		next if ($line eq '');
		
		my @values = split("\t", $line);
		
		my $M_acc					= $values[0];
		my $HG19_chromosome			= $values[1];
		my $HG19_coordSTART			= $values[2];
		my $HG19_coordEND			= $values[3];
		my $HG19_strand				= $values[4];
		my $HG19_upstreamFLANK		= uc($values[5]);
		my $HG19_downstreamFLANK	= uc($values[6]);
		my $MN_wildBASE				= uc($values[7]);
		my $MN_mutBASE				= uc($values[8]);
		my $MN_hgvs					= $values[9];
		my $M_disease				= $values[10];
		my $M_gene					= $values[11];
		my $M_tag					= $values[12];
		
		if ($MN_wildBASE =~ m/^[ATCG]$/i && $MN_mutBASE =~ m/^[ATCG]$/i) {
			
			if ($MN_hgvs =~ m/([ATCG]{1})>([ATCG]{1})$/i) {
				
				my $hgvs_ref = uc($1);
				my $hgvs_alt = uc($2);
				
				unless ($MN_wildBASE eq $hgvs_ref && $MN_mutBASE eq $hgvs_alt) {
					
					$logger->fatal('ERROR1 REF <-> ALT mismatch on line ' . $line_counter);
					$logger->fatal('Ref: ' . $MN_wildBASE . ' // Alt: ' . $MN_mutBASE . ' // HGVS: ' . $MN_hgvs . '.');
					$logger->fatal($line);
					exit 1;
					
				}
				
			} elsif ($MN_hgvs =~ m/^(\\N)?$/i) {
				
				$logger->debug('HGMD HGVS annotation missing on line ' . $line_counter);
				$logger->trace($line);
				
			} else {
				
				$logger->fatal('HGMD HGVS annotation in unsupported format on line ' . $line_counter);
				$logger->fatal('Ref: ' . $MN_wildBASE . ' // Alt: ' . $MN_mutBASE . ' // HGVS: ' . $MN_hgvs . '.');
				$logger->fatal($line);
				exit 1;
				
			}
			
			my $coordinate = $HG19_chromosome . '_' . $HG19_coordSTART;
			
			#
			# Sanity check: check if reference a.k.a wild base is the same 
			# if multiple variant records are present for the same genomic position.
			#
			if (exists($hgmd_records{$coordinate})) {
				
				foreach my $record (@{$hgmd_records{$coordinate}}) {
					
					#my $MN_wildBASE_existing = ${$hgmd_records{$coordinate}[$record]{'MN_wildBASE'}};
					my $MN_wildBASE_existing = ${$record}{'MN_wildBASE'};
					
					$logger->debug('DUPLICATE mutation coordinates: ' . $coordinate . 
									' with wild base ' . $MN_wildBASE_existing . '=' . ${$record}{'MN_wildBASE'} . '.');
					
					if (${$record}{'HG19_strand'} ne $HG19_strand) {
						
						$MN_wildBASE_existing =~ tr/ATCG/TAGC/;
						
						if ($MN_wildBASE ne $MN_wildBASE_existing) {
							
							$logger->error('DUPLICATE mutation coordinates: ' . $coordinate . '!');
							$logger->error('Existing line: ' . ${$record}{'line'});
							$logger->error('Current line:  ' . $line);
							exit 1;
							
						}
					}
				}
			}
			
			push(@{$hgmd_records{$coordinate}}, {
				'M_acc' 				=> $M_acc, 
				'HG19_strand'			=> $HG19_strand, 
				'HG19_downstreamFLANK'	=> $HG19_downstreamFLANK, 
				'HG19_upstreamFLANK'	=> $HG19_upstreamFLANK, 
				'MN_wildBASE'			=> $MN_wildBASE, 
				'MN_mutBASE'			=> $MN_mutBASE, 
				'MN_hgvs'				=> $MN_hgvs, 
				'M_disease'				=> $M_disease, 
				'M_gene'				=> $M_gene, 
				'M_tag'					=> $M_tag, 
				'line'					=> $line}
			);
			
		}
	}
	
	close($file_path_in_fh);
	
	$logger->info('Parsed input.');
	
	return (\%hgmd_records);
	
}

sub _LoadRef {
	
	my ($fasta) = @_;
	
	$logger->info('Parsing FastA ' . $fasta . '...');
	
	my $file_path_fh;
	my $line_counter;
	my $sequence;

	eval {
		open($file_path_fh, "<$fasta");
	};
	if ($@) {
		$logger->fatal('Cannot read input file: ' . $@);
		exit;
	}
	
	LINE: while (my $line = <$file_path_fh>) {

		$line_counter++;
		$logger->trace('Parsing record/line ' . $line_counter . '.');
		
		# Skip header and meta-data lines.
		next if ($line =~ m/^>/);
		$line =~ s/[\r\n\f]+//g;
		$sequence .= $line;
		
	}
	
	close($file_path_fh);
	
	return(\$sequence);
}

sub _Check {
	
	my ($file_path_vcf, $mutations, $refseq, $file_path_out) = @_;

	$logger->info('Parsing VCF ' . $file_path_vcf . '...');

	my $file_path_vcf_fh;
	my $file_path_out_fh;
	my $line_counter = 0;
	my $shared_ref_allele_matches = 0;
	my $shared_variants = 0;
	my $ref_allele_mismatches = 0;
	my $ref_allele_putative_swap = 0;
	
	eval {
		open($file_path_vcf_fh, "<$file_path_vcf");
	};
	if ($@) {
		$logger->fatal('Cannot read input file: ' . $@);
		exit;
	}
	eval {
		open($file_path_out_fh,	">$file_path_out");
	};
	if ($@) {
		$logger->fatal('Cannot write converted file: ' . $@);
		exit;
	} else {
		# Write header.
		print $file_path_out_fh "#CHROM\tPOS\tREF\tALT\tHGMD_ACC\tHGMD_TAG\tHGMD_GENE\tHGMD_DISEASE\tHGMD_HGVS\n";
	}
	
	# VCF columns:
	#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
	
	LINE: while (my $line = <$file_path_vcf_fh>) {

		$line_counter++;
		$logger->debug('Parsing record/line ' . $line_counter . '.');
		
		# Skip header and meta-data lines.
		next if ($line =~ m/^#/);
		$line =~ s/[\r\n\f]+//g;
		
		my @values = split("\t", $line);
		
		my $vcf_chr = $values[0];
		my $vcf_pos = $values[1];
		my $vcf_ref = $values[3];
		my $vcf_alt = $values[4];
		
		my $vcf_coordinate = $vcf_chr . '_' . $vcf_pos;
		
		if (defined(${$mutations}{$vcf_coordinate})) {
			
			foreach my $record (@{${$mutations}{$vcf_coordinate}}) {
				
				my $hgmd_acc		= ${$record}{'M_acc'};
				my $hgmd_strand		= ${$record}{'HG19_strand'};
				my $hgmd_flank_up	= ${$record}{'HG19_upstreamFLANK'};
				my $hgmd_flank_down	= ${$record}{'HG19_downstreamFLANK'};
				my $hgmd_wild		= ${$record}{'MN_wildBASE'};
				my $hgmd_mut		= ${$record}{'MN_mutBASE'};
				my $hgmd_hgvs		= ${$record}{'MN_hgvs'};
				my $hgmd_disease	= ${$record}{'M_disease'};
				my $hgmd_gene		= ${$record}{'M_gene'};
				my $hgmd_tag		= ${$record}{'M_tag'};
				my $hgmd_line		= ${$record}{'line'};
				
				my $hgmd_wild_fwd 		= $hgmd_wild;
				my $hgmd_mut_fwd 		= $hgmd_mut;
				my $hgmd_flank_up_fwd	= $hgmd_flank_up;
				
				#
				# Reverse complement if mutation was located on the reverse strand.
				#
				if ($hgmd_strand eq '-') {
					
					$hgmd_wild_fwd =~ tr/ATCG/TAGC/;
					$hgmd_mut_fwd  =~ tr/ATCG/TAGC/;
					$hgmd_flank_up_fwd = $hgmd_flank_down;
					$hgmd_flank_up_fwd =~ tr/ATCG/TAGC/;
					$hgmd_flank_up_fwd = scalar(reverse($hgmd_flank_up_fwd));
					
				}
				
				my $flank_length = length($hgmd_flank_up_fwd);
				my $refseq_flank = substr(${$refseq}, ($vcf_pos - 1 - $flank_length), $flank_length);
				my $refseq_base  = substr(${$refseq}, ($vcf_pos - 1), 1);
				
				if ($vcf_ref ne $refseq_base) {
					$logger->fatal('Different ref alleles for coordinate: ' . $vcf_coordinate . '!');
					$logger->fatal('HGMD wild on + = ' . $hgmd_wild_fwd . ' // VCF Ref = ' . $vcf_ref . ' Refseq base = ' . $refseq_base . '.');
					$logger->fatal('HGMD line: ' . $hgmd_line);
					$logger->fatal('VCF  line: ' . $line);
					exit 1;
				}
				if ($hgmd_flank_up_fwd ne $refseq_flank) {
					$logger->fatal('Different flanking sequence for coordinate: ' . $vcf_coordinate . '!');
					$logger->fatal('HGMD wild on + = ' . $hgmd_wild_fwd . ' // VCF Ref = ' . $vcf_ref . ' Refseq base = ' . $refseq_base . '.');
					$logger->fatal('HGMD flank up on + = ' . $hgmd_flank_up_fwd . ' Refseq flank up on + = ' . $refseq_flank . '.');
					$logger->fatal('HGMD line: ' . $hgmd_line);
					$logger->fatal('VCF  line: ' . $line);
					exit 1;
				}
				
				# CHROM	POS	REF	ALT	HGMD_ACC	HGMD_TAG	HGMD_GENE	HGMD_DISEASE	HGMD_HGVS
				my $new_line = join("\t", $vcf_chr, $vcf_pos, $vcf_ref, $vcf_alt, $hgmd_acc, $hgmd_tag, $hgmd_gene, $hgmd_disease, $hgmd_hgvs);
				
				if ($vcf_ref eq $hgmd_wild_fwd) {
					
					$shared_ref_allele_matches++;
					
					if ($vcf_alt eq $hgmd_mut_fwd) {
					
						#
						# Both ref and alt allelel are shared for this genomic coordinate.
						#
						$shared_variants ++;
						print $file_path_out_fh $new_line . "\n";
						
						
					} else {
						
						$logger->debug('Different alt allele for matching ref allele for coordinate: ' . $vcf_coordinate . '!');
						$logger->debug('HGMD wild on + = ' . $hgmd_wild_fwd . ' // VCF Ref = ' . $vcf_ref . '.');
						$logger->debug('HGMD line: ' . $hgmd_line);
						$logger->debug('VCF  line: ' . $line);
						
					}
					
				} else {
					
					$logger->debug('Reference allele mismatch for coordinate: ' . $vcf_coordinate . '!');
					$logger->debug('HGMD wild on + = ' . $hgmd_wild_fwd . ' // VCF Ref = ' . $vcf_ref . '.');
					$logger->debug('HGMD line: ' . $hgmd_line);
					$logger->debug('VCF  line: ' . $line);
					$ref_allele_mismatches++;
					
					if (($vcf_ref eq $hgmd_mut_fwd) && ($vcf_alt eq $hgmd_wild_fwd)) {
						$ref_allele_putative_swap++;
						print $file_path_out_fh $new_line . "\n";
					}
				}
			}
			
		}
	}
	
	close($file_path_vcf_fh);
	close($file_path_out_fh);
	
	$logger->info('Parsed VCF.');
	$logger->fatal('Found ' . $shared_ref_allele_matches . ' shared variants with identical reference alleles.');
	$logger->fatal("\t" . $shared_variants . '/' . $shared_ref_allele_matches . ' of these shared variants also with identical alt alleles.');
	$logger->fatal('Found ' . $ref_allele_mismatches . ' shared variants with mismatching reference alleles.');
	$logger->fatal("\t" . $ref_allele_putative_swap . '/' . $ref_allele_mismatches . ' of these shared variants with putative ref<->alt allele swaps.');
	
}

#
# Usage
#

sub _Usage {

	print STDERR "\n"
	  . 'FilterHGMD.pl:' . "\n\n"
	  . '   TODO.' . "\n\n"
	  . 'Usage:' . "\n\n"
	  . '   FilterHGMD.pl options' . "\n\n"
	  . 'Available options are:' . "\n\n"
	  . '   -i [file]   Input file in tab delimited format. (HGMD SQL query into outfile)' . "\n"
	  . '   -v [file]   VCF file whose mutations should be compared to those of the input file.' . "\n"
	  . '   -f [file]   FastA file with the reference genome used for both the Input and VCF file.' . "\n"
	  . '               Is used to double check the variants using the flanking sequence ' . "\n"
	  . '               reported in the Input file.' . "\n"
	  . '   -o [file]   Output file in tab delimited format.' . "\n"
	  . '               Lists variants the Input and VCF file have in common.' . "\n"
	  . '               Can be used with vcftools to annotate the VCF file.' . "\n"
	  . '   -l [LEVEL]  Log4perl log level. One of: ALL, TRACE, DEBUG, INFO (default), WARN, ERROR, FATAL or OFF.' . "\n"
	  . "\n";
	exit;

}
