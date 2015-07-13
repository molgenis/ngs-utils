#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use File::Glob ':glob';
use Getopt::Long;
local $| = 1; # autoflush stdout

#############################
##USAGE: perl collect_alignment_stats.pl -w <working_directory> -o <output_file>
#############################
# WBKoetsier
# Feb / Mar 2011
# v3
##
# This script collects metrics produced by the Picard tools CollectAlignment-
# SummaryMetrics (ASM), MarkDuplicates and, if applicable, CalculateHsMetrics
# (HSM). The metrics are assembled into a comma separated QC report (with similar
# columns as the BGI report) and written to the given outputfile.
# To print only the Picard documentation add the -P flag.
## 
# This script has been optimised for the in-house data that now (Feb/Mar 2011)
# exist in /data/wbkoetsier/data/Results/, which means that the working directory
# has the following structure: it contains at least one 'outputs' directory,
# that can have any name, as long as it starts with 'output'. The outputs
# directories contain lane directories (the in-house data prior to 110114 is
# one sample per lane) with the data analysis results. The script can be given
# a run directory or one or more of it's outputs directories.
# For (almost) all samples, FastQC has been run prior to the alignment. After 
# alignment, ASM and HSM have been run. Then follows MarkDuplicates, resulting
# in dedup metrics. After realignment and after recalibration ASM and HS were 
# repeated. This script collects ASM and HSM from after recalibration, allthough 
# this can be easily changed (look for the $version and $vversion variables).
# The output is comma-separated, because this is the most straight-forward text 
# format to import into Excel (for Win, 2003, the default for most clients). It
# can be changed by checking out all separators when joining the metrics and 
# header arrays.
##
# "All right, but apart from the sanitation, medicine, education, wine, public order, irrigation,
# roads, the fresh water system and public health, what have the Romans ever done for us?"
# (John Cleese as Reg in 'Life of Brian')
##

##### commandline input

## options
my ( $help, $working_dir, $output_file, $picard_docs);
GetOptions(
		"h"		=> \$help,
		"w=s"	=> \$working_dir,
		"o=s"	=> \$output_file,
		"P"		=> \$picard_docs
		);
		
# -P: print Picard docs only

## read command line
usage() and exit(0) if $help;
print_picard_documentation() and exit(0) if $picard_docs;
usage() and exit(0) unless $working_dir;
usage() and exit(0) unless $output_file;

## prepare the input
# chomp, allthough... is this neccesary?
chomp $output_file;
chomp $working_dir;
# strip trailing slash
$working_dir =~ s!\/$!!;
# check if the dir is really a dir
unless ( -d $working_dir ) { print "$working_dir is not a directory." and exit(0); };
# check (crudely) if we're dealing with working- or outputs- dir
$working_dir = $working_dir . '/output*' if $working_dir !~ /output/;

##### global vars
## the only global var :P
my $printable_metrics = '';


##### the script
# the LANEDIR foreach loop iterates over each lane directory. Per dir the samplename 
# and dirname are extracted using (the first) .sai file. Then the metrics file names 
# are extracted and counted. 
# If a file is missing, it is reported (STDOUT), the program continues nontheless. 
# Next, the data are gathered by calling a subroutine for each metrics file. The
# header and data are assembled using separate subroutines. Make sure to change both
# if columns need to change.
##
LANEDIR: foreach my $lanedir ( split "\n", `ls -d $working_dir/*/` ) {

	## ignore jobs and logs dirs (if encountered in the first place)
	next if ( $lanedir =~ /jobs/ || $lanedir =~ /logs/ );
	print "Lanedir: $lanedir, ";


	## list the files for this lane and extract samplename
	my $ls = `ls $lanedir`;
	# get samplename, and pe or se lane?
	my $samplename = '';
	my $pe = 0;
	if ( $ls =~ /(.+?)\.(.+?)(\.\d)*\.sai\n/ ) {
		$samplename = $1;
		print "Sample: $samplename\n";
		$pe = 1 if $3;
	} else {
		# this can happen!
		print "\tno data\n";
		next LANEDIR;
	}
	
	## get all metrics file names
	# fastqc and dedup occur only at one place in the pipeline
	# FastQC file(s): for paired end data, there are two fastqcs
	my ( $fastqsumm_file, $fastqsumm_file_pair1, $fastqsumm_file_pair2, $dedupmet_file ) = 
		( "", "", "", "" );

	if ( $ls =~ /(.+s_\d_sequence_fastqcsummary.txt)/ ) { # single end
		# runinfo_Llane_s_lane_sequence_fastqcsummary.txt
		$fastqsumm_file = $lanedir . $1;

	} elsif ( $ls =~ /(.+s_\d_\d_fastqcsummary.txt)/ ) { # paired end
		# runinfo_Llane_s_lane_pair_fastqcsummary.txt
		$fastqsumm_file_pair1 = $lanedir . $1 if ( $ls =~ /(.+s_\d_1_fastqcsummary.txt)/ ); # first of pair
		$fastqsumm_file_pair2 = $lanedir . $1 if ( $ls =~ /(.+s_\d_2_fastqcsummary.txt)/ ); # second of pair

	} # non of the above? missing fastqc.
	
	# dedup metrics file
	$dedupmet_file = $lanedir . $1 if $ls =~ /(.+dedup.metrics)/;
	
	# check numbers, tell if fastqc or dedup metrics is missing
	my $count = 0;
	my $miss = "\tMissing:";
	$miss .= " fastqcsummary" and $count++ 
		unless ( !$fastqsumm_file != !( !(!$fastqsumm_file_pair1 && !$fastqsumm_file_pair2) ) );
		# xor: I invite you to draw a truth table for a!=(!(b&&c)) (okay, whatever, it works)
	$miss .= " dedupmetrics" and $count++ if !($dedupmet_file);
	
	# there could be more than one version of ASM or HSM
	# add 'v' for version to the hash key to avoid an empty string as hash key
	my ( %asmetrics_filenames, %hsmetrics_filenames ) = ( (), () ); 
	foreach my $version ( "", "realigned", "recal" ) {

		my ( $asummet_file, $hsmet_file ) = ( "", "" );
		$asummet_file = $lanedir . $1 if $ls =~ /(.+$version.AlignmentSummaryMetrics)/;
		$asmetrics_filenames{ 'v' . $version } = $asummet_file;
		
		$hsmet_file  = $lanedir . $1 if $ls =~ /(.+$version.HsMetrics)/;
		$hsmetrics_filenames{ 'v' . $version } = $hsmet_file;
		
		# check numbers, tell if something is missing
		$miss .= " $version.AlignmentSummaryMetrics" and $count++ if !($asummet_file);	
		$miss .= " $version.HsMetrics" and $count++ if !($hsmet_file);
	}

	$miss .= "\n";
	# print missing files
	print $miss if $count;
	
	## speedup?? ls is no longer needed
	undef ( $ls );
		
	# location and sample
	my %sampleinfo = (); my %sampleinfo1 = (); my %sampleinfo2 = ();
	$lanedir =~ s|\/$||;
	$sampleinfo{ 'LOCATION' } = $lanedir;
	$sampleinfo{ 'SAMPLE' } = '';
	if ( $pe ) {
		$sampleinfo1{ 'LOCATION' } = 'first of pair';
		$sampleinfo1{ 'SAMPLE' } = '';
		$sampleinfo2{ 'LOCATION' } = 'second of pair';
		$sampleinfo2{ 'SAMPLE' } = '';
	}

	# fastqc (3 lines for pe)
	my %fqc_hash = (); my %fqc_hash1 = (); my %fqc_hash2 = ();
	if ( $pe ) {
		my $f = 0;
		if ( $fastqsumm_file_pair1 ) {
			$f += 1;
			%fqc_hash1 = get_fastqc_metrics( $fastqsumm_file_pair1 );
		} 
		if ( $fastqsumm_file_pair2 ) {
			$f += 2;
			%fqc_hash2 = get_fastqc_metrics( $fastqsumm_file_pair2 );
		}
		# missing columns?
		if ( $f == 1 ) { #only pair1
		
		} elsif ( $f == 2 ) { #only pair2
		
		} elsif ( $f == 0 ) { #neither
		
		} # f3-both

	} elsif ( $fastqsumm_file ) {
		%fqc_hash = get_fastqc_metrics( $fastqsumm_file );
	}
	
	# dedup (1 line)
	my %dedup_hash = ();
	%dedup_hash = get_dedup_metrics( $dedupmet_file ) if $dedupmet_file;
	
	my $vversion = 'vrecal'; # ( remember, a v precedes the version)
	
	# asm (3 lines for pe)
	my ( $asm_hash_ref, $asm_hash1_ref, $asm_hash2_ref );
	my %asm_hash = (); my %asm_hash1 = (); my %asm_hash2 = ();
	if ( $pe ) {
		( $asm_hash_ref, $asm_hash1_ref, $asm_hash2_ref ) = get_as_metrics( $pe, $asmetrics_filenames{ $vversion } ) 
			if $asmetrics_filenames{ $vversion } ne "";
		
		%asm_hash = %$asm_hash_ref;
		%asm_hash1 = %$asm_hash1_ref;
		%asm_hash2 = %$asm_hash2_ref;

	} else {
		%asm_hash = get_as_metrics( $pe, $asmetrics_filenames{ $vversion } ) 
			if $asmetrics_filenames{ $vversion } ne "";
	}
	
	# hsm (1 line)
	my %hsm_hash = ();
	%hsm_hash = get_hs_metrics( $hsmetrics_filenames{ $vversion } )
		if $hsmetrics_filenames{ $vversion } ne "";
	
	## assemble the lines
	if ( $pe ) {

		$printable_metrics .= assemble_columns( \%sampleinfo,
												\%fqc_hash,
												\%dedup_hash,
												\%asm_hash,
												\%hsm_hash,
												\%sampleinfo1, \%sampleinfo2,
												\%fqc_hash1, \%fqc_hash2,
												\%asm_hash1, \%asm_hash2);
												
	} else {

		$printable_metrics .= assemble_columns( \%sampleinfo,
												\%fqc_hash,
												\%dedup_hash,
												\%asm_hash,
												\%hsm_hash );	
	}
	

} # end foreach LANEDIR

##### the output

## assemble printable two-line header
my $printable_header = assemble_header();

## print data to output file
open (OUTPUT, ">$output_file") or die "unable to open $output_file for reading: $!\n";
print OUTPUT "$printable_header\n$printable_metrics";
print OUTPUT "\n\n\n";
close OUTPUT;

##### END #####

##### subs

## assemble columns
# takes hashes in fixed order
# each hash containing key-value pairs for metrics read from a metrics file
# places a value (based on the key) at specific positions in arrays for 
# each line (1 line for se, but 3 for pe!)
# prints NA for missing values, or an empty string for data that wasn't
# supposed to be present anyway.
# returns the line(s) as a string ending with \n
sub assemble_columns {

	my @pmetrics = ();
	my $metrics = '';
	my $pe = 0;

	my %sampleinfo = (); my %sampleinfo1 = (); my %sampleinfo2 = ();
	my %fqc_hash = (); my %dedup_hash = ();	my %asm_hash = (); my %hsm_hash = ();
	my %fqc_hash1 = ();	my %fqc_hash2 = ();	my %asm_hash1 = ();	my %asm_hash2 = ();

	## recieve hashes
	%sampleinfo = %{shift @_};
	%fqc_hash = %{shift @_};
	%dedup_hash = %{shift @_};
	%asm_hash = %{shift @_};
	%hsm_hash = %{shift @_};
	if ( scalar(@_)) {
		# if these hashes are present, we're dealing with pe data
		$pe = 1;
		%sampleinfo1 = %{shift @_};
		%sampleinfo2 = %{shift @_};
		%fqc_hash1 = %{shift @_};
		%fqc_hash2 = %{shift @_};
		%asm_hash1 = %{shift @_};
		%asm_hash2 = %{shift @_};
	}
	
	## some data may be missing!
	# keys:
	my ( $ks, $ks1, $ks2, $kf, $kf1, $kf2, $kd, $ka, $ka1, $ka2, $kh ) = 
		( 0, 0, 0, 0, 0, 0, 0, 0, 0 );
	$ks = 1 if keys( %sampleinfo );
	$ks1 = 1 if keys( %sampleinfo1 );
	$ks2 = 1 if keys( %sampleinfo2 );
	$kf = 1 if keys( %fqc_hash );
	$kf1 = 1 if keys( %fqc_hash1 );
	$kf2 = 1 if keys( %fqc_hash2 );
	$kd = 1 if keys( %dedup_hash );
	$ka = 1 if keys( %asm_hash );
	$ka1 = 1 if keys( %asm_hash1 );
	$ka2 = 1 if keys( %asm_hash2 );
	$kh = 1 if keys( %hsm_hash );
	
	
	## line 1
	## main
	# sampleinfo
	if ( $ks ) {
		push @pmetrics, ( $sampleinfo{ 'LOCATION' },
						  $sampleinfo{ 'SAMPLE' } );
	} else {
		push @pmetrics, ( "NA", "NA" );
	}
	# fastqc
	if ( $kf ) {
		push @pmetrics, ( $fqc_hash{ 'MEAN_QUALITY_LWR_27' },
						  $fqc_hash{ 'GC_CONTENT' },
						  $fqc_hash{ 'NON_UNIQUE_READS' } );
	} else {
 		push @pmetrics, ( "NA", "NA", "NA" ) if !$pe;
 		push @pmetrics, ( " ", " ", " " ) if $pe;
	}
	# asm/dedup/hs
	if ( $kh ) { 
		push @pmetrics, $hsm_hash{ 'TARGET_TERRITORY' };
	} else {
		push @pmetrics, "NA";
	}
	if ( $ka ) {
		push @pmetrics, ( $asm_hash{ 'PF_READS' },
						  ( $asm_hash{ 'PF_READS' } * $asm_hash{ 'MEAN_READ_LENGTH' } ),
						  $asm_hash{ 'PF_HQ_ALIGNED_READS' } );
	} else {
		push @pmetrics, ( "NA", "NA", "NA" );
	}
	if ( $kh && $ka ) {
		my $otb = $hsm_hash{ 'ON_TARGET_BASES' };
		my $mrl = $asm_hash{ 'MEAN_READ_LENGTH' };
		if ( $otb eq "NA" ||
			 $otb eq 0 ||
			 $mrl eq "NA" ||
			 $mrl eq 0 )
		{
			push @pmetrics, "NA";
		} else {
			push @pmetrics, ( $otb / $mrl );
		}
	} else {
		push @pmetrics, "NA";
	}
	if ( $kh ) {
		push @pmetrics, ( $hsm_hash{ 'ON_TARGET_BASES' },
						  $hsm_hash{ 'MEAN_TARGET_COVERAGE' } );
		my $otb = $hsm_hash{ 'ON_TARGET_BASES' };
		my $tat = $hsm_hash{ 'TARGET_TERRITORY' };
		if ( $otb eq "NA" ||
			 $otb eq 0 ||
			 $tat eq "NA" ||
			 $tat eq 0 )
		{
			push @pmetrics, "NA";
		} else {
			push @pmetrics, ( $otb / $tat );
		}
	} else {
		push @pmetrics, ( "NA", "NA", "NA" );
	}
	if ( $ka ) {
		push @pmetrics, ( $asm_hash{ 'MEAN_READ_LENGTH' },
						  $asm_hash{ 'PF_HQ_ERROR_RATE' } );
	} else {
		push @pmetrics, ( "NA", "NA" );
	}
	if ( $kh ) {
		push @pmetrics, ( $hsm_hash{ 'PCT_TARGET_BASES_2X' },
						  $hsm_hash{ 'PCT_TARGET_BASES_10X' },
						  $hsm_hash{ 'PCT_TARGET_BASES_20X' },
						  $hsm_hash{ 'PCT_TARGET_BASES_30X' },
						  $hsm_hash{ 'PCT_SELECTED_BASES' } );
	} else {
		push @pmetrics, ( "NA", "NA", "NA", "NA", "NA" );
	}
	if ( $kh && $ka ) {
		my $obb = $hsm_hash{ 'ON_BAIT_BASES' };
		my $nbb = $hsm_hash{ 'NEAR_BAIT_BASES' };
		my $mrl = $asm_hash{ 'MEAN_READ_LENGTH' };
		if ( $obb eq "NA" || $obb eq 0 ||
			 $nbb eq "NA" || $nbb eq 0 ||
			 $mrl eq "NA" || $mrl eq 0 )
		{
			push @pmetrics, "NA";
		} else {
			push @pmetrics, ( ($obb + $nbb) / $mrl );
		}
	} else {
		push @pmetrics, "NA";
	}
	if ( $kd ) {
		push @pmetrics, $dedup_hash{ 'UNPAIRED_READ_DUPLICATES' };
	} else {
		push @pmetrics, "NA";
	}
	## rest
	if ( $ka ) {
		push @pmetrics, ( $asm_hash{ 'CATEGORY' },
						  $asm_hash{ 'PF_NOISE_READS' },
						  $asm_hash{ 'PF_HQ_MEDIAN_MISMATCHES' },
						  $asm_hash{ 'READS_ALIGNED_IN_PAIRS' },
						  $asm_hash{ 'PCT_READS_ALIGNED_IN_PAIRS' },
						  $asm_hash{ 'BAD_CYCLES' },
						  $asm_hash{ 'STRAND_BALANCE' },
						  $asm_hash{ 'PCT_CHIMERAS' },
						  $asm_hash{ 'PCT_ADAPTER' } );
	} else {
		push @pmetrics, ( "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA" );
	}
	if ( $kh ) {
		push @pmetrics, ( $hsm_hash{ 'ZERO_CVG_TARGETS_PCT' },
						  $hsm_hash{ 'PCT_USABLE_BASES_ON_TARGET' } );
	} else {
		push @pmetrics, "NA";
	}
	if ( $kd ) {
		push @pmetrics, ( $dedup_hash{ 'UNPAIRED_READS_EXAMINED' },
						  $dedup_hash{ 'READ_PAIRS_EXAMINED' },
						  $dedup_hash{ 'READ_PAIR_DUPLICATES' },
						  $dedup_hash{ 'READ_PAIR_OPTICAL_DUPLICATES' } );
		if ( exists $dedup_hash{ 'ESTIMATED_LIBRARY_SIZE' } ) {
 			push @pmetrics, $dedup_hash{ 'ESTIMATED_LIBRARY_SIZE' };
		} else {
			push @pmetrics, "NA";
		}
	} else {
		push @pmetrics, ( "NA", "NA", "NA", "NA", "NA" );
	}

	$metrics .= $_ . "," foreach @pmetrics;
	$metrics .= "\n";
	
	@pmetrics = ();
	
	if ( $pe ) {
		
		## line 2
		## main
		# sampleinfo
		if ( $ks1 ) {
			push @pmetrics, ( $sampleinfo1{ 'LOCATION' },
							  $sampleinfo1{ 'SAMPLE' } );
		} else {
			push @pmetrics, ( "NA", "NA" );
		}
		# fastqc
		if ( $kf1 ) {
			push @pmetrics, ( $fqc_hash1{ 'MEAN_QUALITY_LWR_27' },
							  $fqc_hash1{ 'GC_CONTENT' },
							  $fqc_hash1{ 'NON_UNIQUE_READS' } );
		} else {
			push @pmetrics, ( "NA", "NA", "NA" );
		}
		# asm/dedup/hs
		# hsm TARGET_TERRITORY
		push @pmetrics, " ";
		if ( $ka1 ) {
			push @pmetrics, ( $asm_hash1{ 'PF_READS' },
							  ( $asm_hash1{ 'PF_READS' } * $asm_hash1{ 'MEAN_READ_LENGTH' } ),
							  $asm_hash1{ 'PF_HQ_ALIGNED_READS' } );
		} else {
			push @pmetrics, ( "NA", "NA", "NA" );
		}
		if ( $kh && $ka1 ) {
			my $otb = $hsm_hash{ 'ON_TARGET_BASES' };
			my $mrl = $asm_hash1{ 'MEAN_READ_LENGTH' };
			if ( $otb eq "NA" ||
				 $otb eq 0 ||
				 $mrl eq "NA" ||
				 $mrl eq 0 )
			{
				push @pmetrics, "NA";
			} else {
				push @pmetrics, ( $otb / $mrl );
			}
		} else {
			push @pmetrics, "NA";
		}
		# hsm ON_TARGET_BASES, MEAN_TARGET_COVERAGE', (ON_TARGET_BASES / TARGET_TERRITORY)
		push @pmetrics, ( " ", " ", " " );
		if ( $ka1 ) {
			push @pmetrics, ( $asm_hash1{ 'MEAN_READ_LENGTH' },
							  $asm_hash1{ 'PF_HQ_ERROR_RATE' } );
		} else {
			push @pmetrics, ( "NA", "NA" );
		}
		# hsm PCT_TARGET_BASES_2X, PCT_TARGET_BASES_10X, PCT_TARGET_BASES_20X,
		# PCT_TARGET_BASES_30X, PCT_USABLE_BASES_ON_TARGET
		push @pmetrics, ( " ", " ", " ", " ", " " );
		if ( $kh && $ka1 ) {
			my $obb = $hsm_hash{ 'ON_BAIT_BASES' };
			my $nbb = $hsm_hash{ 'NEAR_BAIT_BASES' };
			my $mrl = $asm_hash1{ 'MEAN_READ_LENGTH' };
			if ( $obb eq "NA" || $obb eq 0 ||
				 $nbb eq "NA" || $nbb eq 0 ||
				 $mrl eq "NA" || $mrl eq 0 )
			{
				push @pmetrics, "NA";
			} else {
				push @pmetrics, ( ($obb + $nbb) / $mrl );
			}
		} else {
			push @pmetrics, "NA";
		}
		# dedup UNPAIRED_READ_DUPLICATES
		push @pmetrics, " ";
		## rest
		if ( $ka1 ) {
			push @pmetrics, ( $asm_hash1{ 'CATEGORY' },
							  $asm_hash1{ 'PF_NOISE_READS' },
							  $asm_hash1{ 'PF_HQ_MEDIAN_MISMATCHES' },
							  $asm_hash1{ 'READS_ALIGNED_IN_PAIRS' },
							  $asm_hash1{ 'PCT_READS_ALIGNED_IN_PAIRS' },
							  $asm_hash1{ 'BAD_CYCLES' },
							  $asm_hash1{ 'STRAND_BALANCE' },
							  $asm_hash1{ 'PCT_CHIMERAS' },
							  $asm_hash1{ 'PCT_ADAPTER' } );
		} else {
			push @pmetrics, ( "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA" );
		}
		# hsm ZERO_CVG_TARGETS_PCT, PCT_SELECTED_BASES
		push @pmetrics, " ";
		# dedup UNPAIRED_READS_EXAMINED, READ_PAIRS_EXAMINED, READ_PAIR_DUPLICATES,
		# READ_PAIR_OPTICAL_DUPLICATES, ESTIMATED_LIBRARY_SIZE
		push @pmetrics, ( " ", " ", " ", " ", " " );
		
		$metrics .= $_ . "," foreach @pmetrics;
		$metrics .= "\n";
		
		@pmetrics = ();
		
		## line 3
		## main
		# sampleinfo
		if ( $ks2 ) {
			push @pmetrics, ( $sampleinfo2{ 'LOCATION' },
							  $sampleinfo2{ 'SAMPLE' } );
		} else {
			push @pmetrics, ( "NA", "NA" );
		}
		# fastqc
		if ( $kf2 ) {
			push @pmetrics, ( $fqc_hash2{ 'MEAN_QUALITY_LWR_27' },
							  $fqc_hash2{ 'GC_CONTENT' },
							  $fqc_hash2{ 'NON_UNIQUE_READS' } );
		} else {
			push @pmetrics, ( "NA", "NA", "NA" );
		}
		# asm/dedup/hs
		# hsm TARGET_TERRITORY
		push @pmetrics, " ";
		if ( $ka2 ) {
			push @pmetrics, ( $asm_hash2{ 'PF_READS' },
							  ( $asm_hash2{ 'PF_READS' } * $asm_hash2{ 'MEAN_READ_LENGTH' } ),
							  $asm_hash2{ 'PF_HQ_ALIGNED_READS' } );
		} else {
			push @pmetrics, ( "NA", "NA", "NA" );
		}
		if ( $kh && $ka2 ) {
			my $otb = $hsm_hash{ 'ON_TARGET_BASES' };
			my $mrl = $asm_hash2{ 'MEAN_READ_LENGTH' };
			if ( $otb eq "NA" ||
				 $otb eq 0 ||
				 $mrl eq "NA" ||
				 $mrl eq 0 )
			{
				push @pmetrics, "NA";
			} else {
				push @pmetrics, ( $otb / $mrl );
			}
		} else {
			push @pmetrics, "NA";
		}
		# hsm ON_TARGET_BASES, MEAN_TARGET_COVERAGE', (ON_TARGET_BASES / TARGET_TERRITORY)
		push @pmetrics, ( " ", " ", " " );
		if ( $ka2 ) {
			push @pmetrics, ( $asm_hash2{ 'MEAN_READ_LENGTH' },
							  $asm_hash2{ 'PF_HQ_ERROR_RATE' } );
		} else {
			push @pmetrics, ( "NA", "NA" );
		}
		# hsm PCT_TARGET_BASES_2X, PCT_TARGET_BASES_10X, PCT_TARGET_BASES_20X,
		# PCT_TARGET_BASES_30X, PCT_USABLE_BASES_ON_TARGET
		push @pmetrics, ( " ", " ", " ", " ", " " );
		if ( $kh && $ka2 ) {
			my $obb = $hsm_hash{ 'ON_BAIT_BASES' };
			my $nbb = $hsm_hash{ 'NEAR_BAIT_BASES' };
			my $mrl = $asm_hash2{ 'MEAN_READ_LENGTH' };
			if ( $obb eq "NA" || $obb eq 0 ||
				 $nbb eq "NA" || $nbb eq 0 ||
				 $mrl eq "NA" || $mrl eq 0 )
			{
				push @pmetrics, "NA";
			} else {
				push @pmetrics, ( ($obb + $nbb) / $mrl );
			}
		} else {
			push @pmetrics, "NA";
		}
		# dedup UNPAIRED_READ_DUPLICATES
		push @pmetrics, " ";
		## rest
		if ( $ka2 ) {
			push @pmetrics, ( $asm_hash2{ 'CATEGORY' },
							  $asm_hash2{ 'PF_NOISE_READS' },
							  $asm_hash2{ 'PF_HQ_MEDIAN_MISMATCHES' },
							  $asm_hash2{ 'READS_ALIGNED_IN_PAIRS' },
							  $asm_hash2{ 'PCT_READS_ALIGNED_IN_PAIRS' },
							  $asm_hash2{ 'BAD_CYCLES' },
							  $asm_hash2{ 'STRAND_BALANCE' },
							  $asm_hash2{ 'PCT_CHIMERAS' },
							  $asm_hash2{ 'PCT_ADAPTER' } );
		} else {
			push @pmetrics, ( "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA" );
		}
		# hsm ZERO_CVG_TARGETS_PCT, PCT_SELECTED_BASES
		push @pmetrics, " ";
		# dedup UNPAIRED_READS_EXAMINED, READ_PAIRS_EXAMINED, READ_PAIR_DUPLICATES,
		# READ_PAIR_OPTICAL_DUPLICATES, ESTIMATED_LIBRARY_SIZE
		push @pmetrics, ( " ", " ", " ", " ", " " );
		
		$metrics .= $_ . "," foreach @pmetrics;
		$metrics .= "\n";
		
	} # end if $pe	
	
	return $metrics;
}

## assemble header
# if you change the header, change the data assembly as well!
# returns header as a string ending with \n
sub assemble_header {

	my @pheader1 = ();
	my @pheader2 = ();
	push @pheader1, ( 'Location',
					'Project/sample',
					'MEAN_QUALITY_LWR_27',
					'GC_CONTENT',
					'NON_UNIQUE_READS' );
	push @pheader1, ( 'TARGET_TERRITORY',
					'PF_READS',
					'PF_READS*MEAN_READ_LENGTH',
					'PF_HQ_ALIGNED_READS',
					'ON_TARGET_BASES / MEAN_READ_LENGTH',
					'ON_TARGET_BASES',
					'MEAN_TARGET_COVERAGE',
					'ON_TARGET_BASES / TARGET_TERRITORY',
					'MEAN_READ_LENGTH',
					'PF_HQ_ERROR_RATE',
					'PCT_TARGET_BASES_2X',
					'PCT_TARGET_BASES_10X',
					'PCT_TARGET_BASES_20X',
					'PCT_TARGET_BASES_30X',
					'PCT_SELECTED_BASES',
					'(ON_BAIT_BASES + NEAR_BAIT_BASES) / MEAN_READ_LENGTH',
					'UNPAIRED_READ_DUPLICATES' );
	push @pheader2, ( ' ', 
					' ',
					'Percentage of reads with base call accuracy < 99.8%',
					'GC content (%)',
					'Percentage of non-unique reads in the first 200000' );
	push @pheader2, ( 'target region (bp)',
					'Raw reads',
					'Raw data yield',
					'Reads mapped to genome',
					'Reads mapped to target region',
					'Data mapped to target region',
					'Mean depth of target region',
					'Coverage of target region (%)',
					'Average read length (bp)',
					'Rate of nucleotide mismatch',
					'Fraction of target covered >= 2X','Fraction of target covered >= 10X',
					'Fraction of target covered >= 20X','Fraction of target covered >= 30X',
					'Capture specificity',
					'Reads mapped to target and flanking regions',
					'Number of reads marked as duplicate' );
		# PCT_USABLE_BASES_ON_TARGET - Capture specificity: The number of aligned, de-duped, on-target bases out of the PF bases available
		# PCT_SELECTED_BASES: On+Near Bait Bases / PF Bases Aligned. - Exome: good selections > 80%
		# ZERO_CVG_TARGETS_PCT: The number of targets that did not reach coverage>=2X over any base. - exomes: typically < 3%
		# coverage: % of target bases (PCT_TARGET_BASES_20X) typically 80% for 20X
		# STRAND_BALANCE: The number of PF reads aligned to the positive strand of the genome divided by the number of PF reads aligned to the genome.
	push @pheader2, ( 'CATEGORY',
					'PF_NOISE_READS',
					'PF_HQ_MEDIAN_MISMATCHES',
					'READS_ALIGNED_IN_PAIRS',
					'PCT_READS_ALIGNED_IN_PAIRS',
					'BAD_CYCLES',
					'STRAND_BALANCE',
					'PCT_CHIMERAS',
					'PCT_ADAPTER',
					'ZERO_CVG_TARGETS_PCT',
					'PCT_USABLE_BASES_ON_TARGET',
					'UNPAIRED_READS_EXAMINED',
					'READ_PAIRS_EXAMINED',
					'READ_PAIR_DUPLICATES',
					'READ_PAIR_OPTICAL_DUPLICATES',
					'ESTIMATED_LIBRARY_SIZE' );

	my $phdr = '';
	$phdr .= "$_," foreach @pheader1;
	$phdr .= "\n";
	$phdr .= "$_," foreach @pheader2;
	return $phdr;
}



## gather FastQC metrics from given fastqcsummary.txt file
# returns hash
sub get_fastqc_metrics { 

	my $fastqc_file = shift;
	
	my %fastq_metrics = ();
	
	open ( FQC, "<$fastqc_file" ) or die "cannot open \"$fastqc_file\": $!\n";
	while ( <FQC> ) {
		$fastq_metrics{ "MEAN_QUALITY_LWR_27" } = $2 if 
			m|(.+?) reads \((.+?)% of total\) have a mean quality score lower than 27|;
		$fastq_metrics{ "GC_CONTENT" } = $1 if
			m|The mean GC content for all positions is (.+?)%|;
		$fastq_metrics{ "NON_UNIQUE_READS" } = $1 if
			m|(.+?)% \((.+?) reads\) of the first 200,000 reads are not unique.|;
	}
	close FQC;
	
	return %fastq_metrics;

}

## get dedup metrics
# in the dedup metrics file, lines starting with '#' are info. 
# Allthough it's very much possible to plot the histogram, we won't.
# (http://sourceforge.net/apps/mediawiki/picard/index.php?title=Main_Page#Q:_What_is_meaning_of_the_histogram_produced_by_MarkDuplicates.3F)
# returns hash
sub get_dedup_metrics {

	my $dedupmet_file = shift;
	
	my %dedup_metrics = ();
	my ( @header, @metrics );

	## read in the data
	open ( DUP, "<$dedupmet_file" ) or die "cannot open \"$dedupmet_file\": $!\n";
	while ( <DUP> ) {
		
		chomp;
		
		# ignore lines:
		next if /^#/;  # that's the file info header
		next if /^\s/; # empty lines start with whitespace

		# catch lines and omit anything after that:
		@header = split( /\s/, $_ ) and next if /^LIBRARY/; # that's the data header
		@metrics = split( /\s/, $_ ) and last;				# last: the metrics

	}
	close DUP;
	
	# if ESTIMATED_LIBRARY_SIZE is missing, stuff goes horribly wrong.
	# Assuming that element is always last in the dedup metrics:
	if ( 
		  ( $header[-1] eq 'ESTIMATED_LIBRARY_SIZE' ) &&
		  ( scalar( @header ) != scalar( @metrics ) )
	   ) { push @metrics, "NA"; };

	# replace other missing data (empty columns) with "NA"	
	map { s/^\s*$/NA/ } @metrics;

	## store in hash
	# use header for counter because metrics values may be missing 
	foreach ( 0..$#header ) {
		$dedup_metrics{ $header[ $_ ] } = $metrics[ $_ ];
	}
	
	return %dedup_metrics;
}

## get ASmetrics for given file
# in case of pe, the file will have three lines of metrics:
# for first of pair, second of pair and pair.
# return a hash for se, but hash refs for pe!!
sub get_as_metrics {

	my $pe = shift;
	my $asummet_file = shift;
	
	my %as_metrics = (); my %as_metrics1 = (); my %as_metrics2 = ();
	my ( @header, @metrics, @metrics1, @metrics2 );
	
	open ( ASM, "<$asummet_file" ) or die "cannot open \"$asummet_file\": $!\n";
	while ( <ASM> ) {

		# ignore lines:
		next if /^#/;  # that's the file info header
		next if /^\s/; # empty lines start with whitespace

		# catch lines:
		@header = split( /\s/, $_ ) and next if /^CATEGORY/; # that's the data header
		@metrics = split( /\s/, $_ ) and last if /^UNPAIRED/; # last: the metrics
		@metrics1 = split( /\s/, $_ ) if /^FIRST_OF_PAIR/;
		@metrics2 = split( /\s/, $_ ) if /^SECOND_OF_PAIR/;
		@metrics = split( /\s/, $_ ) if /^PAIR/;
	}
	close ASM;
	
	# replace missing data (empty columns) with "NA"
	map { s/^\s*$/NA/ } @metrics;
	map { s/^\s*$/NA/ } @metrics1;
	map { s/^\s*$/NA/ } @metrics2;
	
	## store in hash
	# use header for counter because metrics values may be missing
	if ( $pe ) {
		foreach ( 0..$#header ) {
			$as_metrics{ $header[ $_ ] } = $metrics[ $_ ];
			$as_metrics1{ $header[ $_ ] } = $metrics1[ $_ ];
			$as_metrics2{ $header[ $_ ] } = $metrics2[ $_ ];
		}
		return ( \%as_metrics, \%as_metrics1, \%as_metrics2 ); # return refs...
	} else {
		foreach ( 0..$#header ) {
			$as_metrics{ $header[ $_ ] } = $metrics[ $_ ];
		}
		return %as_metrics;
	}
		
}

## get HsMetrics for given file
# like asmetrics
sub get_hs_metrics {

	my $hsmet_file = shift;
	
	my %hs_metrics = ();
	my ( @header, @metrics );

	open ( HSM, "<$hsmet_file" ) or die "cannot open \"$hsmet_file\": $!\n";
	while ( <HSM> ) {
	
		# ignore lines:
		next if /^#/;  # that's the file info header
		next if /^\s/; # empty lines start with whitespace

		# catch lines:
		@header = split( /\s/, $_ ) and next if /^BAIT_SET/; # that's the data header
		@metrics = split( /\s/, $_ ) and last;				 # last: the metrics
	}
	close HSM;
	
	# replace missing data (empty columns) with "NA"	
	map { s/^\s*$/NA/ } @metrics;
	
	## store in hash
	# use header for counter because metrics values may be missing 
	foreach ( 0..$#header ) {
		$hs_metrics{ $header[ $_ ] } = $metrics[ $_ ];
	}
	
	return %hs_metrics;
}


## commify large numbers. Do not use in combination with comma
# separated output...
# usage: commify( $whatever_var_containing_large_number )
# http://www.perlmonks.org/?node_id=653 by 'faq_monk'
sub commify {
        local $_  = shift;
        1 while s/^(-?\d+)(\d{3})/$1,$2/;
        return $_;
}

	
# print docs on the Picard metrics specifically
sub print_picard_documentation {
print OUTPUT <<'EOF';
## Alignment Information ##
# metrics file: AlignmentSummaryMetrics
CATEGORY: whether the reads are paired or unpaired and, if applicable, which pair
TOTAL_READS: total number of reads in the sample
PF_READS: all reads in the sample that passed the Illumina filter (PF)
PCT_PF_READS: percentage of the total
PF_NOISE_READS: number of PF reads that are composed entirly out of As and/or Ns
PF_READS_ALIGNED: total number of aligned PF reads
PCT_PF_READS_ALIGNED: percentage of the PF reads
PF_HQ_ALIGNED_READS: number of PF reads that aligned with a Phred mapping quality of 20 or higher (or a chance of 1/100 or less that the aligner is wrong)
PF_HQ_ALIGNED_BASES: number of PF bases (with >Q20) that aligned to the reference genome. When significantly different from PF_HQ_ALIGNED_READS * MEAN_READ_LENGTH, many reads may be aligned with gaps
PF_HQ_ALIGNED_Q20_BASES: subset of the above with Phred base calling quality of 20 or higher
PF_HQ_MEDIAN_MISMATCHES: median number of mismatches in PF_HQ_ALIGNED_READS
PF_HQ_ERROR_RATE: percentage of bases that mismatch in PF_HQ_ALIGNED_READS
MEAN_READ_LENGTH: mean length of all reads (should be the read length when metrics are from a single lane)
READS_ALIGNED_IN_PAIRS: number of aligned reads that have a mate that is also aligned (0 when single end)
PCT_READS_ALIGNED_IN_PAIRS: percentage of the above
BAD_CYCLES: number of instrument cycles where 80% or more base calls were no-calls
STRAND_BALANCE: the proportion of PF reads aligned to the positive and negative strand (PF reads aligned to the positive strand / total number of aligned PF reads)
PCT_CHIMERAS: percentage of reads that map outside of a maximum insert size (usually 100kb) or that have the two ends mapping to different chromosomes
PCT_ADAPTER: percentage of unaligned PF reads that match to a known adapter sequence right from the start of the read.

## Duplicates Information ##
# metrics file: dedup.metrics
# duplicate reads are reads that have identical 5' end coordinates
UNPAIRED_READS_EXAMINED: number of mapped reads without a mapped mate pair the tool considered
READ_PAIRS_EXAMINED: number of mapped read pairs the tool considered
UNMAPPED_READS: number of unmapped reads the tool considered
UNPAIRED_READ_DUPLICATES: number of fragments marked as duplicate
READ_PAIR_DUPLICATES: number of read pairs marked as duplicate
READ_PAIR_OPTICAL_DUPLICATES: subset of the above that are optical duplicates (sequences from one cluster that were identified by the Illumina software as originating from an adjacent cluster)
PERCENT_DUPLICATION: percentage of reads that are marked as duplicate
ESTIMATED_LIBRARY_SIZE: estimated number of unique sequences in the library (paired end only)

## Hybrid selection and coverage Information ##
# metrics file: HsMetrics
GENOME_SIZE: number of bases in the reference sequence
BAIT_TERRITORY: number of bases with one or more bait on top of them
TARGET_TERRITORY: unique number of target (e.g. exons) bases
BAIT_DESIGN_EFFICIENCY: TARGET_TERRITORY/BAIT_TERRITORY
ON_BAIT_BASES: number of PF bases that mapped to a baited region of the genome
NEAR_BAIT_BASES: number of bases that mapped to within a fixed interval of the baited regions
OFF_BAIT_BASES: number of bases that mapped neither on or near a baited region
ON_TARGET_BASES: number of bases that mapped to a targeted region of the genome
PCT_SELECTED_BASES: number of on and near bait bases / PF aligned bases (>80% for a good exome set)
ON_BAIT_VS_SELECTED: percentage on bait bases of all selected bases
MEAN_BAIT_COVERAGE: mean coverage of all baits
MEAN_TARGET_COVERAGE: mean coverage of targets with a coverage of 2 or more at at least one base
PCT_USABLE_BASES_ON_BAIT: number of aligned, non duplicate, on bait PF bases
PCT_USABLE_BASES_ON_TARGET: number of aligned, non duplicate, on target PF bases
FOLD_ENRICHMENT: fold by which the baited region has been amplified above genomic background
ZERO_CVG_TARGETS_PCT: number of targets with no bases that have a coverage of 2 or more (<3% for a good exome set)
FOLD_80_BASE_PENALTY: measure of non-uniformity. Number of units of sequencing necessary to raise 80% of the bases to the mean coverage (<5, between 3 and 4 for a good exome set)
PCT_TARGET_BASES_NX: percentage of target bases with a coverage of at least N (typically 80% for 20X)
HS_LIBRARY_SIZE: as dedup metrics ESTIMATED_LIBRARY_SIZE
HS_PENALTY_NX: hybrid selection penalty. Penalty needed to raise 80% of the bases to a coverage of N. Measure for the amount of input sequence needed to get a target base to a coverage of N. Should be <10. When 0 or empty, the library cannot be sequenced to NX coverage

# Information:
# http://picard.sourceforge.net/picard-metric-definitions.shtml
# http://www.broadinstitute.org/files/shared/mpg/nextgen2010/nextgen_cibulskis.pdf

EOF
}

## usage
# print usage
sub usage {
	print <<EOF;
collect_alignment_stats_v3.pl collects metrics produced by the Picard tools 
CollectAlignmentSummaryMetrics, MarkDuplicates and, if applicable, CalculateHsMetrics.
The metrics are assembled into a comma separated QC report and written to the
given outputfile.
To print only the Picard documentation add the -P flag.

## USAGE: perl collect_alignment_stats_v3.pl -w <working_dir> -o <output_file>

To print only the Picard documentation add the -P flag.

EOF
 
}


