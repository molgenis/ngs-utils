#!/usr/bin/perl

# fastqc_report.pl
# this script takes the fastqc_data.txt that FastQC creates and translates it into a short, easily readable and useful text file.
# WBKoetier / WBruins Nov 2010

# changlelog v1: added per seq base content module, made minor changes to output of some other modules and changed the zip command to only extract the fastqc_data.txt file (replace if exists).
# changelog v2.1: Now gathers data in string and prints at once. Now contains some initial code to print to the end of the summary file some tab separated columns that report FastQCs module status for some modules. Added condinitional that causes the script to not mind wether -r has an extension or not.
# In the pipeline, this script is called 'v1' because of backwards compatability.

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Switch;

#### vars
# option vars
my ($help, $reportname, $path_to_archive, $output, $logfile);
# other global vars
my $lfh;
my $total_num_of_seqs = 0;
my $summary = '';
my $columns = ''; # per_base_seq_qual\tper_base_GC_content
####

#### get options
GetOptions(
		"h"		=> \$help,
		"r=s"	=> \$reportname,
		"p=s"	=> \$path_to_archive,
		"o=s"	=> \$output,
		"l=s"	=> \$logfile
		);
# help!
usage() and exit(0) if $help;
# mandatory args
usage() and exit(0) unless $reportname;
usage() and exit(0) unless $path_to_archive;
usage() and exit(0) unless $output;
# optional arg: logfile
if ( $logfile ) { open ($lfh, ">$logfile") or die "could not open $logfile [$!]\n"; }
# else { $lfh = \*STDOUT; } # replace with open /dev/null!
else { open ($lfh, ">/dev/null") or die "could not open /dev/null [$!]\n"; }
# and that's all, folks
speak ( "The following options remain unprocessed:" )  if $ARGV[0];
foreach (@ARGV) {
  print $lfh "$_\n";
}
####

#### input files and paths
speak ( "Verbose output for debug:" );
speak ( "The FastQC input reads file: $reportname" );

# FastQC will strip the extension off the input filename and append _fastqc
# to be save, strip path of reportname, strip trailing slash off archive path
if ($path_to_archive =~ /\/\z/ ) { $path_to_archive =~ s/\/\z//; };

# dir that FastQC writes
my $dirname;
if ( basename ($reportname) =~ /(.+)\.(\w+?)$/ ) { $dirname = "$1_fastqc"; }
else { $dirname = basename ($reportname); };

my $input = "$path_to_archive/$dirname";
# the archive
if ( -f "$input.zip" ) { print $lfh "The FastQC archive is $input.zip\n"; }
else { speak ( "There is no archive called $input.zip. Please check your input.") and exit(1) };
# did someone or something already unzip it for you?
my $unpacked = 0;
if ( -d $input ) { $unpacked = 1; speak ( "Archive $input.zip already unpacked."); }
# unzipping using the -o flag: this unzip overwrites!!
else { $unpacked = 0; print $lfh `unzip -u $input.zip $dirname/fastqc_data.txt -d $path_to_archive`; }

# check if the default datafile exists
if (!( -f "$input/fastqc_data.txt" )) {
	speak ( "the -f test failed on $input/fastqc_data.txt." );
}
####

#### reading the -r input file
# open the data file for reading.
open ( DATA, '<', "$input/fastqc_data.txt") or die "[FATAL] could not open $input/fastqc_data.txt for reading $!\n";

# it's not really a very big file. Will easily fit in one array.
speak ( "Reading in data..." );
my @data_modules;
{
	# set $/ () to the desired delimiter (for this block only)
	local $/ = ">>END_MODULE";
	@data_modules = <DATA>;
	# result: each element in @data_modules now contains >> followed by whatever data followed by >>END_MODULE
}

# remove first line of the fastqc_data.txt file
$data_modules[0] =~ s/##FastQC\s+(.+?)\n//;

# iterate over the modules: generate data for each.
foreach my $module ( @data_modules ) {

	my $key = $1 if $module =~ />>(.+?)\t/;
	next if !$key;
	
	switch ( $key ) {
		case /Basic Statistics/i { basic_statistics( $module ) }
		case /Per base sequence quality/i { perbase_seq_qual( $module ) } # status to columns
		case /Per sequence quality scores/i { perseq_qual( $module ) }
		case /Per base sequence content/i { perbase_seq_content( $module ) }
		case /Per base GC content/i { perbase_GC_content( $module ) } # status and % to columns
		case /Per sequence GC content/i { perseq_GC_content( $module ) } # currently unsupported
		case /Per base N content/i { perbase_N_content( $module ) }
		case /Sequence Length Distribution/i { seq_length_distr( $module ) }
		case /Sequence Duplication Levels/i { seq_duplication_levels( $module ) }
		case /Overrepresented sequences/i { overrepresented_seqs( $module ) }
		case /Kmer Content/i { kmer_content( $module ) } # currently unsupported
		else { speak ( "The switch failed. Was a new module added to the FastQC report, or did an existing module change names?" ) }

	}
}

close ( DATA );
speak ( "$input/fastqc_data.txt closed." );
####

#### writing to the -o output file

# open the outputfile for writing.
open ( OUT, '>', $output ) or die "[FATAL] could not open $output for writing. $!\n";

print OUT $summary;
# print OUT "\n\n";
# print OUT $columns;
####

close ( OUT );

speak ( "$output closed." );
####

#### Clean up after yourself...
# if this script was responsible for unpacking the FastQC report archive, remove the unpacked directory.
if ( $unpacked == 0 ) {
	`rm -rf $input`;
	speak ( "Removing $input." );
}
####

#### Yep, done.
speak ( "\nAll done!\n" );
close($lfh);

##### subroutines

## module data generation
# Each module starts with the module name and the status. The second line contains a header. All that follows (excluding the last line) is the tab delimited data. The last line defines the end of the module by stating >>END_MODULE.
# Provide this sub an extra variable (1 or whatever if true, optionally 0 or an empty string if false) if you only wish to report anything at warn or fail
sub initialise_module {
	my @lines = split "\n", shift;
	my $do_not_report = shift;
	
	# remove first element of @lines if it is an empty string
	shift @lines if $lines[0] eq '';

	# extract module name and status:
	my $module_start = shift @lines;
	my ($module, $status) = ($1, $2)
		if ( $module_start =~ />>(.+?)\t(.+?)\z/ );

	if ( $do_not_report ) {
		return ( ) # absolutely nothing if the module passed
			if ( $status =~ /pass/i );
	}

	speak ( "$module: $status" );
	$summary .= "\n### $module ###\n";
	#$summary .= "Module status: $status\n";
	
	# extract header (if it exists) and if it exists, the info line preceding the header
	# should this fail, check output and adapt this script
	my $header = shift @lines if ( $lines[0] =~ /^#/ );
	$header .= "\n" . shift @lines if ( $lines[0] =~ /^#/ );

	# the rest is the data and whatever follows it:
	my $data = join "\n", @lines;
	# return all but module name:
	return ( $header, $data, $status );
}

# basic statistics:
# print in readable format
sub basic_statistics {
 	my $module_string = shift;
 	my ( $header, $data, $status ) = initialise_module ( $module_string );
 	
	$summary .= "Original reads file name: $1\n"
		if $data =~ /Filename\t(.+?)\n/;
	$summary .= "Base call type: $1\n"
		if $data =~ /File type\t(.+?)\n/;
	
	# set global variable for number of sequences
	$total_num_of_seqs = $1
		if $data =~ /Total.+\t(.+?)\s*\n/;
	
	$summary .= "Total number of sequences in reads file: $total_num_of_seqs\n";
	$summary .= "Read length: $1 bases\n"
		if $data =~ /length\t(.+?)\s*\n/;
	$summary .= "Overall GC content: $1%\n"
		if $data =~ /GC\t(.+?)\s*\n/;
}

# per base sequence quality:
sub perbase_seq_qual {
 	my $module_string = shift;
 	my ( $header, $data, $status ) = initialise_module ( $module_string );
 	
 	$columns .= "$status\t";

	my @table = split "\n", $data;
	my ( @positions_q10, @positions_m25, @positions_q5, @positions_m20 ) = ((),(),(),());
	foreach my $row ( @table ) {
		my @r = split "\t", $row;
		last if $r[0] =~ />>END_MODULE/;
		# all cases: lower quartile less than 10
		push ( @positions_q10, $r[0] ) if $r[3] < 12;
		# warn: median score less than 25
		if ( $status =~ /warn/i ) {
			push ( @positions_m25, $r[0] ) if $r[2] < 25;
		}
		# fail: median score less than 20 and lower quartile less than 5
		if ( $status =~ /fail/i ) {
			push ( @positions_q5, $r[0] ) if $r[3] < 5;
			push ( @positions_m20, $r[0] ) if $r[2] < 20;
		}
	}

	my $str = "25% of all reads have a score of less than 10 at position";
	$str .= range_notation ( join " ", @positions_q10 );
	$summary .= "$str.\n";

	if ( $status =~ /warn/i ) {
	 	$str = "The median score is less than 25 at position";
		$str .= range_notation ( join " ", @positions_m25 );
 	 	$summary .= "$str (this causes FastQC to issue the 'warn' status).\n";
 	}
 	
 	if ( $status =~ /fail/i ) {
 		my $q5 = range_notation ( join " ", @positions_q5 );
 		my $m20 = range_notation ( join " ", @positions_m20 );
 		
 		if ( $q5 ) {
 		$str = "25% of all reads have a score of less than 5 at position";
	 	$str .= range_notation ( join " ", @positions_q5 );
 	 	$summary .= "$str (this causes FastQC to issue the 'fail' status).\n";
 	 	}
 	 	if ( $m20 ) {
 		$str = "The median score is less than 20 at position";
	 	$str .= range_notation ( join " ", @positions_m20 );
 	 	$summary .= "$str (this causes FastQC to issue the 'fail' status).\n";
 	 	}
 	}
}

# per sequence quality scores:
sub perseq_qual {
 	my $module_string = shift;
 	my ( $header, $data, $status ) = initialise_module ( $module_string );

	# create hash of the data: score as key, count as value
	my %data_hash;
	foreach ( split "\n", $data ) {
			$data_hash{ $1 } = $2 if $_ =~ /(\d+)\t(\d+\.*\d*)/;
	}
	
	# sort the values of the hash
	my @keys_sorted_by_counts = sort { $data_hash{$b} <=> $data_hash{$a} } keys %data_hash;
	
	# filter out the most observed score
	my $most_obs_score = $keys_sorted_by_counts[0];
	$summary .= "The most observed mean quality score is " . 
		$most_obs_score;

	# use this score to calculate error rate (10^(Q/-10)*100%)
	my $err_rate = ( 10 ** ( $most_obs_score / -10 ) ) * 100;
	$summary .= " (this equates to an error rate of ";
	$summary .= sprintf ( "%.2f", $err_rate );
	$summary .= "%)\n";

	# sum the counts for score < 27 (and no, there are only 33 entries in the hash)
	my $sum;
	foreach ( keys %data_hash ) {
		$sum += $data_hash{ $_ } if $_ < 27;
	}
	$summary .= "$sum reads (";
	$summary .= sprintf "%.2f", ( ($sum / $total_num_of_seqs)*100 );
	$summary .= "% of total) have a mean quality score lower than 27.\n";


}

# per base sequence content (A, C, T, G)
sub perbase_seq_content {
	my $module_string = shift;
	my ( $header, $data, $status ) = initialise_module ( $module_string, 1);

	my ( $warn, $fail ) = ( 0, 0 );

	if ( $data ) {

		my ( @positions_GC_10, @positions_AT_10, @positions_GC_20, @positions_AT_20 );
		foreach ( split "\n", $data ) {
			next if not /^\d/;

			# pos	G	A	T	C
 			my @bases = split "\t";
 			# discrepancies in the total number of bases may exist due to Ns
	 		my $total = $bases[1] + $bases[2] + $bases[3] + $bases[4];
	 		my $GC = ( $bases[1] / $total ) - ( $bases[4] / $total );
	 		my $AT = ( $bases[2] / $total ) - ( $bases[3] / $total );


			# please note that FastQC rounds off the percentage before checking the below, which may cause the status to be different than expected. Because of this rounding off, I will stick with the numbers and skip the status.
			
			# warn
			if ( $GC > .1 || $GC < -.1 ) {
				push @positions_GC_10, $bases[0];
				$warn = 1;
			}
			if ( $AT > .1 || $AT < -.1 ) {
				push @positions_AT_10, $bases[0];
				$warn = 1;
			}
		
			# fail
			if ( $GC > .2 || $GC < -.2 ) {
				push @positions_GC_20, $bases[0];
				$fail = 1;
			}
			if ( $AT > .2 || $AT < -.2 ) {
				push @positions_AT_20, $bases[0];
				$fail = 1;
			}
		}

		# warn
		my $str = "The difference between G and C is more than 10% of the total at position";
		$str .= range_notation ( join " ", @positions_GC_10 );
 		$summary .= "$str.\n" 
	 		if $warn == 0;
		$str = "The difference between A and T is more than 10% of the total at position";
		$str .= range_notation ( join " ", @positions_AT_10 );
	 	$summary .= "$str.\n" 
 			if $warn == 1;

		# fail
		$str = "The difference between G and C is more than 20% of the total at position";
		$str .= range_notation ( join " ", @positions_GC_20 );
	 	$summary .= "$str.\n"
 			if $fail == 1;
		$str = "The difference between A and T is more than 20% of the total at position";
		$str .= range_notation ( join " ", @positions_AT_20 );
 		$summary .= "$str.\n"
 			if $fail ==1;
 	
 	}

}

# per base GC content
sub perbase_GC_content {
 	my $module_string = shift;
  	my ( $header, $data, $status ) = initialise_module ( $module_string );
  	
  	$columns .= "$status(";
  	
  	# create hash of the data: position as key, mean GC content as value
	my %data_hash;
	foreach ( split "\n", $data ) {
			$data_hash{ $1 } = $2 if $_ =~ /(\d+)\t(\d+\.*\d*)/;
	}
	
	# mean GC content for all positions: mean of all rows
	my $sum;
	foreach ( values %data_hash ) { $sum += $_; };
	my $sum_mean = $sum / scalar ( keys %data_hash );
	$summary .= "The mean GC content for all positions is ";
	$summary .= sprintf "%.2f", $sum_mean;
	$summary .= "%\n";
	$columns .= sprintf "%.2f", $sum_mean;
	$columns .= ")\n";
	
	# collect positions, if any
	my ( @positions_5, @positions_10 );
	foreach ( keys %data_hash ) {
		my $pos = $_;
		if ( $data_hash{ $pos } > $sum_mean + 5 ||
			$data_hash{ $pos } < $sum_mean - 5 ) {
				push ( @positions_5, $pos );
		}
		if ( $data_hash{ $pos } > $sum_mean + 10 ||
			$data_hash{ $pos } < $sum_mean - 10 ) {
				push ( @positions_10, $pos );
		}
	}
	# warn
	my $str = "The GC content deviates more than 5% from the overall mean at position";
	$str .= range_notation ( join " ", @positions_5 );
 	$summary .= "$str (this causes FastQC to issue the 'warn' status).\n" 
 		if $status =~ /warn/i;
	# fail
	$str = "The GC content deviates more than 10% from the overall mean at position";
	$str .= range_notation ( join " ", @positions_5 );
 	$summary .= "$str (this causes FastQC to issue the 'warn' status).\n"
 		if $status =~ /fail/i;
}

# per sequence GC content
sub perseq_GC_content {
	# more than [15%|30%] of the reads have a GC% that deviates from the theoretical normal distribution
}

# per base N content
sub perbase_N_content {
 	my $module_string = shift;
  	my ( $header, $data, $status ) = initialise_module ( $module_string );

  	# create hash of the data: position as key, N content (is a percentage!) as value
	my %data_hash;
	foreach ( split "\n", $data ) {
			$data_hash{ $1 } = $2 if $_ =~ /(\d+)\t(\d+\.*\d*)/;
	}
	
	# sort the values of the hash
	my @keys_sorted_by_values = sort { $data_hash{$b} <=> $data_hash{$a} } keys %data_hash;
	
	# find maximum value
	$summary .= "The highest observed N content is ";
	$summary .= sprintf "%.2f", $data_hash{ $keys_sorted_by_values[0] };
	$summary .= "%, the lowest observed N content is ";
	$summary .= sprintf "%.2f", $data_hash{ $keys_sorted_by_values[-1] };
	$summary .= "%\n";
	
	# warn: An N content greater than 5% is found at positions ... (generates warn)
	# fail: 20%
	my ( @positions_5, @positions_20 );
	foreach ( keys %data_hash ) {
		push ( @positions_5, $_ ) if $data_hash{ $_ } > 5;
		push ( @positions_20, $_ ) if $data_hash{ $_ } > 20;
	}
	
	# warn
	my $str = "The N content is greater than 5% at position";
	$str .= range_notation ( join " ", @positions_5 );
 	$summary .= "$str (this causes FastQC to issue the 'warn' status).\n"
 		if $status =~ /warn/i;
	# fail
	$str = "The N content is greater than 20% at position";
	$str .= range_notation ( join " ", @positions_20 );
 	$summary .= "$str (this causes FastQC to issue the 'warn' status).\n"
 		if $status =~ /fail/i;

}

# sequence length distribution
# will only print in case of fail or warn
# and in that case, perhaps change the script to process those results? For now I have no example of a warn or fail.
sub seq_length_distr {
 	my $module_string = shift;
  	my ( $header, $data, $status ) = initialise_module ( $module_string, 1 );

	# I picked 'header' because it is the first. Doesn't matter.
	# if it's not defined, the module status is 'pass' and there's nothing to print.
	if ( $header ) {
		$summary .= "Not all sequences have the same length. This causes FastQC to issue the 'warn' status. Please check the full FastQC report for more details.\n" if $status =~ /warn/i;
		
		$summary .= "Sequences of length 0 have been detected. This causes FastQC to issue the 'fail' status. Please check the full FastQC report for more details.\n" if $status =~ /fail/i;
	}
}

# sequence duplication levels
sub seq_duplication_levels {
  	my $module_string = shift;
   	my ( $header, $data, $status ) = initialise_module ( $module_string );

	# ugly but save. Extract the percentage from the first line of the header.
	my $dupl_level = ( split "\t", ( split "\n", $header )[0] )[1];
	$summary .= sprintf "%.2f", $dupl_level;
	$summary .= "% (" . ( $dupl_level / 100 ) * 200000;
	$summary .= " reads) of the ";
	if ( !( $total_num_of_seqs < 200000 ) ) { $summary .= "first 200,000 "; }
	$summary .= "reads are not unique.\n";
}

# overrepresented sequences.
# no data if pass!
# if fail or warn, please adapt this script. I have no example.
sub overrepresented_seqs {
  	my $module_string = shift;
   	my ( $header, $data, $status ) = initialise_module ( $module_string );

	if ( $status =~ /pass/i ) {
		$summary .= "No overrepresented sequences have been detected in the first 200,000 reads.\n";
	} else {
		$summary .= "Overrepresented sequences have been detected in the first 200,000 reads. This causes FastQC to issue the 'warn' or 'fail' status. Please check the full FastQC report for more details.\n";
	}
}

# kmer content. Leave out for now.
sub kmer_content {
#  	my $module_string = shift;
#   my ( $header, $data, $status ) = initialise_module ( $module_string );
}

### other subs

# create a range notation 60, 62-90, 101
# assumes sorted range!
# numbers have to be in space separated string format (1 2 3 4 6 8)
sub range_notation {
	my @numbers = split " ", shift;
	my $range = '';
	
	# check if empty
# 	if ( !( @numbers ) ) { return " -no positions-"; };
	if ( !( @numbers ) ) { return 0; };
	
	# sort the positions first
	my @sorted = sort { $a <=> $b } @numbers;
	
	# if there is more than one number, create some complicated range notation
	if ( scalar ( @sorted ) > 1 ) {
		$range .= 's ';
		my $prev = shift @sorted;
		$range .= $prev;
		foreach my $curr ( @sorted ) {
			if ( ( $curr - $prev ) == 1 ) {
				$range .= "-$curr" if $range !~ /-(\d*)$/;
				$range =~ s/-$prev/-$curr/;
			} else {
				$range .= ", $curr";
			}
			$prev = $curr;
		}
	} else {
		# ... then there's just one number to print.
		$range .= " " . $sorted[0];
	}
	return $range;
}

# short for logging
sub speak {
	my $mess = shift;
	print $lfh "$mess\n";
	# add stuff to print at different debug levels: speak ( WARN, "message" )
}

sub usage {
	print <<EOF;
#########################################################################################
FastQC writes its report to a zip file and, if -Dfastqc.unzip=false was not used, to a
directory as well. The name of the archive and the directory is the same as the FastQC
input reads file with its extension stripped and _fastqc appended. This script takes
that report and translates the data file (always /path/to/report_fastqc/fastqc_data.txt)
into an easy to interpret text file for a quick QC assessment.
#########################################################################################
Usage: ./fastqc_report_v1.pl
\t-r /path/to/sequence.fastq\t[mandatory] fastq file name that served as input
\t\t\t\t\tfor the FastQC report (leave path out)
\t-p /path/to/FastQC_archive/\t[mandatory] path to the directory the FastQC
\t\t\t\t\tarchive lives in
\t-o /path/to/summary.txt\t\t[mandatory] file to write the summary to (append
\t\t\t\t\tis false!)
\t-l /path/to/fastqc_report.log\t[optional] file to write a verbose log to. If
\t\t\t\t\tgiven, verbose output is written to this file. If
\t\t\t\t\tnot given, only exit status is printed to STDERR
#########################################################################################
EOF
 
}