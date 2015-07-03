#!/usr/bin/perl 
#-d:ptkdb
use strict;
use warnings;
#ARGUMENTS:
use Getopt::Std;

my $usagetxt = "Usage: vcftools_diff_sites_concordance.pl [-options] \"in_diff_files\" <out>\nExample: vcftools_diff_concordance \"./*.diff.sites\" concordance_result\nIMPORTANT: do not forget quotes around the first argument if you want to use wildcards (*).\nOptions:\n-t: outputs the result in a tab-delimited way rather than textual way.\n-n: ignores sex chromosomes\n-l <snps_position_file>: only considers snps which positions are in the file. The file should be tab-delimited with 2 columns: chromosome position. Note that all additional columns are ignored but will not disturb the process (f.ex. a VCF file is a valid input).\n-L \"suffix\": Same as -l but uses a list of SNPs per file. The files have to be named like the sites files and use a static suffix. Note that this is ignored if -s is used. If used in conjunction with -l, union of the files is considered. If a file is missing, it is ignored.\n-p \"pattern\": a perl pattern that will be applied to the filenames to extract the sample name for reporting. The sample name is defined as the first capturing group. Ignored if -s is used.\n-s: output snp-major list instead of sample-major list. Note that output using -s are always displayed in a tab-delimited file. Also by default, it only outputs SHARED SNPs. You can override this by using the -C option.\n-b <bim_file>: Can only be used along with -s option. uses the BIM file provided to output the SNPids along with other SNP information.\n-f <plink_frq_file>: Only valid along with -s AND -b options. Uses the FRQ file provided to output the maf along with other SNP information.\n-C: Output all SNPs. Only valid along with -s option (ignored otherwise). This is both memory and cpu intensive (it sorts everything) and is therefore not recommended.\n-S: Do not sort output file. Used along with -s option. Dramatically reduces memory/cpu/time.\n";

#Disable STDOUT buffer
$|=1;

#parse options
my %options;
die($usagetxt) if !getopts("tnl:sCSb:f:L:p:", \%options);

#Check options dependencies
die($usagetxt) if(defined $options{"b"} && !defined $options{"s"});
die($usagetxt) if(defined $options{"f"} && !defined $options{"b"});

#Check arguments
die($usagetxt) if $#ARGV < 1;

my $filename;
my @infiles;
my $infile;
my @concordance;

#list files to be processed
foreach $filename(<$ARGV[0]>){
	if(-f $filename){
		push @infiles, $filename;
	}	
	elsif(-d $filename){
		print "Directories are not supported. Skipping $filename\n";
	}
	else{
		print "$filename not found. Skipped\n";
	}
}

print "Found ".(@infiles)." input files.\n";

die "No input file found. Operation canceled.\n" if !@infiles;

#Open and parse snp list file if any
my %snplist;
my $use_snplist = 0;
if(defined $options{"l"}){
	my $snpfile;
	die("Could not open snp list file (-l option): $options{'l'}\n") if !open($snpfile, "<", $options{"l"});
	while(<$snpfile>){
		if($_=~/^([1234567890xyXYmtMT]+)\t(\d+)/){
			$snplist{$1."-".$2} = 1;
		}
	}
	$use_snplist=1;
}

#Prepare variables for sample-SNP list 
my %samplesnplist;
my $use_samplesnplist = 0;
my $samplesnpfile;
if(defined $options{"L"}){
	$use_samplesnplist = 1;
}


#Open and parse bim file if provided
my %snps;
my $use_bim = 0;
if(defined $options{"b"}){
	#Open and parse frq file if provided
	my %frq;
	my $use_frq = 0;
	if(defined $options{"f"}){
		print "Reading frequency file..\n";
		my $frqfile;
		my @frqfields;
		die("Could not open FRQ file (-f option): $options{'f'}\n") if !open($frqfile, "<", $options{"f"});
		#skip header line
		<$frqfile>;
		while(<$frqfile>){
			#trim/chomp line
			$_ =~ s/^\s+|\s$//;
			@frqfields = split(/\s+/,$_);
			#Check for duplicates
			if(defined $frq{$frqfields[1]}){
				$frq{$frqfields[1]} = "AMBIG";
			}
			else{
				$frq{$frqfields[1]} = $frqfields[4];
			}
		}
		$use_frq=1;
	}
	
	print "Reading BIM file...\n";
	my $bimfile;
	my @bimfields;
	my $snp_pos;
	die("Could not open BIM file (-b option): $options{'b'}\n") if !open($bimfile, "<", $options{"b"});
	while(<$bimfile>){
		@bimfields = split(/\t/,$_);
		#Check for duplicates; keep rs reference if duplicate
		$snp_pos = getChrName($bimfields[0])."-".$bimfields[3];
		if(!defined $snps{$snp_pos} || $bimfields[1] =~ /rs/){
			$snps{$snp_pos} = {'snpid' => $bimfields[1],'total'=>0, 'only_1' => 0, 'only_2' => 0, 'both'=>0, 'matching_alleles' => 0, 'common_called' => 0, 'discordant' => 0, 'discordance' => 0};
			#Add frq is set to
			if($use_frq && $frq{$bimfields[1]}){
				#push(@{$snps{$snp_pos}},$frq{$bimfields[1]});
				$snps{$snp_pos}->{'maf'} = $frq{$bimfields[1]};
			}	
		}
	}
	$use_bim=1;
}

#foreach my $snp(keys %snps){
#	print $snp.":".$snps{$snp}->{'snpid'}.",".$snps{$snp}->{'maf'}."\n";
#}
#die;

#Starting Concordance aggregation
print "Aggregating concordance files...\n";

#Check if ignore sex chromosomes
my $sex=1;
if ( defined $options{"n"} ) {
      $sex = 0;
}

#Calculate site concordance for each file.
foreach $filename(@infiles){
	if(!open($infile, "<", $filename)){
		print "Could not open file $filename. File skipped.";
		next;
	}
	print "Processing file $filename...\n";
	if(defined $options{"s"}){
		calcPerSNP();
	}
	else{

		if($use_samplesnplist){
			if(open($samplesnpfile, "<", $filename.$options{"L"})){
				%samplesnplist = ();
				while(<$samplesnpfile>){
					if($_ =~ /^([1234576890xXyYmtMT]+)\t(\d+)\t/){
						$samplesnplist{$1."-".$2}=1;
					}
				}	
			}
		}
		#Get the sample name
		if(defined $options{"p"} && $filename =~ $options{"p"}){;
			calcPerSample($1);
		}
		else{
			calcPerSample($filename);
		}
	}
	close($infile);	
}

#Output concordance


my $outfile;
open($outfile, ">", $ARGV[1]) or die("Could not open output file: $ARGV[1]. Process aborted.\n");

if(defined $options{"s"}){
	if(!defined $options{"C"}){
		#Remove all non-common SNPs
		foreach my $snp (keys %snps){		
			delete $snps{$snp} if($snps{$snp}->{'both'} < 1);
		}	
	}	
	if(!defined $options{"S"}){
		print "Sorting results...\n";
		@concordance = sort{chromPosSort($a,$b)}keys %snps;
	}
	else{
		@concordance = keys %snps;
	}

	print "Printing concordance report...\n";
	printPerSNP();
}
else{
	print "Printing concordance report...\n";
	printPerSample();
}

print "Done.\n";

#### SUBROUTINES ###

#Calculates the concordance per sample for 1 file
#ARGUMENTS:
#Filename or sample name
#Alos, @concordance, %snplist, $sex and $use_snplist global variables have to be present
sub calcPerSample{
	
	my @cols;
	my %filesums = ('file'=>shift,'total'=>0, 'only_1' => 0, 'only_2' => 0, 'both'=>0, 'matching_alleles' => 0, 'common_called' => 0, 'discordant' => 0, 'discordance' => 0);
	<$infile>; #skip line 1

	#Sum columns
	while(<$infile>){
		@cols = split(/\t/,$_);
		#If the nosex option is on, skip sex chromosomes
		if($sex || $cols[0]=~/\d+/){
			#if the snp list option is on, skip positions not present in the file
			if(!$use_samplesnplist && !$use_snplist || ($use_samplesnplist && defined $samplesnplist{$cols[0]."-".$cols[1]}) || defined $snplist{$cols[0]."-".$cols[1]}){
				++$filesums{'total'};
				if($cols[2] eq "B"){
					++$filesums{'both'};
					$filesums{'matching_alleles'} += $cols[3];
					$filesums{'common_called'} += $cols[4];
					$filesums{'discordant'} += $cols[5];
					$filesums{'discordance'} += $cols[6];
					
				}
				elsif($cols[2] == 1){
					++$filesums{'only_1'};
				}
				elsif($cols[2] == 2){
					++$filesums{'only_2'};
				}
			}
		}
	}
	
	#save file concordance in general concordance
	push @concordance, {%filesums};
}

#Prints out the sample-based results
#ARGUMENTS
#None but @concordance should be properly formated
sub printPerSample{
	if (defined $options{"t"} ) {
		print $outfile "file\tloci in file 1 only\tloci in file only 1(%)\tloci in file 2 only\tloci in file 2 only(%)\tloci in both files\tloci in both files(%)\tmatching alleles\tmatching alleles(%)\tcommon called\tcommon called(%)\tdiscordant\tdiscordant(%)\tdiscordance\tdiscordant(%)\n";
			
		for my $i (0..$#concordance){
			print $outfile $concordance[$i]{'file'}."\t";
			print $outfile sprintf("%d\t%.3f\t", $concordance[$i]{'only_1'}, ($concordance[$i]{'only_1'}/$concordance[$i]{'total'})*100);
			print $outfile sprintf("%d\t%.3f\t", $concordance[$i]{'only_2'}, ($concordance[$i]{'only_2'}/$concordance[$i]{'total'})*100);
			print $outfile sprintf("%d\t%.3f\t", $concordance[$i]{'both'}, ($concordance[$i]{'both'}/$concordance[$i]{'total'})*100);
			print $outfile sprintf("%d\t%.3f\t", $concordance[$i]{'matching_alleles'}, ($concordance[$i]{'matching_alleles'}/$concordance[$i]{'both'})*100);
			print $outfile sprintf("%d\t%.3f\t", $concordance[$i]{'common_called'}, ($concordance[$i]{'common_called'}/$concordance[$i]{'both'})*100);
			print $outfile sprintf("%d\t%.3f\t", $concordance[$i]{'discordant'}, ($concordance[$i]{'discordant'}/$concordance[$i]{'both'})*100);
			print $outfile sprintf("%d\t%.3f\n", $concordance[$i]{'discordance'}, ($concordance[$i]{'discordance'}/$concordance[$i]{'both'})*100);
		}
		
	}
	else{
		for my $i (0..$#concordance){
			print $outfile $concordance[$i]{'file'}.":\n";
			print $outfile sprintf("Loci in file 1 only: %d (%.3f%%)\n", $concordance[$i]{'only_1'}, ($concordance[$i]{'only_1'}/$concordance[$i]{'total'})*100);
			print $outfile sprintf("Loci in file 2 only: %d (%.3f%%)\n", $concordance[$i]{'only_2'}, ($concordance[$i]{'only_2'}/$concordance[$i]{'total'})*100);
			print $outfile sprintf("Loci in both files: %d (%.3f%%)\n", $concordance[$i]{'both'}, ($concordance[$i]{'both'}/$concordance[$i]{'total'})*100);
			print $outfile sprintf("\tMatching Alleles: %d (%.3f%%)\n", $concordance[$i]{'matching_alleles'}, ($concordance[$i]{'matching_alleles'}/$concordance[$i]{'both'})*100);
			print $outfile sprintf("\tCommon called: %d (%.3f%%)\n", $concordance[$i]{'common_called'}, ($concordance[$i]{'common_called'}/$concordance[$i]{'both'})*100);
			print $outfile sprintf("\tDiscordant: %d (%.3f%%)\n", $concordance[$i]{'discordant'}, ($concordance[$i]{'discordant'}/$concordance[$i]{'both'})*100);
			print $outfile sprintf("\tDiscordance: %d (%.3f%%)\n", $concordance[$i]{'discordance'}, ($concordance[$i]{'discordance'}/$concordance[$i]{'both'})*100);
			print $outfile "----------------------\n\n";
		}
	}
}

#Calculates the concordance per SNP for 1 file
#ARGUMENTS:
#None, but %bim, %snplist, $sex and $use_snplist and @concordance global variables have to be present
sub calcPerSNP{
	
	my @cols;
	my $snp_pos;
	<$infile>; #skip line 1

	#Sum columns
	while(<$infile>){
		@cols = split(/\t/,$_);
		#If the nosex option is on, skip sex chromosomes
		if($sex || $cols[0]=~/\d+/){
			#if the snp list option is on, skip positions not present in the file
			$snp_pos = $cols[0]."-".$cols[1];
			if(!$use_snplist || defined $snplist{$snp_pos}){
				if(!defined $snps{$snp_pos}){
					$snps{$snp_pos} = {'total'=>0, 'only_1' => 0, 'only_2' => 0, 'both'=>0, 'matching_alleles' => 0, 'common_called' => 0, 'discordant' => 0, 'discordance' => 0};
				}
				$snps{$snp_pos}->{"total"} += 1;
				if($cols[2] eq "B"){
					$snps{$snp_pos}->{'both'} += 1;
					$snps{$snp_pos}->{'matching_alleles'} += $cols[3];
					$snps{$snp_pos}->{'common_called'} += $cols[4];
					$snps{$snp_pos}->{'discordant'} += $cols[5];
					$snps{$snp_pos}->{'discordance'} += $cols[6];
					
				}
				elsif($cols[2] == 1){
					$snps{$snp_pos}->{'only_1'} += 1;
				}
				elsif($cols[2] == 2){
					$snps{$snp_pos}->{'only_2'} += 1;
				}
				#Use @concordance to keep the keys sorted in the same way the input
				#file is
				push(@concordance,$snp_pos);
			}
		}
	}
}

#Prints out the SNP-based results
#ARGUMENTS
#None but @concordance should be properly formated
sub printPerSNP{
	my $optional_cols = "";
	my $out_common = defined($options{"C"});
	if(defined $options{"b"}){
		$optional_cols .= "\tid";
		if(defined $options{"f"}){
			$optional_cols .= "\tmaf";
		}	
	}
	print $outfile "chrom\tposition$optional_cols\tloci in file 1 only\tloci in file only 1(%)\tloci in file 2 only\tloci in file 2 only(%)\tloci in both files\tloci in both files(%)\tmatching alleles\tmatching alleles(%)\tcommon called\tcommon called(%)\tdiscordant\tdiscordant(%)\tdiscordance\tdiscordant(%)\n";
	
	foreach my $snp (@concordance){
		#Check that the format of the position is OK
		if($snp =~ /^(.+)-(.+)$/){
			print $outfile $1."\t".$2."\t";
			#Print optional columns if any
			if(defined $options{"b"}){
				if(defined $snps{$snp}->{'snpid'}){
					print $outfile sprintf("%s\t",$snps{$snp}->{'snpid'});
				}
				else{
					print $outfile "MISSING\t";
				}
				if(defined $options{"f"}){
					if(defined $snps{$snp}->{'maf'}){
						print $outfile sprintf("%s\t", $snps{$snp}->{'maf'});
					}
					else{
						print $outfile "-1\t";
					}
				}
			}
			#Print standard columns
			print $outfile sprintf("%d\t%.3f\t", $snps{$snp}->{'only_1'}, ($snps{$snp}->{'only_1'}/$snps{$snp}->{'total'})*100);
			print $outfile sprintf("%d\t%.3f\t", $snps{$snp}->{'only_2'}, ($snps{$snp}->{'only_2'}/$snps{$snp}->{'total'})*100);
			if($snps{$snp}->{'both'} > 0){
				print $outfile sprintf("%d\t%.3f\t", $snps{$snp}->{'both'}, ($snps{$snp}->{'both'}/$snps{$snp}->{'total'})*100);
				print $outfile sprintf("%d\t%.3f\t", $snps{$snp}->{'matching_alleles'}, ($snps{$snp}->{'matching_alleles'}/$snps{$snp}->{'both'})*100);
				print $outfile sprintf("%d\t%.3f\t", $snps{$snp}->{'common_called'}, ($snps{$snp}->{'common_called'}/$snps{$snp}->{'both'})*100);
				print $outfile sprintf("%d\t%.3f\t", $snps{$snp}->{'discordant'}, ($snps{$snp}->{'discordant'}/$snps{$snp}->{'both'})*100);
				print $outfile sprintf("%d\t%.3f\n", $snps{$snp}->{'discordance'}, ($snps{$snp}->{'discordance'}/$snps{$snp}->{'both'})*100);
			}
			else{
					
				print $outfile sprintf("%d\t0\t", $snps{$snp}->{'both'});
				print $outfile sprintf("%d\t0\t", $snps{$snp}->{'matching_alleles'});
				print $outfile sprintf("%d\t0\t", $snps{$snp}->{'common_called'});
				print $outfile sprintf("%d\t0\t", $snps{$snp}->{'discordant'});
				print $outfile sprintf("%d\t0\n", $snps{$snp}->{'discordance'});
			}
		}		
	}	
}

#Sorting by chromosome/position function based on a string chrom-position
sub chromPosSort{
	my @chrompos1 = split(/-/,$_[0]);
	my @chrompos2 = split(/-/,$_[1]);
	
	#If not in the right format, sort last
	return 1 if($#chrompos1 < 1);
	return -1 if($#chrompos2 < 1);

	#sort by chrom first
	return 1 if(getChrNum($chrompos1[0]) > getChrNum($chrompos2[0]));
	return -1 if(getChrNum($chrompos2[0]) > getChrNum($chrompos1[0]));

	#Then by position
	return 1 if($chrompos1[1] > $chrompos2[1]);
	return -1 if($chrompos2[1] > $chrompos1[1]);

	return 0;
}

#Returns a number for a given chrom incl. X,Y,XY and MT.
sub getChrNum{
	my $ch = shift;
	return $ch if($ch=~/^\d+$/);
	return 23 if lc($ch) eq "x";
	return 24 if lc($ch) eq "y";
	return 25 if lc($ch) eq "xy";
	return 26 if lc($ch) eq "mt";
	return 27; 
}

#Returns X/Y chrom
sub getChrName{
	my $ch = shift;
	return "X" if $ch == 23;
	return "Y" if $ch == 24;
	return "XY" if $ch == 25;
	return "MT" if $ch == 26;
	return $ch;
}

