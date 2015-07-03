#!/usr/bin/perl


use strict;
use warnings;

my $usagetxt = "Usage: select_snps <in.bim> \"files_pattern\" \"sample_in_pattern\" \"files_out_pattern\"\nOutputs all the snp positions that are in the bim file AND appear in both files the vcftools.diff.sites files returned by the pattern.\nArguments:\nfiles_pattern: bash pattern to select all files (i.e. what you would use with the ls unix command). ex: \"/data/*.diff.sites\"\nsample_in_pattern: a perl regexp that is applied to filenames. It should contain a capturing group that is the sample unique name. ex: \"(A\\d+[abc])\\.\" could be used for filename like A1a.diff.sites, A2b.diff.sites, A1c.diff.sites, etc.\nfiles_out_pattern: A string for the outfilenames where \$1 will be replaced by the corresponding sample name. ex: \"\$1.diff.sites.list\" will give you files like sample1.diff.sites.list, sample2.diff.sites.list\n";

die $usagetxt if $#ARGV < 3;

#list files to be processed
my $filename;
my @infiles;
my $infile;
foreach $filename(<$ARGV[1]>){
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


die "No input file found. Operation canceled.\n" if !@infiles;


#Read bim file
my %snplist;
my $use_snplist = 0;
my $snpfile;
my @line;
die("Could not open bm file: $ARGV[0]\n") if !open($snpfile, "<", $ARGV[0]);
while(<$snpfile>){
	@line = split(/\t/,$_);
	if($#line > 2){
		#$snplist{getChrName($line[0])."-".$line[3]} = 0;
		$snplist{getChrName($line[0])."-".$line[3]} = $line[1];
	}
}
close($snpfile);

#Read vcftools.diff.sites files
my $outfile;
my $samplename;
foreach $filename(@infiles){
	next if(!open($infile, "<", $filename));
	$filename =~ $ARGV[2];
	$samplename = $1;
	$filename = $ARGV[3];
	$filename =~ s/X1/$samplename/g;
	next if(!open($outfile,">",$filename));	
	
	while (<$infile>){
		@line = split(/\t/,$_);		
		if($line[2] eq "B" && defined($snplist{$line[0]."-".$line[1]})){
			#$snplist{$line[0]."-".$line[1]} = 1;
			print $outfile $line[0]."\t".$line[1]."\t".$snplist{$line[0]."-".$line[1]}."\n";
		}	

	}

	close($infile);	
	close($outfile);
}

##Sort results
#foreach my $snp (keys %snplist){		
#	delete $snplist{$snp} if($snplist{$snp} < 1);
#}	
#my @keys = sort{chromPosSort($a,$b)}keys %snplist;
#
##Print results
#my $snp;
#foreach $snp (@keys){
#	@line = split(/-/,$snp);
#	print join("\t",@line)."\n";
#}


###SUBS

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

