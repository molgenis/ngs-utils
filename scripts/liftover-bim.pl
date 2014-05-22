#!/usr/bin/perl -d:ptkdb
use strict;
use warnings;

die("This tool lifts over a bim file using a VCF file to determine the new positions.\nFiles format:\n* IDs have to match in both files (only matching IDs are kept)\n* All IDs not found in BOTH files will be dropped\nIMPORTANT1: The algorithm used is not memory/cpu optimized and therefore it will work fine on reasonably big SNP files (up to a few million SNPs). Using this on larger number of SNPs can eventually crash Perl.\nIMPORTANT2: If SNPs are dropped, the BIM file will not correspond to its previous FAM/BED files. Best is to then use the result BIM as a SNP filter in PLINK (plink --extract) and use the resulting FAM/BED files with the lifted-over BIM.\nIMPORTANT3: The resulting file is sorted. SNPs which are not sorted after the liftover are dropped altogether as to permit direct use of the FAM/BED files with the BIM.\nUsage: liftover_bim  <in_bim_file> <in.vcf_file> <out_bim>\n") if($#ARGV<2);

if (!(-e $ARGV[1])){ die ("File $ARGV[1] does not exist.") }

open(BIM, "<", $ARGV[0]) or die ("Could not open .bim file $ARGV[0]");
open(VCF, "<", $ARGV[1]) or die ("Could not open VCF file $ARGV[1]");
open(OUT, ">", $ARGV[2]) or die ("Could not create or open target file $ARGV[2]");

print "Starting liftover. Note that this assumes both input files are sorted by chromosome and position.\n";

my $bim;
my $vcf;
my @bim_fields; 
my @previous_bim_fields = (0,0,0,0);
my @vcf_fields;
my %vcf_snps;
my %ambig_snps;
my $snps_not_found = 0;
my $snps_not_ordered = 0;
my $snps_ambig = 0;
my $allele_mismatch = 0;

#Subroutine to allow for non-numeric chr
sub getChrNum{
	my $ch = shift;
	return $ch if($ch=~/\d+/);
	return 23 if lc($ch) eq "x";
	return 24 if lc($ch) eq "y";
	return 25 if lc($ch) eq "xy";
	return 26 if lc($ch) eq "mt";
	return 27; 
}

#Create a big hash with all the VCF SNPs
while($vcf = <VCF>){

	#Skip header lines
	if(!($vcf=~/^#/)){
		@vcf_fields = split(/\t/,$vcf);	
		#Skip / tag ambiguous SNPs IDs
		if(!(defined $ambig_snps{$vcf_fields[2]})){
			if(defined $vcf_snps{$vcf_fields[2]}){
				$ambig_snps{$vcf_fields[2]} = 1;
				delete $vcf_snps{$vcf_fields[2]}	
			}	
			else{
				$vcf_snps{$vcf_fields[2]} = [$vcf_fields[0],$vcf_fields[1],$vcf_fields[3],$vcf_fields[4]];
			}
		}
	}

}
close(VCF);

#Read the bim file and write to bim out file with lifted-over positions
while($bim = <BIM>){
	my $mismatch=0;
	chomp($bim);	
	@bim_fields = split(/\t/,$bim);
	
	if( defined $vcf_snps{$bim_fields[1]}){
		$bim_fields[0] = getChrNum($vcf_snps{$bim_fields[1]}[0]);
		$bim_fields[3] = $vcf_snps{$bim_fields[1]}[1];
		#Check the SNP ordering has not changed
		if((getChrNum($bim_fields[0]) == getChrNum($previous_bim_fields[0]) && $bim_fields[3] > $previous_bim_fields[3]) || getChrNum($bim_fields[0]) > getChrNum($previous_bim_fields[0])){
			#Check that the alleles did not change and add ref allele where necessary
			if($bim_fields[4] eq "0"){
				if($vcf_snps{$bim_fields[1]}[2] eq $bim_fields[5]){
					if($vcf_snps{$bim_fields[1]}[3] ne "."){
						$bim_fields[4] = $vcf_snps{$bim_fields[1]}[3];	
					}
				}
				elsif($vcf_snps{$bim_fields[1]}[3] eq $bim_fields[5]){
					$bim_fields[4] = $vcf_snps{$bim_fields[1]}[2];	
				}
				else{
					$allele_mismatch++;				
					$mismatch++;
				}
					
			}
			elsif($bim_fields[5] eq "0"){
				if($vcf_snps{$bim_fields[1]}[2] eq $bim_fields[4]){
					if($vcf_snps{$bim_fields[1]}[3] ne "."){
						$bim_fields[5] = $vcf_snps{$bim_fields[1]}[3];	
					}
				}
				elsif($vcf_snps{$bim_fields[1]}[3] eq $bim_fields[4]){
					$bim_fields[5] = $vcf_snps{$bim_fields[1]}[2];	
				}
				else{
					$allele_mismatch++;				
					$mismatch++;
				}

			}
			elsif(!($bim_fields[4] eq $vcf_snps{$bim_fields[1]}[2] && $bim_fields[5] eq $vcf_snps{$bim_fields[1]}[3]) && !($bim_fields[4] eq $vcf_snps{$bim_fields[1]}[3] && $bim_fields[5] eq $vcf_snps{$bim_fields[1]}[2])){
				$mismatch++;
				$allele_mismatch++;
			}
			print OUT join("\t",@bim_fields)."\n" if $mismatch == 0;
			@previous_bim_fields = @bim_fields;
		}
		else{
			$snps_not_ordered++;
		}
	}
	elsif(defined $ambig_snps{$bim_fields[1]}){
		$snps_ambig++;
	}
	else{
		$snps_not_found++;
	}			
}

if($snps_not_found){
	print "CAUTION: $snps_not_found SNPs not found in the VCF file. These were removed from the result BIM file.\n";
}
if($snps_not_ordered){
	print "CAUTION: $snps_not_ordered SNPs were excluded since they would have disturbed the ordering.\n";
}
if($snps_ambig){
	print "CAUTION: $snps_ambig SNPs were excluded since they were ambiguous (duplicate ID in VCF file).\n";
}
if($allele_mismatch){
	print "CAUTION: $allele_mismatch SNPs were excluded since there alleles couldn't be matched between the BIM and VCF files\n";
}

print "liftover completed successfully.\n";

close(BIM);
close(OUT);


