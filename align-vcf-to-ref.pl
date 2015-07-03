#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;


print "Important note: the tab-delimited FASTA file refers to UCSC tab-delimited file format coming out of fastaFromBED.\n";
my $usagetxt = "Usage: align-vcf-to-ref.pl [-f] <in.vcf> <in.tab-delimited-fasta> <out.vcf>\nOptions:\n-f: Flips the strand. Note that this will discard all A/T and C/G SNPs as they are ambiguous.\n";

#parse options
my %options;
die($usagetxt) if !getopts("fb", \%options);

die($usagetxt) if($#ARGV<2);

open(VCF, "<", $ARGV[0]) or die ("Could not open input vcf file $ARGV[0]");
open(FASTA, "<", $ARGV[1]) or die ("Could not open input fasta-tab file $ARGV[1]");
open(VCF_OUT, ">", $ARGV[2]) or die ("Could not create or open target vcf file $ARGV[2]");


print "Starting vcf alleles alignment with reference. Note that this assumes both input files are sorted by chromosome and position.\n";

#Initialize vars
my $vcf;
my $fasta;
my @vcf_fields;
my @fasta_fields = (0,0,0,0); #Initialized to be at virtual Ch0 pos 0.
my $vcf_ln = 0;
my $fasta_ln = 0;

my $total_unchanged = 0;
my $total_changed = 0;
my $total_warnings = 0;
my $total_fasta_warnings = 0;
my $total_errors = 0;
my $total_strand_flipped = 0;
my $total_ambiguous = 0;

#Flipped Strand conversion
my %strandFlip = ('A'=>'T','T'=>'A','C'=>'G','G'=>'C','.'=>'.');

#Temp vars used in the loop
my $match_found;
my $i;
my $temp;


#Loop over input VCF file
#The inner loop over the FASTA file is bound by position
while($vcf = <VCF>){
	$match_found = 1;
	$vcf_ln += 1;
	chomp($vcf);
	#Skip header lines
	if ($vcf =~ /^[#@]/) {
		print VCF_OUT "$vcf\n";
		next;
	}	
	#Parse nicely formatted lines
	elsif($vcf =~ /^(\S+)\s(\d+)\s(\S+)\s([ACGT0.])\s([ACGT0.])(.+)/){
	
		#Check VCF field number
		@vcf_fields = split(/\t/, $vcf);
		if($#vcf_fields < 9){
			print "ERROR: Line $vcf_ln of the VCF file does not contain enough fields (only $#vcf_fields found where at least 9 are expected). Line skipped.\n";
			$total_errors++;
			next;
		}
		
		#Recode 0's
		$vcf_fields[3] = "." if($vcf_fields[3] eq "0");
		$vcf_fields[4] = "." if($vcf_fields[4] eq "0");

		#Loop over FASTA file
		while(!eof FASTA && ((getChrNum($fasta_fields[0]) == getChrNum($vcf_fields[0]) && $fasta_fields[2] < $vcf_fields[1]) || getChrNum($fasta_fields[0]) < getChrNum($vcf_fields[0]))){
			$fasta = <FASTA>;
			$fasta_ln += 1;
			#Check and parse FASTA 
			if($fasta =~ /^(\S+):(\d+)-(\d+)\t([ACTG])/){
				$fasta_fields[0] = $1;
				$fasta_fields[1] = $2;
				$fasta_fields[2] = $3;
				$fasta_fields[3] = uc($4);	
			}	
			else{
				print "WARNING: Line $fasta_ln of the FASTA file is not formated correctly. Line skipped.\n";
				$total_fasta_warnings++;
			}
		}	
		#Check that the position was found								
		if(!($fasta_fields[2] eq $vcf_fields[1])){		
			print "WARNING: Position $vcf_fields[0] $vcf_fields[1] in VCF file at line $vcf_ln not found in FASTA file. Line unchanged. \n";
			$total_warnings++;
		}
		#If the input file is not only on + strand, then verify that SNPs are not ambiguous
		elsif(defined $options{"f"} && ($vcf_fields[3] eq $strandFlip{$vcf_fields[4]})){
			$total_ambiguous++;
			$match_found =0;
		}
		#Flip alleles where necessary
		elsif($fasta_fields[3] ne $vcf_fields[3]){
			#Check that reverse alleles match
			if($fasta_fields[3] eq $vcf_fields[4]){
				
				#Flip alleles
				$temp = $vcf_fields[3];
				$vcf_fields[3] = $vcf_fields[4];
				$vcf_fields[4] = $temp;
				#Modify homozygous individuals appropriately
				for($i = 8; $i<=$#vcf_fields; $i++){
					if($vcf_fields[$i]=~/(.*)0\/0(.*)/){
						$vcf_fields[$i] = $1."1/1".$2;
					}
					elsif($vcf_fields[$i]=~/(.*)1\/1(.*)/){
						$vcf_fields[$i] = $1."0/0".$2;
					}
				}
				$vcf = join("\t",@vcf_fields);

				$total_changed++;
			}
			#If the option to flip strand is activated, try flipping strand
			elsif(defined $options{"f"}){
				#Case strand needs to be flipped
				if($strandFlip{$vcf_fields[3]} eq $fasta_fields[3]){
					$vcf_fields[3] = $strandFlip{$vcf_fields[3]};
					$vcf_fields[4] = $strandFlip{$vcf_fields[4]};
					$vcf = join("\t",@vcf_fields);
					$total_strand_flipped++;
				}
				#Case strand and alleles need to be flipped
				elsif($strandFlip{$vcf_fields[4]} eq $fasta_fields[3]){
					
					#Flip alleles
					$temp = $strandFlip{$vcf_fields[3]};
					$vcf_fields[3] = $strandFlip{$vcf_fields[4]};
					$vcf_fields[4] = $temp;
					#Modify homozygous individuals appropriately
					for($i = 8; $i<=$#vcf_fields; $i++){
						if($vcf_fields[$i]=~/(.*)0\/0(.*)/){
							$vcf_fields[$i] = $1."1/1".$2;
						}
						elsif($vcf_fields[$i]=~/(.*)1\/1(.*)/){
							$vcf_fields[$i] = $1."0/0".$2;
						}
					}
					$vcf = join("\t",@vcf_fields);
					$total_strand_flipped++;
				}
				else{
					$match_found = 0;			
					print "WARNING: No allele match for SNP at position $vcf_fields[0] $vcf_fields[1] at line $vcf_ln. Line skipped.\n";
					$total_warnings++;
				}
			}
			#Report non-matching alleles
			else{
				$match_found = 0;
				print "WARNING: No allele match for SNP at position $vcf_fields[0] $vcf_fields[1] at line $vcf_ln. Line skipped.\n";
				$total_warnings++;
			}			
		}
		else{
			$total_unchanged++;
		}
		if($match_found){	
			print VCF_OUT $vcf . "\n";
		}
	}
	#Report ill formated lines
	else{
		print "ERROR: Line $vcf_ln of the VCF is not formated correctly. Line skipped.\n";
		$total_errors++;
	}

}

print "VCF file update completed.\n\n$total_unchanged SNPs did match the reference previously.\n$total_changed SNPs updated.\n";
print "$total_strand_flipped SNPs were strand-flipped\n$total_ambiguous SNPs were ambiguous (A/T, C/G), could not be strand-flipped and were skipped.\n" if(defined $options{'f'});
print "$total_warnings SNPs did not match the reference or were lacking in the reference and could not be updated.\n$total_errors SNPs skipped due to incorrect formatting\n$total_fasta_warnings FASTA lines skipped\n";

close(VCF);
close(FASTA);
close(VCF_OUT);


#### SUBROUTINES ###


#Returns a number for a given chrom incl. X,Y,XY and MT.
#Arguments:
#[0]String: Chrom ID
#Returns:
#[0]int: Chrom Number
sub getChrNum{
	my $ch = shift;
	return $ch if($ch=~/^\d+$/);
	return 23 if lc($ch) eq "x";
	return 24 if lc($ch) eq "y";
	return 25 if lc($ch) eq "xy";
	return 26 if lc($ch) eq "mt";
	return 27; 
}

