#!/usr/bin/perl
use strict;
use warnings;

die("Usage: flip-vcf-snp.pl <in.vcf> <in.bim_flipped_file> <out.vcf>\n") if($#ARGV<2);

open(VCF, "<", $ARGV[0]) or die ("Could not open .vcf file $ARGV[0]");
open(BIM, "<", $ARGV[1]) or die ("Could not open bim file $ARGV[1]");
open(VCF_OUT, ">", $ARGV[2]) or die ("Could not create or open target file $ARGV[2]");

print "Starting flip. Note that this assumes both input files are sorted by chromosome and position.\n";

#Reverse complement
my %reverse_complement = ('A'=>'T', 'T'=>'A', 'C'=>'G','G'=>'C','.'=>'.');

my $vcf;
my $bim;
my @vcf_fields;
my @bim_fields = (0,0,0,0); #Initialized to be at virtual Ch0 pos 0.
my $num_flipped = 0;
my $num_not_flipped = 0;
my $num_no_match = 0;
my $num_not_flipped_alleles_freq_change = 0;

while($vcf = <VCF>){
	if ($vcf =~ /^[#@]/) {
		print VCF_OUT $vcf;
		next;
	}	
	elsif($vcf =~ /^(\S+)\s(\d+)\s(\S+)\s([ACGT.])\s([ACGT.])(.+)/){
		
		while(!eof BIM && (($bim_fields[0] == $1 && $bim_fields[3] < $2) || $bim_fields[0] < $1)){
			$bim = <BIM>;
			@bim_fields = split(/\t/,$bim);	
		}	
			chomp($bim_fields[5]);									
			if(!($3 eq $bim_fields[1])){		
				print "WARNING: SNP id Mismatch. VCFid=$3 - Immunoid=$bim_fields[1]\n";
				#next;
			}

			#Format missing information from BIM to VCF format
			if($bim_fields[4] eq '0'){
				$bim_fields[4] = '.';
			}
			if($bim_fields[5] eq '0'){
				$bim_fields[5] = '.';
			}

			#Check if flip needed
			#Complete match
			if(($bim_fields[5] eq $4 || $bim_fields[5] eq '.' ) && ($bim_fields[4] eq $5 || $bim_fields[4] eq '.' )){
				$num_not_flipped+=1;
			}
			else{
				
				#Check that alleles match before flipping
				if(($reverse_complement{$4} eq $bim_fields[5] && $reverse_complement{$5} eq $bim_fields[4]) ||
				   ($reverse_complement{$4} eq $bim_fields[4] && $reverse_complement{$5} eq $bim_fields[5]) ||
				   ($reverse_complement{$4} eq $bim_fields[5] && ($bim_fields[4] eq '.' || $5 eq '.' )) ||
				   ($reverse_complement{$4} eq $bim_fields[5] && ($bim_fields[4] eq '.' || $5 eq '.' ))){
					$num_flipped+=1;
					$vcf = "$1\t$2\t$3\t$reverse_complement{$4}\t$reverse_complement{$5}$6\t\n";
					
				}
				#No match, but perfect reverse match (not flipped)
				elsif($4 eq $bim_fields[4] && $5 eq $bim_fields[5]){
					$num_not_flipped_alleles_freq_change+=1;
				}
				#No match
				else{
					$num_no_match +=1;
					print "WARNING: SNP $1:$2 alleles did not match. VCF: $4/$5\tBIM:$bim_fields[5]/$bim_fields[4]\n";
				}

				if($bim_fields[0] != $1 || $bim_fields[3] != $2){
					print "WARNING: SNP id=$bim_fields[1]: VCF coords=(Ch$1,$2) - Immuno coords=(Ch$bim_fields[0],$bim_fields[3])\n";
				}
			}	
		
		print VCF_OUT $vcf;
	}

}
print "$num_no_match SNPs did not match, $num_flipped SNPs flipped, $num_not_flipped SNPs already on the same strand, $num_not_flipped_alleles_freq_change SNPs already on the same strand but with different allele freq.\n";
print "Bim strand flipping completed.\n";


close(VCF);
close(BIM);
close(VCF_OUT);


