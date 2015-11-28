#!/usr/bin/perl
use strict;
use warnings;

die("Usage: reannotate_immuno <in.immunochip_annotation_file> <in_out bim_file>\n") if($#ARGV<1);

if (!(-e $ARGV[1])){ die ("File $ARGV[1] does not exist.") }

#open(BIM, "<", $ARGV[0]) or die ("Could not open .vcf file $ARGV[0]");
open(ANNOTATIONS, "<", $ARGV[0]) or die ("Could not open immunochip annotation file $ARGV[1]");

my @annotations;
my $old;
my $new;
while(<ANNOTATIONS>){
	@annotations=split(/\t/,$_);
	$old = ($annotations[2]);
	$new = ($annotations[3]);
	chomp($old);
	chomp($new);
	if($new ne "NA" && $old ne $new){
	#	print $old." ".$new."\n";
		system("sed -i 's/$old/$new/g' $ARGV[1]");
		if(($? >> 8) != 0){
			die("ERROR while executing sed: $!");
		}
	}
}


#open(VCF_OUT, ">", $ARGV[2]) or die ("Could not create or open target file $ARGV[2]");
#
#print "Starting annotation update. Note that this assumes both input files are sorted by chromosome and position.\n";
#
#my $vcf;
#my $annotation = <ANNOTATIONS>; #Skip header line
#my @vcf_fields;
#my @annotation_fields;
#
#my %from_plink_sex_chr = (23=>"X", 24=>"Y", 25=>"XY", 26=>"MT");
#my %to_plink_sex_chr   = ("X"=>23, "Y"=>24, "XY"=>25, "MT"=>26);
#
#while($vcf = <VCF>){
#	if ($vcf =~ /^[#@]/) {
#		print VCF_OUT $vcf;
#		next;
#	}	
#	elsif($vcf =~ /^(\S+)\s(\d+)\s(\S+)(.+)/){
#		if(!eof ANNOTATIONS){
#		do{
#			$annotation = <ANNOTATIONS>;
#			@annotation_fields = split(/\t/,$annotation);
#			chomp($annotation_fields[3]);									
##			print "VCF: $1 $2 $3 -- IMM: ".join(" ",@annotation_fields)."\n";
#
#			($annotation_fields[0] = $to_plink_sex_chr{$annotation_fields[0]}) if($to_plink_sex_chr{$annotation_fields[0]});
#			
#			if($3 eq $annotation_fields[2] && $annotation_fields[3] ne "NA"){			
#				$vcf = "$1\t$2\t$annotation_fields[3]$4\n";
#				if($annotation_fields[0] != $1 || $annotation_fields[1] != $2){
#					print "WARNING: SNP id=$annotation_fields[3]: VCF coords=(Ch$1,$2) - Immuno coords=(Ch$annotation_fields[0],$annotation_fields[1])\n";
#				}
#				
#			}
#		}
# 		while(!eof ANNOTATIONS && $annotation_fields[0] <= $1 && $annotation_fields[1] < $2);
#		}
#
#		if($vcf =~ /(\S+)\s(.+)/){
#			$vcf = "$from_plink_sex_chr{$1}\t$2\n" if($from_plink_sex_chr{$1});
#		}
#
#		print VCF_OUT $vcf;
#	}
#
#}
#print "Annotation update completed successfully.\n";
#
#
#close(VCF);
close(ANNOTATIONS);
#close(VCF_OUT);
#
#
