#!/bin/sh

##
# bed_filter_batch_generator.sh <path to folder containing vcf files>
##
# B Miedema
##
# This script generates bash scripts for filtering vcf files based on locations in canonical_exon_sequences_filtered_headered_sorted.bed file.
# It also bgzips and generates Tabix index files for these vcf files, the original filtered vcf is also maintained.
# 
# Requires a template script located in the same folder as this script. The template script is left intact.
# 
# output: bash scripts to be deployed run for filtering vcf files on locations in canonical_exon_sequences_filtered_headered_sorted.bed bedfilefunction.
# These scripts can be run on the cluster. The output scripts are named in this form: bed_filtered.<VCFfile name>.sh.
##

currentDir=${PWD}
pathToVcfFiles=$1
bedFile=$2

if [ "$1" = "" ]; then 
	echo "This script generates bash scripts for filtering vcf files based on genomic locations specified in supplied BED file."
	echo "It also bgzips and generates Tabix index files for the filtered vcf files, the original filtered vcf is also maintained."
	echo "Requires a template script located in the same folder as this script. The template script is left intact.\n"
	echo "Also requires and installation of VCFtools (at least 1.12b), tabix and bgzip. "
	echo "output: bash scripts with a bed filtering command specified by the commandline option bed file."
	echo "usage: bed_filter_batch_generator.sh <path to folder containing vcf files> <bedfile>"

	exit
fi

cd $pathToVcfFiles
for item in $( ls -1 *.vcf.gz.gz ); do
	cd $currentDir
	
	strippedFileName=${item/\//_}  
	
	BgZipCmd="bgzip -c ./filtered.${strippedFileName}.recode.vcf > ./filtered.${strippedFileName}.recode.vcf.gz"
	TabixCmd="tabix -s 1 -b 2 -e 3 ./filtered.${strippedFileName}.recode.vcf.gz"
	bedFilterCommand="vcftools --gzvcf ${pathToVcfFiles}${item} --bed ${bedFile} --out filtered.${strippedFileName} --recode --remove-indels --recode-INFO-all > ${outputPath}filtered.${item}.log"
	
	sed "s|<line>|${bedFilterCommand}|g" bedfilter.template > bed_filtered.${item}.sh
	sed -i "s|<bgzip>|${BgZipCmd}|g" bed_filtered.${item}.sh
	sed -i "s|<tabix>|${TabixCmd}|g" bed_filtered.${item}.sh
	sed -i "s|<INPUT>|err_log|g"  bed_filtered.${item}.sh
	
	bedFilterComand=""
    done
echo "DONE! wrote output to $outputPath"
	


