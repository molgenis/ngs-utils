#!/bin/bash

set -e
set -u


function showHelp() {
	#
	# Display commandline help on STDOUT.
	#
	cat <<EOH
===============================================================================================================
Script to compare a vcf file with the VKGL data
Usage:
	$(basename "${0}") OPTIONS
Options:
	-h   Show this help.
	-i   inputFile
	-t   workDir working directory (default= CURRENTDIR)
	-b   which build, ucsc (chr1, chr2,chrX etc) or regular (1,2,X etc) default is regular
	-o   outputFolder (default:\${workDir}/output/)

===============================================================================================================
EOH
	trap - EXIT
	exit 0
}


function giabvshc(){
	giabSample="${1}"
	sampleName=$(basename "${giabSample}")
	sampleName=${sampleName%%.*}

	if [ "${buildType}" == "regular" ]
	then
		namingConvention=""
		refGenome="/apps/data/1000G/phase1/human_g1k_v37_phiX.fasta"
		bedfile="/apps/data/Agilent/Exoom_v1/human_g1k_v37/captured.merged.bed"
	elif [ "${buildType}" == "ucsc" ]
	then
		namingConvention="ucsc"
		refGenome="/apps/data/GRC/GCA_000001405.14_GRCh37.p13_full_analysis_set.fna"
		bedfile="/apps/data/Agilent/Exoom_v1-ucsc/human_g1k_v37/captured.merged.bed"
	else
		echo "buildType ${buildType} unknown"
	fi
		
	##capture union bedfile
	bedtools intersect \
		-header \
		-a "${giabSample}" \
		-b "/apps/data/NIST/${namingConvention}/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.19_2mindatasets_5minYesNoRatio_noMT.bed" \
		> "${outputFolder}/tmp/${sampleName}.union.vcf"

	#split giabFile into INDEL and SNP
	# shellcheck disable=SC2154	
	gatk --java-options "-XX:ParallelGCThreads=1 -Xmx5g" SelectVariants \
	-R "${refGenome}" \
	--variant "${outputFolder}/tmp/${sampleName}.union.vcf" \
	-O "${outputFolder}/tmp/${sampleName}_INDEL.union.vcf" \
	--select-type-to-include INDEL \
	-sn "${sampleName}"

	echo "${sampleName}_INDEL.vcf done"

	#Select SNPs and MNPs
	# shellcheck disable=SC2154
	gatk --java-options "-XX:ParallelGCThreads=1 -Xmx5g" SelectVariants \
	-R "${refGenome}" \
	--variant "${outputFolder}/tmp/${sampleName}.union.vcf" \
	-O  "${outputFolder}/tmp/${sampleName}_SNP.union.vcf" \
	--select-type-to-exclude INDEL \
	-sn "${sampleName}"

	echo "${sampleName}_SNP.vcf done"
	
	## capture 20bp bedfile SNP
	bedtools intersect \
		-header \
		-a "${outputFolder}/tmp/${sampleName}_SNP.union.vcf" \
		-b "${bedfile}" \
		> "${outputFolder}/tmp/${sampleName}_SNP.union.20bp.vcf"
	
	## capture 20bp bedfile INDEL
	bedtools intersect \
		-header \
		-a "${outputFolder}/tmp/${sampleName}_INDEL.union.vcf" \
		-b "${bedfile}" \
		> "${outputFolder}/tmp/${sampleName}_INDEL.union.20bp.vcf"

	## Do comparison per type
	firstLine="false"
	for type in SNP INDEL
	do
		grep "^#" "${outputFolder}/tmp/${sampleName}_${type}.union.20bp.vcf" > "${outputFolder}/tmp/${sampleName}_${type}.union.20bp_filtered.vcf"
		grep -v "^#" "${outputFolder}/tmp/${sampleName}_${type}.union.20bp.vcf" |  grep -v "0/0" | grep -v "\./\." | grep "PASS" | grep -v "1/0" >> "${outputFolder}/tmp/${sampleName}_${type}.union.20bp_filtered.vcf"
		echo "comparing ${type}: ${sampleName}vsHC"
		## sample vs HC callset (precision)
		# shellcheck disable=SC2154
		vcf-compare-precision-sensitivity.sh \
		-1 "${outputFolder}/tmp/${sampleName}_${type}.union.20bp_filtered.vcf" \
		-2 "/apps/data/NIST/${namingConvention}/GIAB_HC.${type}_20bp.vcf" \
		-o "${outputFolder}/tmp/precision-${sampleName}vsHC_union_${type}/"


		if [[ "${firstLine}" == "false" ]]
		then
			awk '{if(NR==1){print "Measurement\tType\t"$0}}' "${outputFolder}/tmp/precision-${sampleName}vsHC_union.${type}/comparison.txt" > "${outputFolder}/precision_output.txt"
		fi
		awk -v type=${type} '{if(NR>1){print "precision\t"type"\t"$0}}' "${outputFolder}/tmp/precision-${sampleName}vsHC_union.${type}/comparison.txt" >> "${outputFolder}/precision_output.txt"
		
		##HC callset vs sample (sensitivity)
		# shellcheck disable=SC2154
		vcf-compare-precision-sensitivity.sh \
		-1 "/apps/data/NIST/${namingConvention}/GIAB_HC_${type}_20bp.vcf" \
		-2 "${outputFolder}/tmp/${sampleName}_${type}.union.20bp_filtered.vcf" \
		-o "${outputFolder}/tmp/sensitivity-HCvs${sampleName}_union_${type}/"
		
		if [[ "${firstLine}" == "false" ]]
		then
			awk '{if(NR==1){print "Measurement\tType\t"$0}}' "${outputFolder}/tmp/sensitivity-HCvs${sampleName}_union_${type}/comparison.txt" > "${outputFolder}/sensitivity_output.txt"
			firstLine="True"
		fi
		awk -v type="${type}" '{if(NR>1){print "sensitivity\t"type"\t"$0}}' "${outputFolder}/tmp/sensitivity-HCvs${sampleName}_union_${type}/comparison.txt" >> "${outputFolder}/sensitivity_output.txt"

	done
}

while getopts "i:o:w:b:h" opt; 
do
	# shellcheck disable=SC2249
	# shellcheck disable=SC2220
	case "${opt}" in h)showHelp;; i)inputFile="${OPTARG}";; w)workDir="${OPTARG}";; o)outputFolder="${OPTARG}";;  b)buildType="${OPTARG}";;
esac 
done

if [[ -z "${inputFile:-}" ]]; then showHelp ; echo "inputFile is not specified" ; fi ; echo "inputFile=${inputFile}"
if [[ -z "${workDir:-}" ]]; then workDir="$(pwd)" ; fi ; echo "workDir=${workDir}"
if [[ -z "${outputFolder:-}" ]]; then outputFolder="${workDir}/output/" ; fi ; echo "outputFolder=${outputFolder}"
if [[ -z "${buildType:-}" ]]; then buildType="regular" ; fi ; echo "buildType=${buildType}" 

ml GATK 
ml HTSlib
ml BEDTools
ml ngs-utils

mkdir -p "${workDir}/input/"
rm -rf "${outputFolder}/tmp/"
mkdir -p "${outputFolder}/tmp/"

echo '' > "${outputFolder}/output.txt"

echo "starting giab hc callset"
mapfile -t giabSample < <(find "${inputFile}")
if [[ "${#giabSample[@]:-0}" -eq '0' ]]
then
	echo "${inputFile} not found" 
	echo "${inputFile} not found" >> "${outputFolder}/output.txt"
else
	echo -e "##OPEN## GIAB vs HC and versa" >> "${outputFolder}/output.txt"
	for i in "${giabSample[@]}"
	do
		echo "processing ${i}.."
		giabvshc "${i}"
	done
fi
	
##precision
head -1 "${outputFolder}/precision_output.txt" > "${outputFolder}/output.txt"

#awk '{if (NR>1){print $1,$2"\t"$5"\t"$4"\t"$6}}' "${outputFolder}/precision_output.txt"
awk '{if (NR>1){print $0}}' "${outputFolder}/precision_output.txt" >> "${outputFolder}/output.txt"
cat "${outputFolder}/precision_output.txt" 
echo -e "\n"  >> "${outputFolder}/output.txt" 

cat "${outputFolder}/sensitivity_output.txt"
awk '{if (NR>1){print $0}}' "${outputFolder}/sensitivity_output.txt" >> "${outputFolder}/output.txt"
