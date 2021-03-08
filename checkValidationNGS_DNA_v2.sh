#!/bin/bash

set -e
set -u


function showHelp() {
	#
	# Display commandline help on STDOUT.
	#
	cat <<EOH
===============================================================================================================
Script to copy (sync) data from a succesfully finished analysis project from tmp to prm storage.
Usage:
	$(basename $0) OPTIONS
Options:
	-h   Show this help.
	-i   inputFolder
	-t   inputType (vcf or vcf.gz) (default= vcf.gz)
	-o   outputFolder (default:\${inputFolder}/output/)
	-v   validationFolder, folder where the vcfs are with the SNPs that should be found back (default=/groups/umcg-gd/prm06/projects/validationVcfs/)
===============================================================================================================
EOH
	trap - EXIT
	exit 0
}

function doVariantEval(){

	folder="${validationFolderTmp}"
	output="${outputFolder}"

	for i in $(ls "${folder}/"*".${inputType}")
	do
		name=$(basename $i ".${inputType}")
		java -jar ${EBROOTGATK}/GenomeAnalysisTK.jar \
		-T VariantEval \
		-R /apps/data/1000G/phase1/human_g1k_v37_phiX.fasta \
		-o "${output}/output.${name}.eval.grp" \
		--eval "${inputFolder}/"*"${name}"*".${inputType}" \
		--comp "${i}"

	done

	for i in $(ls "${folder}/"*".${inputType}")
	do

		name=$(basename "${i}" ".${inputType}")
		check=$(awk '{if (NR==5){if ($11 == "100.00"){print "correct"}}}' "${output}/output.${name}.eval.grp")
		if [ "${check}" == "correct" ]
		then
			zcat "${i}" | awk -v sample="${name}" 'BEGIN {OFS="  "}{if ($1 !~ /^#/){print sample,$1,$2,$4,$5,"FOUND BACK"}}'
		else
			zcat "${i}" | awk -v sample="${name}" 'BEGIN {OFS="  "}{if ($1 !~ /^#/){print sample,$1,$2,$4,$5,"Not 100% concordant!"}}' 
			exit 1
		fi
	done

}

function doComparisonFiltered (){
	refCall="${1}"
	if [ "${refCall}" == "referenceCall" ]
        then
		folder="${validationFolderTmp}/filtered/referenceCall"
                output="${outputFolder}/filtered/referenceCall"
        else
		folder="${validationFolderTmp}/filtered/"
                output="${outputFolder}/filtered/"
        fi
	mkdir -p "${output}"

for i in $(ls "${folder}/"*".${inputType}")
do

	name=$(basename "${i}" ".${inputType}")

        if [[ "${inputType}" == "vcf.gz" ]]
        then
		validationSample=$(zcat ${i} | grep -v '^#' | awk '{print $1"-"$2"-"$4"-"$5"-"$7}')
                inputSample=$(zcat "${inputFolder}/"*"${name}"*".${inputType}" | grep 18598089 | awk '{print $1"-"$2"-"$4"-"$5"-"$7}')

                if [ "${validationSample}" == "${inputSample}" ]
                then
			if [ "${refCall}" == "referenceCall" ]
                        then
				zcat "${i}" | awk -v sample="${name}" 'BEGIN {OFS="  "}{if ($1 !~ /^#/){print sample,$1,$2,$4,$5,$7,"FOUND BACK,REF CALL"}}'
                        else
                                zcat "${i}" | awk -v sample="${name}" 'BEGIN {OFS="  "}{if ($1 !~ /^#/){print sample,$1,$2,$4,$5,$7,"FOUND BACK"}}'
                        fi
                else
			zcat "${i}" | awk -v sample="${name}" 'BEGIN {OFS="  "}{if ($1 !~ /^#/){print sample,$1,$2,$4,$5,$7,"Not 100% concordant!"}}' 
                        zcat "${inputFolder}/"*"${name}"*".${inputType}" | grep 18598089
                        zcat "${i}"
                        exit 1
                fi
        fi
done

}

function checkFrankenstein() {

        knownVariants=${validationFolderTmp}/ValidationSet.annotated.vcf.gz
        knownMissingVariants=${validationFolderTmp}/knownMissingVariants.txt

        totalNoOfKnownVariants=$(zcat "${knownVariants}" | grep -v '^#' | wc -l)

        echo "Testing Frankenstein on ${totalNoOfKnownVariants} known variants"

	tmpFolder="${outputFolder}/tmp"
	input="${inputFolder}/Frankenstein_1.final.vcf.gz"

	echo "comparing \"true\" set: /groups/umcg-atd/tmp01/ValidationSet.annotated.vcf.gz"
	echo "with new: ${input}"
	echo "working folder: ${tmpFolder}"

	bedtools intersect -header -b ${input} -a ${knownVariants} > ${tmpFolder}/SameAsFrank.vcf
	bgzip -c ${tmpFolder}/SameAsFrank.vcf > ${tmpFolder}/SameAsFrank.vcf.gz

	totalNoOfKnownVariants=$(zcat "${knownVariants}" | grep -v '^#' | wc -l)
	totalSame=$(cat "${tmpFolder}/SameAsFrank.vcf" | grep -v '^#' | wc -l)

	echo "found ${totalSame} of in total ${totalNoOfKnownVariants} back in the ${input}"

	tortilla.sh -v -1  ${tmpFolder}/SameAsFrank.vcf.gz -2 ${knownVariants} -o ${outputFolder}


	if diff -q "${knownMissingVariants}" ${outputFolder}/notInVcf1.txt
	then
	    	echo "Everything is OK, the missing variants are the same as the already missing variants"
	else
	    	echo "diff ${knownMissingVariants} ${outputFolder}/notInVcf1.txt"
	        echo "OOPS, there are differences between the known missing variants and the sample!"

	        diff -y "${knownMissingVariants}" ${outputFolder}/notInVcf1.txt
	fi
}


while getopts "i:o:v:t:h" opt; 
do
	case $opt in h)showHelp;; i)inputFolder="${OPTARG}";; o)outputFolder="${OPTARG}";; v)validationFolderPrm="${OPTARG}";; t)inputType="${OPTARG}";;
esac 
done

if [[ -z "${inputFolder:-}" ]]; then showHelp ; echo "inputFolder is not specified" ; fi ; echo "inputFolder=${inputFolder}"
if [[ -z "${outputFolder:-}" ]]; then mkdir -p "${inputFolder}/output/tmp" ; outputFolder="${inputFolder}/output/" ; fi ; echo "outputFolder=${outputFolder}"
if [[ -z "${inputType:-}" ]]; then inputType="vcf.gz" ; fi ; echo "inputType=${inputType}"
if [[ -z "${validationFolderPrm:-}" ]]; then validationFolderPrm="/groups/umcg-gd/prm06/projects/validationVcfs/" ; fi ; echo "validationFolderPrm=${validationFolderPrm}"

ml GATK 
ml HTSlib
ml BEDTools
ml ngs-utils

mkdir -p "${outputFolder}"


if [ $(hostname) != "calculon" ]
then
	validationFolderTmp=${inputFolder}/validationVcfs/
	mkdir -p "${outputFolder}/filtered/"

	echo "copying validationVcfs"
	if [ -f "${validationFolderTmp}/DNA087244.${inputType}" ]
	then
		echo "already copied, skipped"
	else
		scp -r calculon.hpc.rug.nl:${validationFolderPrm}/ "${validationFolderTmp}/"
		#scp calculon.hpc.rug.nl:${validationFolderPrm}/referenceCall/*{.gz,tbi} "${validationFolderTmp}/referenceCall/"
		#scp -r calculon.hpc.rug.nl:${validationFolderPrm}/filtered/* "${validationFolderTmp}/filtered/"
	fi
else
	validationFolderTmp=${inputFolder}/validationVcfs/
	if [ -f "${validationFolderTmp}/DNA087244.${inputType}" ]
        then
		echo "already copied, skipped"
        else
		cp -r calculon.hpc.rug.nl:${validationFolderPrm}/ "${validationFolderTmp}/"
	fi
fi

#doVariantEval
#doComparisonFiltered "no"
#doComparisonFiltered "referenceCall"
checkFrankenstein
