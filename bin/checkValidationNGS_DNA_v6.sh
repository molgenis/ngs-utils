#!/bin/bash

set -e
set -u

function showHelp() {
	#
	# Display commandline help on STDOUT.
	#
	cat <<EOH
===============================================================================================================
Script to do the "old" method of validating the NGS_DNA, this script should be combined with compareWithVKGL.sh to have a complete validation
Usage:
	$(basename "${0}") OPTIONS
Options:
	-h   Show this help.
	-i   inputFolder
	-t   inputType (vcf or vcf.gz) (default= vcf.gz)
	-o   outputFolder (default:\${workDir}/validationFolder/)
	-v   validationFolder, folder where the vcfs are with the SNPs that should be found back (default=/groups/umcg-gd/prm06/projects/validationVcfs/)
	-l   validationLevel (all|1|2) default is all
					1	is old validation (finding back some SNPs in 11 samples) + checkChromosomes
					2	frankenstein
					3	checkChromosomes
					all	is running both option 1 and 2 
===============================================================================================================
EOH
	trap - EXIT
	exit 0
}

function checkAllChromosomes(){
	## we need only first file, all files are coming from same project, thus same variants (not genotypes of course)
	firstInputFile=$(ls "${inputFolder}/"*".${inputType}" | head -1)
	chromosomesInFile=($(zcat "${firstInputFile}" | grep -v '^#' | awk '{print $1}' | sort -V | uniq))
	exitAfterLoopFinished="no"
	for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 MT NC_001422.1 X Y
	do
		if [[ ! " ${chromosomesInFile[*]} " == *" ${i} "* ]]
		then
			echo -e "\nCHROMOSOME: ${i} is not in the data!\n"
			exitAfterLoopFinished="yes"
		fi
	done
	# exit when there is a chromosome missing
	if [[ "${exitAfterLoopFinished}" == "yes" ]]
	then
		echo "there is/are chromosomes missing from this list: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 MT NC_001422.1 X Y"
		exit 1
	else
		echo -e "\ncheckAllChromosomes: All chromosomes are found back!\n"
	fi

}

function doVariantEval(){

	folder="${validationFolderTmp}"

	mapfile -t validationFiles < <(find "${folder}" -maxdepth 1 -name "*.${inputType}")
	if [[ "${#validationFiles[@]:-0}" -eq '0' ]]
	then
		echo "There are no files found in: ${folder}"
		exit 1
	else
		echo "##OPEN## COMPARING vcf files, all variants should be found back!##" >>  "${outputFolder}/output.txt"

		for i in "${validationFiles[@]}"
		do
			name=$(basename "${i}" ".${inputType}")
			inputFile=$(ls "${inputFolder}/"*"${name}"*".${inputType}")
			bgzippedInput=${inputFile%.*}.bgz

			zcat "${inputFile}" | bgzip -c > "${bgzippedInput}"
			tabix -p vcf "${bgzippedInput}"
			inputFile="${bgzippedInput}"

			# shellcheck disable=SC2154
			java -jar "${EBROOTGATK}/GenomeAnalysisTK.jar" \
			-T VariantEval \
			-R '/apps/data/1000G/phase1/human_g1k_v37_phiX.fasta' \
			-o "${outputFolder}/output.${name}.eval.grp" \
			--eval "${inputFile}" \
			--comp "${i}"


			check=$(set -e ; awk '{if (NR==5){if ($11 == "100.00"){print "correct"}}}' "${outputFolder}/output.${name}.eval.grp")
			if [ "${check}" == "correct" ]
			then
				zcat "${i}" | awk -v sample="${name}" 'BEGIN {OFS="  "}{if ($1 !~ /^#/){print sample,$1,$2,$4,$5,"FOUND BACK"}}' >> "${outputFolder}/output.txt"
			else
				zcat "${i}" | awk -v sample="${name}" 'BEGIN {OFS="  "}{if ($1 !~ /^#/){print sample,$1,$2,$4,$5,"Not 100% concordant!"}}' >> "${outputFolder}/output.txt"
			fi
		done
		echo "##CLOSE## COMPARING vcf files, all variants should be found back!##" >>  "${outputFolder}/output.txt"
	fi

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

	mapfile -t inputFiles < <(find "${folder}" -maxdepth 1 -name "*.${inputType}")
	if [[ "${#inputFiles[@]:-0}" -eq '0' ]]
	then
		echo "There are no files found in: ${folder}"
	else
		for i in "${inputFiles[@]}"
		do
			name=$(basename "${i}" ".${inputType}")

			if [[ "${inputType}" == "vcf.gz" ]]
			then
				validationSample=$(zcat "${i}" | grep -v '^#' | awk '{print $1"-"$2"-"$4"-"$5"-"$7}')
				echo "SAMPLE: ${inputFolder}/*${name}*.${inputType}"
				if zcat "${inputFolder}/"*"${name}"*".${inputType}" | grep 18598089 
				then
					echo "found it, proceed"
				else
					echo "position not found, exiting"
					exit 1
				fi
				inputSample=$(set -e ; zcat "${inputFolder}/"*"${name}"*".${inputType}" | grep 18598089 | awk '{print $1"-"$2"-"$4"-"$5"-"$7}')
				if [ "${validationSample}" == "${inputSample}" ]
				then
					if [ "${refCall}" == "referenceCall" ]
					then
						zcat "${i}" | awk -v sample="${name}" 'BEGIN {OFS="  "}{if ($1 !~ /^#/){print sample,$1,$2,$4,$5,$7,"FOUND BACK,REF CALL"}}' >> "${outputFolder}/output.txt"
					else
						zcat "${i}" | awk -v sample="${name}" 'BEGIN {OFS="  "}{if ($1 !~ /^#/){print sample,$1,$2,$4,$5,$7,"FOUND BACK"}}' >> "${outputFolder}/output.txt"
					fi
				else
					# shellcheck disable=SC2129
					echo -e "\n##OPEN## ${refCall}: Expected variant####" >>  "${outputFolder}/output.txt"
					zcat "${i}" | awk -v sample="${name}" 'BEGIN {OFS="  "}{if ($1 !~ /^#/){print sample,$1,$2,$4,$5,$7,"Not 100% concordant!"}}' >> "${outputFolder}/output.txt"
					zcat "${inputFolder}/"*"${name}"*".${inputType}" | grep '18598089' >> "${outputFolder}/output.txt"
					zcat "${i}" >> "${outputFolder}/output.txt"
					echo -e "\n##CLOSE## ${refCall}: Expected variant####" >>  "${outputFolder}/output.txt"

				fi
			fi

		done
	fi
}

function checkFrankenstein() {

	knownVariants="${validationFolderTmp}/ValidationSet.annotated.vcf.gz"
	knownMissingVariants="${validationFolderTmp}/knownMissingVariants.txt"

	totalNoOfKnownVariants=$(zcat "${knownVariants}" | grep -c -v '^#')

	echo "Testing Frankenstein on ${totalNoOfKnownVariants} known variants"

	tmpFolder="${outputFolder}/tmp"
	input="${inputFolder}/Frankenstein_1.final.vcf.gz"

	echo "comparing \"true\" set: /groups/umcg-atd/tmp01/ValidationSet.annotated.vcf.gz"
	echo "with new: ${input}"
	echo "working folder: ${tmpFolder}"

	bedtools intersect -header -b "${input}" -a "${knownVariants}" > "${tmpFolder}/SameAsFrank.vcf"
	bgzip -c "${tmpFolder}/SameAsFrank.vcf" > "${tmpFolder}/SameAsFrank.vcf.gz"

	totalNoOfKnownVariants=$(zcat "${knownVariants}" | grep -c -v '^#')
	totalSame=$(grep -c -v '^#' "${tmpFolder}/SameAsFrank.vcf")

	echo "found ${totalSame} of in total ${totalNoOfKnownVariants} back in the ${input}"

	tortilla.sh -v -1  "${tmpFolder}/SameAsFrank.vcf.gz" -2 "${knownVariants}" -o "${outputFolder}"


	if diff -q "${knownMissingVariants}" "${outputFolder}/notInVcf1.txt"
	then
		echo "Everything is OK, the missing variants are the same as the already missing variants"
	else
		echo "diff ${knownMissingVariants} ${outputFolder}/notInVcf1.txt"
		echo "OOPS, there are differences between the known missing variants and the sample!"
		diff -y "${knownMissingVariants}" "${outputFolder}/notInVcf1.txt"
	fi
}

while getopts "i:o:v:t:l:h" opt; 
do
	# shellcheck disable=SC2249
	# shellcheck disable=SC2220
	case "${opt}" in h)showHelp;; i)inputFolder="${OPTARG}";; o)outputFolder="${OPTARG}";; v)validationFolderPrm="${OPTARG}";; t)inputType="${OPTARG}";; l)validationLevel="${OPTARG}";;
esac 
done

if [[ -z "${inputFolder:-}" ]]; then showHelp ; echo "inputFolder is not specified" ; fi ; echo "inputFolder=${inputFolder}"
if [[ -z "${outputFolder:-}" ]]; then outputFolder="${inputFolder}/validationFolder/" ; mkdir -p "${outputFolder}" ; fi ; echo "outputFolder=${outputFolder}"
if [[ -z "${inputType:-}" ]]; then inputType="vcf.gz" ; fi ; echo "inputType=${inputType}"	
if [[ -z "${validationFolderPrm:-}" ]]; then validationFolderPrm="/groups/umcg-gd/prm06/projects/validationVcfs/" ; fi ; echo "validationFolderPrm=${validationFolderPrm}"
if [[ -z "${validationLevel:-}" ]]; then validationLevel="all" ; fi ; echo "validationLevel=${validationLevel}" 

ml GATK 
ml HTSlib
ml BEDTools
ml ngs-utils

rm -rf "${outputFolder}/tmp/"
mkdir -p "${outputFolder}/tmp/"

echo '' > "${outputFolder}/output.txt"

if [[ "${validationLevel}" != "all" && "${validationLevel}" != "1" && "${validationLevel}" != "2" && "${validationLevel}" != "3" ]] 
then
	echo "this is an unknown validationLevel [${validationLevel}]"
	echo "bye bye"
	exit 1
fi
if [[ "${validationLevel}" == "all" || "${validationLevel}" == "1" ]]
then

	whichHost=$(hostname -s)

	if [[ "${whichHost}" == "leucine-zipper"  || "${whichHost}" == "zinc-finger" ]]
	then
		validationFolderTmp="${outputFolder}/input/validationVcfs/"
		mkdir -p "${validationFolderTmp}/filtered/"

		echo "copying validationVcfs"
		if [[ -f "${validationFolderTmp}/DNA087244.${inputType}" ]]
		then
			echo "already copied, skipped"
		else
			rsync -av "chaperone:${validationFolderPrm}/" "${validationFolderTmp}/"
		fi
	else
		echo "please run on leucine-zipper or zinc-finger"
	fi

	doVariantEval
	doComparisonFiltered "findVariant"
	doComparisonFiltered "referenceCall"
	checkAllChromosomes
fi

if [[ "${validationLevel}" == "all" || "${validationLevel}" == "2" ]]
then
	validationFolderTmp="${outputFolder}/input/validationVcfs/Frankenstein/"
	checkFrankenstein
fi

if [[ "${validationLevel}" == "all" || "${validationLevel}" == "3" ]]
then
	validationFolderTmp="${outputFolder}/input/validationVcfs/Frankenstein/"
	checkAllChromosomes
fi
