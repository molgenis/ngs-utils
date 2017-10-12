#!/bin/bash
set -eu


function showHelp(){

cat <<EOH
===============================================================================================================

Script requires one initial argument:

	-m|--makeSamplesheet	Creating an (external) samplesheet based on a inputfolder containing e.g. FastQ files (makeSamplesheet.sh)
	-b|--bamout		Recreating the bam file for a certain region where the variant calling is based on (bamout.sh)
	-c|--countCoverage	Counting coverage (avg,med,sd,percentage 10/20/30/50/100x coverage) per Gene and target based on the panel that is given  (countCoverage.sh)
	-v|--vcfCompare		Comparing 2 vcf files with eachother, this will output the differences + a vcf stats file (vcf-compare_2.0.sh)
	-n|--validateNGS	Script to check the known SNPs back in the NGS_DNA_Verification_test (checkValidationNGS_DNA.sh)
	-r|--revertBamToFastQ	go back from bam to fastq (paired end only)
===============================================================================================================
EOH
trap - EXIT
exit 0
}

ml ngs-utils

if [ -z ${1:-} ]
then
	showHelp
fi

if [[ "${1}" == "--makeSamplesheet" || "${1}" == "-m" ]]
then
	shift
	${EBROOTNGSMINUTILS}/makeSamplesheet.sh ${@}

elif [[ "${1}" == "--bamout" || "${1}" == "-b" ]]
then
	shift
	${EBROOTNGSMINUTILS}/bamout.sh ${@}

elif [[ "${1}" == "--countCoverage" || "${1}" == "-c" ]]
then
	shift
	${EBROOTNGSMINUTILS}/countCoverage.sh ${@}

elif [[ "${1}" == "--vcfCompare" || "${1}" == "-v" ]]
then
	shift
	${EBROOTNGSMINUTILS}/vcf-compare_2.0.sh ${@}
elif [[ "${1}" == "--validateNGS" || "${1}" == "-n" ]]
then
	shift
	${EBROOTNGSMINUTILS}/checkValidationNGS_DNA.sh ${@}
elif [[ "${1}" == "--revertBamToFastQ" || "${1}" == "-r" ]]
then
	${EBROOTNGSMINUTILS}/revertFromBamToFastQ.sh ${@}
fi
