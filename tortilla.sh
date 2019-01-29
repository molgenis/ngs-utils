#!/bin/bash
set -eu


function showHelp(){

cat <<EOH
===============================================================================================================

Script requires one initial argument:

	-m|--makeSamplesheet        Creating an samplesheet for inhouse sequencing runs, no FastQ files are made jet. (makeSamplesheet_inhouse_research.sh)
	-s|--makeSamplesheetExternal    Creating an (external) samplesheet based on a inputfolder containing e.g. FastQ files (makeSamplesheet_external_samples.sh)
	-b|--bamout		Recreating the bam file for a certain region where the variant calling is based on (bamout.sh)
	-c|--countCoverage	Counting coverage (avg,med,sd,percentage 10/20/30/50/100x coverage) per Gene and target based on the panel that is given (countCoverage.sh)
	-v|--vcfCompare		Comparing 2 vcf files with eachother, this will output the differences + a vcf stats file (vcf-compare_2.0.sh)
	-n|--validateNGS	Script to check the known SNPs back in the NGS_DNA_Verification_test (checkValidationNGS_DNA.sh)
	-r|--revertBamToFastQ	go back from bam to fastq (paired end only)
	-cc|--calculateCoverage	CoveragePerBase or per Target calculations for a specific targetpanel
	-d|--cramToBam		converting cram files to bam(CramConversion.sh)
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
	${EBROOTNGSMINUTILS}/makeSamplesheet_inhouse_research.sh ${@}

elif [[ "${1}" == "--makeSamplesheetExternal" || "${1}" == "-s" ]]
then
	shift
	${EBROOTNGSMINUTILS}/makeSamplesheet_external_samples.sh ${@}

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
elif [[ "${1}" == "--calculateCoverage" || "${1}" == "-cc" ]]
then
	${EBROOTNGSMINUTILS}/coverage_calc.sh ${@}
elif [[ "${1}" == "--cramToBam" || "${1}" == "-d" ]]
then
	${EBROOTNGSMINUTILS}/CramConversion.sh ${@}
fi
