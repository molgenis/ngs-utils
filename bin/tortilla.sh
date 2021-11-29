#!/bin/bash
set -eu


function showHelp(){

cat <<EOH
===============================================================================================================

Script requires one initial argument:

    -m|--makeSamplesheet        Go to tortilla_makeSamplesheet, for different make samplesheet options. (tortilla_makeSamplesheet.sh)
    -b|--bamout                 Recreating the bam file for a certain region where the variant calling is based on (bamout.sh)
    -c|--countCoverage          Counting coverage (avg,med,sd,percentage 10/20/30/50/100x coverage) per Gene and target based on the panel that is given (countCoverage.sh)
    -v|--vcfCompare             Comparing 2 vcf files with eachother, this will output the differences + a vcf stats file (vcf-compare_2.0.sh)
    -n|--validateNGS            Script to check the known SNPs back in the NGS_DNA_Verification (checkValidationNGS_DNA_v6.sh)
    -cv|--compareWithVKGL       script that compares vcf files with the VKGL standard
    -r|--revertBamToFastQ       go back from bam to fastq (paired end only)
    -cc|--calculateCoverage     CoveragePerBase or per Target calculations for a specific targetpanel
    -d|--cramToBam              Converting cram files to bam(CramConversion.sh)
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
	${EBROOTNGSMINUTILS}/bin/tortilla_makeSamplesheet.sh ${@}

elif [[ "${1}" == "--bamout" || "${1}" == "-b" ]]
then
	shift
	${EBROOTNGSMINUTILS}/bin/bamout.sh ${@}

elif [[ "${1}" == "--countCoverage" || "${1}" == "-c" ]]
then
	shift
	${EBROOTNGSMINUTILS}/bin/countCoverage.sh ${@}

elif [[ "${1}" == "--vcfCompare" || "${1}" == "-v" ]]
then
	shift
	${EBROOTNGSMINUTILS}/bin/vcf-compare_2.0.sh ${@}
elif [[ "${1}" == "--validateNGS" || "${1}" == "-n" ]]
then
	shift
	${EBROOTNGSMINUTILS}/bin/checkValidationNGS_DNA_v6.sh ${@}
elif [[ "${1}" == "--compareWithVKGL" || "${1}" == "-cv" ]]
then
	shift
	${EBROOTNGSMINUTILS}/bin/compareWithVKGL.sh ${@}
elif [[ "${1}" == "--revertBamToFastQ" || "${1}" == "-r" ]]
then
	${EBROOTNGSMINUTILS}/bin/revertFromBamToFastQ.sh ${@}
elif [[ "${1}" == "--calculateCoverage" || "${1}" == "-cc" ]]
then
	${EBROOTNGSMINUTILS}/bin/coverage_calc.sh ${@}
elif [[ "${1}" == "--cramToBam" || "${1}" == "-d" ]]
then
	${EBROOTNGSMINUTILS}/bin/CramConversion.sh ${@}
fi
