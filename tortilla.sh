#!/bin/bash
set -eu


function showHelp(){

cat <<EOH
===============================================================================================================

Script requires one initial argument:

	MakeSamplesheet		Creating an (external) samplesheet based on a inputfolder containing e.g. FastQ files
	Bamout			Recreating the bam file for a certain region where the variant calling is based on
	CountCoverage		Counting coverage (avg,med,sd,percentage 10/20/30/50/100x coverage) per Gene and target based on the panel that is given 
	VcfCompare		Comparing 2 vcf files with eachother, this will output the differences + a vcf stats file

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

if [[ "${1}" == "MakeSampleSheet" ]]
then
	shift
	/home/umcg-rkanninga/makeSamplesheet.sh ${@}
elif [[ "${1}" == "Bamout" ]]
then
	shift
	$EBROOTNGSMINUTILS/bamout.sh ${@}

elif [[ "${1}" == "CountCoverage" ]]
then
	shift
	$EBROOTNGSMINUTILS/CountCoverage.sh ${@}

elif [[ "${1}" == "VcfCompare" ]]
then
	shift
	$EBROOTNGSMINUTILS/vcf-compare_2.0.sh
fi
