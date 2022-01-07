#!/bin/bash
set -eu


function showHelp(){

cat <<EOH
===============================================================================================================

Script requires one initial argument:

    -m|--makeSamplesheetSimple          Create an simple sample sheet, output with DummyStartdate, DummySequencer, DummyRun (makeSamplesheet.sh)
    -w|--makeSamplesheetWithFastQ       Creating a samplesheet for (external) samples based on a inputfolder containing e.g. FastQ files 
                                        Information is received from the FastQ file. (makeSamplesheet_externalSamplesWithFastQ.sh)
    -n|--makeSamplesheetNoFastQ         Creating a samplesheet for inhouse sequencing runs, no FastQ files information is needed. (makeSamplesheet_InhouseNoFastQinfo.sh)
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

if [[ "${1}" == "--makeSamplesheetSimple" || "${1}" == "-m" ]]
then
	shift
	"${EBROOTNGSMINUTILS}/makeSamplesheet.sh" "${@}"

elif [[ "${1}" == "--makeSamplesheetWithFastQ" || "${1}" == "-w" ]]
then
	shift
	"${EBROOTNGSMINUTILS}/makeSamplesheet_externalSamplesWithFastQ.sh" "${@}"

elif [[ "${1}" == "--makeSamplesheetNoFastQ" || "${1}" == "-n" ]]
then
	shift
	"${EBROOTNGSMINUTILS}/makeSamplesheet_InhouseNoFastQinfo.sh" "${@}"
fi


