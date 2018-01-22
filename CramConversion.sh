#!/bin/bash
set -e
set -u

function showHelp() {
        #
        # Display commandline help on STDOUT.
        #
	cat <<EOH
===============================================================================================================
Script to convert cram back to bam
Usage:
	$(basename $0) OPTIONS
Options:
        -h	Show this help.

        Required:
        -i	inputFolder (where all the crams are)

        Optional:
	-w	path to workDir; creating a directory called BAM (default is this directory)
	-r	reference file (default: /apps/data/1000G/phase1/human_g1k_v37_phiX.fasta)

Output will be written in workDir/BAM directory with the name of the cram without the cram but with .bam extension
===============================================================================================================
EOH
	trap - EXIT
        exit 0
}

while getopts "i:w:r:h" opt; 
do
	case $opt in h)showHelp;; w)workDir="${OPTARG}";; i)inputFolder="${OPTARG}";; r)reference="${OPTARG}";;
        esac
done

if [[ -z "${inputFolder:-}" ]]; then showHelp ; echo "inputFolder is not specified" ; fi ; echo "inputFolder=${inputFolder}"
if [[ -z "${workDir:-}" ]]; then workDir=$(pwd)/BAM; mkdir -p ${workDir} ; fi ; echo "workDir=${workDir}"
if [[ -z "${reference:-}" ]]; then reference="/apps/data/1000G/phase1/human_g1k_v37_phiX.fasta"; fi ; echo "reference=${reference}"

module load io_lib/1.14.6-foss-2015b


#To convert from CRAM -> BAM do:
for cramFile in $(ls "${inputFolder}/"*.cram)
do
	fileWithoutExtension="$(basename "${cramFile%.*}")"
	echo "starting to convert ${cramFile}"
	scramble \
	-I cram \
	-O bam \
	-r "${reference}" \
	-m \
	-t 8 \
	"${cramFile}" \
	"${workDir}/${fileWithoutExtension}"

	echo "${fileWithoutExtension} converted to ${workDir}/${fileWithoutExtension}"
done
