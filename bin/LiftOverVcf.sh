#!/bin/bash
set -eu


function showHelp() {
	#
	# Display commandline help on STDOUT.
	#
	cat <<EOH
===============================================================================================================
Script to Liftover a vcf file
Usage:
	$(basename "${0}") OPTIONS
Options:
	-h   Show this help.

    Required:
	-i   specify inputFile or inputFolder
	
    Optional:
	-o   outputFolder (default:creating a directory \${CURRENTDIR}/Liftover)
	outputfiles will be \${inputfilename%%.*}.LIFTOVER.vcf and \${inputfilename%%.*}.REJECT.vcf
	-c   which chain file to use (default is /apps/data/GRC/Hg19toHumanG1kV37.chain)   
	-r   What is target reference (default is /apps/data/1000G/phase1/human_g1k_v37_phiX.fasta)
	-s   recursive search for vcfs (only in combination when a directory is the input) (disabled by default) 
===============================================================================================================
EOH
	trap - EXIT
	exit 0
}

declare recursiveSearch='false'
while getopts "i:c:r:o:sh" opt;
do
	# shellcheck disable=SC2249
	# shellcheck disable=SC2220
	case "${opt}" in h)showHelp;; i)input="${OPTARG}";; o)outputFolder="${OPTARG}";;  r)reference="${OPTARG}";; c)chainFile="${OPTARG}";; s)recursiveSearch='true';;
esac
done

if [[ -z "${input:-}" ]]; then showHelp ; echo "input is not specified" ; fi ; echo "input=${input}"
if [[ -z "${outputFolder:-}" ]]; then mkdir -p "$(pwd)/Liftover" ; outputFolder="$(pwd)/Liftover/" ; fi ; echo "outputFolder=${outputFolder}"
if [[ -z "${reference:-}" ]]; then reference="/apps/data/1000G/phase1/human_g1k_v37_phiX.fasta" ; fi ; echo "reference=${reference}"
if [[ -z "${chainFile:-}" ]]; then chainFile="/apps/data/GRC/Hg19toHumanG1kV37.chain" ; fi ; echo "chainFile=${chainFile}"
if [[ "${recursiveSearch}" == 'true' ]]
then
	 
	echo "recursive search turned on"
	search=""
else 
	echo "no recursive folder search for vcf files"
	search="-maxdepth 1"
fi

module load picard

if [[ -d "${input}" ]]
then
	echo "${input} is a directory"
	# shellcheck disable=SC2086 
	mapfile -t inputFiles < <(find "${input}" ${search} -name "*.vcf" -o -name "*.vcf.gz")
elif [[ -f "${input}" ]]
then
    echo "${input} is a file"
	mapfile -t inputFiles < <(find "${input}")
else
    echo "${input} is not valid"
    exit 1
fi

if [[ "${#inputFiles[@]:-0}" -eq '0' ]]
then	
	echo "no Inputfiles found here: ${input}"
	exit 1
fi
for inputFile in "${inputFiles[@]}"
do
	fileName=$(basename "${inputFile}")
	outputFilePrefix="${fileName%%.*}"
	
	# shellcheck disable=SC2154 
	java -jar -Xmx15g "${EBROOTPICARD}/picard.jar" LiftoverVcf \
	I="${inputFile}" \
	O="${outputFolder}/${outputFilePrefix}.LIFTOVER.vcf" \
	REJECT="${outputFolder}/${outputFilePrefix}.REJECT.vcf" \
	R="${reference}" \
	CHAIN="${chainFile}"

	echo -e "Input=${inputFile}\nReference=${reference}\nChain=${chainFile}\noutputFile=${outputFolder}/${outputFilePrefix}.LIFTOVER.vcf" >> "${outputFolder}/runparameters.txt"
done

echo "output can be found: ${outputFolder}"
