#!/usr/bin/bash

set -e
set -u

function showHelp() {
	#
	# Display commandline help on STDOUT.
	#
	cat <<EOH
===============================================================================================================
Script to make Bed files for Diagnostics.
Usage:
	$(basename $0) OPTIONS
Options:
	-h   Show this help.
	-b   Name of the BED file
	-n   Name of the new BED file
	-w   Bedfile directory --> where to write to (default: this directory)
	-e   Making BED file for exomekit true/false [default=false]
===============================================================================================================
EOH
	trap - EXIT
	exit 0
}


while getopts "b:n:e:w:h" opt;
do
	case $opt in h)showHelp;; b)bedfile="${OPTARG}";; n)name="${OPTARG}";; w)workDir="${OPTARG}";; e)exome="${OPTARG}";;
	esac
done

if [[ -z "${bedfile:-}" ]]
then
	echo -e '\nERROR: Must specify a BED file!\n'

	showHelp
	exit 1
fi

if [[ -z "${name:-}" ]]
then
	echo -e '\nERROR: Must specify a Name for the new Bed file!\n'

	showHelp
	exit 1
fi

if [[ -z "${exome:-}" ]]
then
	exome="false"
fi
if [[ -z "${workdir:-}" ]]
then
	workdir=$(pwd)
fi


if [ -d "/apps/data/Agilent/${name}" ]
then
	echo "WARNING: /apps/data/Agilent/${name} already exists, removing old one"
	rm -rf "/apps/data/Agilent/${name}"
fi	
if [ -d "/apps/data/UMCG/Diagnostics/${name}" ]
then
	echo "WARNING: /apps/data/UMCG/Diagnostics/${name} already exists, removing old one"
	rm -rf "/apps/data/UMCG/Diagnostics/${name}"
fi

umcgDir=/apps/data/UMCG/Diagnostics/

mkdir -p "${workdir}/${name}/human_g1k_v37/"
echo "created ${name}/human_g1k_v37/"
cp "${bedfile}" "${workdir}/${name}"/
echo "copied ${bedfile} ${workdir}/${name}/"
## navigate to folder
cd ${workdir}/"${name}"


cp "${bedfile}" "human_g1k_v37/captured.bed"
echo "copied ${bedfile} to human_g1k_v37/captured.bed"

module load ngs-utils
cd 'human_g1k_v37' || exit

## Run the prepare step

if [[ "${exome}" == 'true' ]]
then
	echo "Creating bedfiles for a new exomekit ${name}"
	"/apps/data/Agilent/prepare_NGS_Bedfiles.sh" -n 'captured'
elif [[ "${exome}" == 'false' ]]
then
	echo "Creating bedfiles for a new kit ${name}"
	"/apps/data/Agilent/prepare_NGS_Bedfiles.sh" -n 'captured' -c 'true' -d 'targeted'
else
	echo "please fill in true or false"
	exit 1
fi

##
cd "${workdir}"
echo "copied ${name} to ${umcgDir}"
cp -r "${name}" ${umcgDir}

cd "${umcgDir}/${name}/human_g1k_v37/"
echo "renaming captured into ${name}"
rename "captured" "${name}" "captured."*

#perbase
cd "${umcgDir}/CoveragePerBase/"
mkdir -p "${name}"


if [[ "${exome}" == 'false' ]]
then
	cd "${umcgDir}/CoveragePerBase/${name}" || exit
	ln -sf "../../${name}"/
fi

#pertarget
cd "${umcgDir}/CoveragePerTarget/"
mkdir -p "${name}"

if [[ "${exome}" == 'false' ]]
then
	cd "${umcgDir}/CoveragePerTarget/${name}" || exit
	ln -sf "../../${name}"/
fi

echo "FINISHED"

