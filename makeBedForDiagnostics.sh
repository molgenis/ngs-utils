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
	-e   Making BED file for exomekit [default=false]
===============================================================================================================
EOH
	trap - EXIT
	exit 0
}


while getopts "b:d:e:h" opt;
do
	case $opt in h)showHelp;; b)bedfile="${OPTARG}";; d)name="${OPTARG}";; e)exome="${OPTARG}";;
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
else
	exome="true"
fi


if [ -d /apps/data/Agilent/"${name}" ]
then
	echo "/apps/data/Agilent/"${name}" already exists"
	exit 1
elif [ -d /apps/data/UMCG/Diagnostics/"${name}" ]
then
	echo "/apps/data/UMCG/Diagnostics/"${name}" already exists"
	exit 1
fi

thisDir=$(pwd)
umcgDir=/apps/data/UMCG/Diagnostics/

mkdir -p "${name}"/human_g1k_v37/
echo "created "${name}"/human_g1k_v37/"
cp "${bedfile}" "${name}"/
echo "copied "${bedfile}" "${name}"/"
## navigate to folder
cd "${name}"


cp "${bedfile}" human_g1k_v37/captured.bed
echo "copied "${bedfile}" to human_g1k_v37/captured.bed"

module load ngs-utils

cd human_g1k_v37/

## Run the prepare step

if [[ ${exome} == 'true' ]]
then
	echo 'Creating Bedfiles for a New ExomeKit'
	sh ${EBROOTNGSMINUTILS}/prepare_NGS_Bedfiles.sh -n captured
else
	echo "Creating Bedfiles for a New kit ${name}"
	sh ${EBROOTNGSMINUTILS}/prepare_NGS_Bedfiles.sh -n captured -c true -d targeted
fi

##
cd $thisDir
echo "copied "${name}" to ${umcgDir}"
cp -r "${name}" ${umcgDir}

cd ${umcgDir}/"${name}"/human_g1k_v37/
echo "renaming captured into "${name}""
rename captured "${name}" captured.*



#perbase
cd ${umcgDir}/CoveragePerBase/
mkdir "${name}"
cd "${name}"
ln -sf ../../"${name}"/

#pertarget
cd ${umcgDir}/CoveragePerTarget/
mkdir "${name}"
cd "${name}"
ln -sf ../../"${name}"/

echo "FINISHED"

