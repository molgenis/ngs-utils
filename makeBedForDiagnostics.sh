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
	-d   Name of the directory
	-e   Making BED file for exomekit [default=false]
===============================================================================================================
EOH
	trap - EXIT
	exit 0
}


while getopts "b:d:e:h" opt;
do
	case $opt in h)showHelp;; b)bedfile="${OPTARG}";; d)directory="${OPTARG}";; e)exome="${OPTARG}";;
	esac
done

if [[ -z "${bedfile:-}" ]]
then
	echo -e '\nERROR: Must specify a BED file!\n'

	showHelp
        exit 1
fi


if [[ -z "${directory:-}" ]]
then
	echo -e '\nERROR: Must specify a directory!\n'

        showHelp
        exit 1
fi

if [[ -z "${exome:-}" ]]
then
	exome="false"
else
	exome="true"
fi


if [ -d /apps/data/Agilent/$directory ]
then
	echo "/apps/data/Agilent/$directory already exists"
	exit 1
elif [ -d /apps/data/UMCG/Diagnostics/$directory ]
then
	echo "/apps/data/UMCG/Diagnostics/$directory already exists"
	exit 1
fi

thisDir=$(pwd)
UMCGDir=/apps/data/UMCG/Diagnostics/

mkdir -p ${directory}/human_g1k_v37/
echo "created ${directory}/human_g1k_v37/"
cp ${bedfile} ${directory}/
echo "copied ${bedfile} ${directory}/"
## navigate to folder
cd ${directory}


cp ${bedfile} human_g1k_v37/captured.bed
echo "copied ${bedfile} to human_g1k_v37/captured.bed"

module load ngs-utils

cd human_g1k_v37/

## Run the prepare step

if [[ ${exome} == 'true' ]]
	then
	echo 'Creating Bedfiles for a New ExomeKit'
	sh ${EBROOTNGSMINUTILS}/prepare_NGS_Bedfiles.sh -n captured
else
	echo "Creating Bedfiles for a New kit ${directory}"
	sh ${EBROOTNGSMINUTILS}/prepare_NGS_Bedfiles.sh -n captured -c true -d targeted
fi

##
cd $thisDir
echo "copied ${directory} to ${UMCGDir}"
cp -r ${directory} ${UMCGDir}

cd ${UMCGDir}/${directory}/human_g1k_v37/
echo "renaming captured into ${directory}"
rename captured ${directory} captured.*



#perbase
cd ${UMCGDir}/CoveragePerBase/
mkdir ${directory}
cd ${directory}
ln -sf ../../${directory}/

#pertarget
cd ${UMCGDir}/CoveragePerTarget/
mkdir ${directory}
cd ${directory}
ln -sf ../../${directory}/

echo "FINISHED"

