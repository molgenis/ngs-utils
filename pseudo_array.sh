#!/bin/bash

set -e
set -u


function showHelp() {
	#
	# Display commandline help on STDOUT.
	#
	cat <<EOH
===============================================================================================================
This script will pseudo anonimize array data from NxClinical 

Usage:
	$(basename $0) [-i FILE] [-f FILE] [-m]

Options:
	-h   Show this help.

   required:
	-i   input folder containing vcf files (not in combination with -s) 
	-f   mapping file (in combination with -i)
	-g   which group (default = umcg-gap)

    optional:
	-w   working directory

===============================================================================================================
EOH
		exit 0
}

while getopts "i:g:d:s:f:mp:h" opt;
do
		case $opt in h)showHelp;; i)input="${OPTARG}";; g)group="${OPTARG}";; f)mapping="${OPTARG}";; w)workDir="${OPTARG}";;
		esac
done

### search db
if [[ -z "${input:-}" ]]
then
	echo "an inputfolder is required"
	showHelp
	exit 1
fi

if [[ -z "${workDir:-}" ]]
then
		workDir="$(pwd)/pseudoTmp/"
else
		workDir="${workDir}/pseudoTmp/"
fi

mkdir -p "${workDir}/"{output,tmp}
echo "WORKDIR=${workDir}"

module load BCFtools
while read line 
do
	# get dna number from first column
	dnaNumber=$(echo "${line}" | awk '{print $1}')
	# get pseudo id from second column
	pseudo=$(echo "${line}" | awk '{print $2}')
	if ls "${input}/${dnaNumber}"*
	then
		arrayFilePath=$(ls "${input}/${dnaNumber}"*)
	else
		echo "There is not a file with ${dnaNumber} found"
		continue
	fi
	dos2unix "${arrayFilePath}"
	arrayFile=$(basename "${arrayFilePath}")
	## get extension
	fileExtension="${arrayFile##*.}"
	## punt extensions(s) eraf
	sampleName="${arrayFile%%.*}"
	
	projectName=$(grep 'Project Name' "${arrayFilePath}" | awk 'BEGIN {FS="= "}{print $2}')
	## pseudo anonimize sample
	perl -p -e "s|${sampleName}|${pseudo}|g" "${arrayFilePath}" > "${workDir}/tmp/${sampleName}.${fileExtension}"
	perl -pi -e "s|${projectName}|dummyProject|g" "${workDir}/tmp/${sampleName}.${fileExtension}"
	
	cp "${workDir}/tmp/${sampleName}.${fileExtension}" "${workDir}/output/${pseudo}.${fileExtension}"

done<"${mapping}"

echo "output can be found here: ${workDir}/output/"

