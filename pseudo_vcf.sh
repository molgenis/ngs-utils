#!/bin/bash

set -e
set -u


function showHelp() {
	#
	# Display commandline help on STDOUT.
	#
	cat <<EOH
===============================================================================================================
Script to pseudo anonimize vcf files. 

There is an option to search the database for samples, or you can give a folder containing all the data to be pseudo anonimized 
Usage:
	$(basename $0) [-s FILE] [-d FILE] [-m]
	$(basename $0) [-i FILE] [-f FILE] [-m]

Options:
	-h   Show this help.

   required:
	-s   search database (in combination with -d) (file containing DNA numbers and mapping in second column, tab seperated) (not in combination with -m)
	-d   database file (containing all the variant vcf files) (default is /groups/umcg-gd/tmp06/pseudo/AllVcfs.txt)
	-i   input folder containing vcf files (not in combination with -s) 
	-f   mapping file (in combination with -i)
	-g   which group (default = umcg-gd)

    optional:
	-m   manta data (default disabled) (manta file should be in the same place as the database file, naming should be {DATABASEFILE}_manta.txt e.g. AllVcfs_manta.txt)
	-p   which prm (e.g. prm06) default is both prm05 and prm06
	-w   working directory


EXTRA INFO:
there is the script for creating a database file (ngs-utils repo; indexing.sh)
Then the file is copied to /groups/umcg-gd/tmp06/pseudo/
===============================================================================================================
EOH
		exit 0
}

while getopts "i:g:d:s:f:mp:h" opt;
do
		case $opt in h)showHelp;; i)input="${OPTARG}";; g)group="${OPTARG}";; s)search="${OPTARG}";; d)database="${OPTARG}";; f)mapping="${OPTARG}";; m)mantaBool='-m';;  p)prm="${OPTARG}";;w)workDir="${OPTARG}";;
		esac
done

if [[ -z "${input:-}" && -z "${search:-}" ]]
then
	showHelp
	echo "neither input as search option is selected, one of them should be chosen" 
elif [[ -n "${input:-}" && -n "${search:-}" ]]
then
	showHelp
	echo "input and search option is selected, one of them should be chosen" 
### inputfolder
elif [[ -z "${search:-}" ]]
then
	pathway=input
	echo "hier komt een pad naar de files ${input}"
	if [[ -z "${mapping:-}" ]]
	then 
		showHelp 
		echo "no mapping file"
		exit 1 
	fi
### search db
elif [[ -z "${input:-}" ]]
then
	echo "search=${search}"
	pathway=search
	if [[ -z "${database:-}" ]]
	then
		database="/groups/umcg-gd/tmp06/pseudo/AllVcfs.txt"
	else
		echo "database=${database}"
	fi

fi

if [[ -z "${workDir:-}" ]]
then
		workDir="$(pwd)/pseudoTmp/"
else
		workDir="${workDir}/pseudoTmp/"
fi

mkdir -p "${workDir}/"{variants,manta}"/"{output,tmp,input}"/"

if [[ "${pathway}" == "input" ]]
then
	module load BCFtools
	while read line 
	do
		# get dna number from first column
		dnaNumber=$(echo "${line}" | awk '{print $1}')
		# get pseudo id from second column
		pseudo=$(echo "${line}" | awk '{print $2}')

		if ls "${input}/"*"${dnaNumber}"*
		then
			vcfFilePath=$(ls "${input}/"*"${dnaNumber}"*)
		else
			echo "There is not a file with ${dnaNumber} found"
			continue
		fi
	
		vcfFile=$(basename "${vcfFilePath}")
		sampleName=${vcfFile%%.*}
		echo "${vcfFilePath}"
		echo "bcftools view -h ${vcfFilePath} > ${workDir}/variants/tmp/header.txt"
		bcftools view -h "${vcfFilePath}" > "${workDir}/variants/tmp/header.txt"
		## pseudo anonimize sample
		perl -pi -e "s|${sampleName}|${pseudo}|g" "${workDir}/variants/tmp/header.txt"
		##removing all paths that starts with /groups/ because it can contain stuff like projectnames we do not want to show 
		sed -e 's#/groups/[^\(\), ]*#dummy#g' "${workDir}/variants/tmp/header.txt" > "${workDir}/variants/tmp/updatedheader.txt"

		bcftools reheader -h "${workDir}/variants/tmp/updatedheader.txt" -o "${workDir}/variants/output/${pseudo}.vcf.gz" "${vcfFilePath}"

done<"${mapping}"


elif [[ "${pathway}" == "search" ]]
then
	module load BCFtools
	## red input file
	while read line
	do

		# get dna number from first column
		dnaNumber=$(echo "${line}" | awk '{print $1}')
		# get pseudo id from second column
		pseudo=$(echo "${line}" | awk '{print $2}')
		echo "${dnaNumber} ${database}"
		
		if grep -m 1 "${dnaNumber}" "${database}"
		then
			vcfFilePath=$(grep -m 1 "${dnaNumber}" "${database}" | awk '{print $1}')
			## first check if there is a Gavin file, if so then use it
			resultsDir=$(dirname "${vcfFilePath}")
			if ssh -n chaperone "test -e ${resultsDir}/GAVIN/*${dnaNumber}*vcf.gz"
			then
				vcfFilePath=$(ssh -n chaperone "ls ${resultsDir}/GAVIN/*${dnaNumber}*vcf.gz")
			fi
		else
			echo "${dnaNumber} cannot be found back in the database (${database})"
			continue 
		fi

		vcfFile=$(basename "${vcfFilePath}")
		sampleName=${vcfFile%%.*}

		## copy file from prm to tmp
		rsync -av "chaperone:${vcfFilePath}" "${workDir}/variants/input/"

		bcftools view -h "${workDir}/variants/input/${vcfFile}" > "${workDir}/variants/tmp/header.txt"
		## pseudo anonimize sample
		perl -pi -e "s|${sampleName}|${pseudo}|g" "${workDir}/variants/tmp/header.txt"
		##removing all paths that starts with /groups/ because it can contain stuff like projectnames we do not want to show 
		sed -e 's#/groups/[^\(\), ]*#dummy#g' "${workDir}/variants/tmp/header.txt" > "${workDir}/variants/tmp/updatedheader.txt"

		bcftools reheader -h "${workDir}/variants/tmp/updatedheader.txt" -o "${workDir}/variants/output/${pseudo}.vcf.gz" "${workDir}/variants/input/${vcfFile}"

		if [[ -n "${mantaBool:-}" ]]
		then
			if [ -f "${database%.txt}_manta.txt" ]
			then
				if grep -s -m 1 "${dnaNumber}" "${database%.txt}_manta.txt" 
				then
					mantaFilePath=$(grep -m 1 "${dnaNumber}" ${database%.txt}_manta.txt | awk '{print $1}')
					mantaFile=$(basename "${mantaFilePath}")
					if ssh -n chaperone "test -e ${mantaFilePath}"
					then
						rsync -av  "chaperone:${mantaFilePath}" "${workDir}/manta/input/"
						bcftools view -h "${workDir}/manta/input/${mantaFile}" > "${workDir}/manta/tmp/header.txt"
						perl -pi -e "s|${sampleName}|${pseudo}|g" "${workDir}/manta/tmp/header.txt"
						##removing all paths that starts with /groups/ because it can contain stuff like projectnames we do not want to show 
						sed -e 's#/groups/[^\(\), ]*#dummy#g' "${workDir}/manta/tmp/header.txt" > "${workDir}/manta/tmp/updatedheader.txt"
						bcftools reheader -h "${workDir}/manta/tmp/updatedheader.txt" -o "${workDir}/manta/output/${pseudo}_diploid.vcf.gz" "${workDir}/manta/input/${mantaFile}"
					else
						echo "original file was not found"
						echo -e "original file was not found" > "${workDir}/manta/output/${pseudo}_diploid.vcf.gz.NOTFOUND"
					fi
				else
					echo "manta file for ${dnaNumber} not found in ${database%.txt}_manta.txt"
				fi
			else
				echo "${database%.txt}_manta.txt not found"
			fi 
		fi
	done<"${search}"

	echo "reheadering done, results can be found here: ${workDir}/variants/output/"

fi
