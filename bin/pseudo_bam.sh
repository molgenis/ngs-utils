#!/bin/bash

set -e
set -u



function reheader(){
	workDir="${1}"
	sampleName="${2}"
	fileName="${3}"
	pseudo="${4}"

	if [ -f "${workDir}/${sampleName}.finished" ]
		then
				echo "skipping ${workDir}/input/${fileName} conversion, already converted"
		else
				samtools view -H "${workDir}/input/${fileName}" > "${workDir}/tmp/${sampleName}.header.sam"
				perl -pi -e s"|$sampleName|$pseudo|"g "${workDir}/tmp/${sampleName}.header.sam"
				echo "starting to reheader ${workDir}/input/${fileName} to ${pseudo}.reheader.bam"
				samtools reheader -P "${workDir}/tmp/${sampleName}.header.sam" "${workDir}/input/${fileName}" >  "${workDir}/tmp/${pseudo}.reheader.bam"
				mv -v "${workDir}/tmp/${pseudo}.reheader.bam" "${workDir}/output/${pseudo}.bam"
				touch "${workDir}/${sampleName}.finished"
		fi

}

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
		-i   input folder containing bam files (not in combination with -s) 
		-f   mapping file (in combination with -i) [
		-g   which group (default = umcg-gd)

	optional:
		-p   which prm (e.g. prm06) default is both prm05 and prm06
	-w   working directory


EXTRA INFO:
there is the script for creating a database file (ngs-utils repo; indexing.sh)
Then the file is copied to /groups/umcg-gd/tmp06/pseudo/
===============================================================================================================
EOH
		exit 0
}

while getopts "i:g:s:d:f:p:w:h" opt;
do
		case $opt in h)showHelp;; i)input="${OPTARG}";; g)group="${OPTARG}";; s)search="${OPTARG}";; d)database="${OPTARG}";; f)mapping="${OPTARG}";;  p)prm="${OPTARG}";; w)workDir="${OPTARG}";;
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
elif [[ -z "${input:-}" ]]
then
	echo "search=${search}"
	pathway="search"

	if [[ -z "${database:-}" ]]
	then
		database="/groups/umcg-gd/tmp06/pseudo/AllVcfs.txt"
	else
		echo "database=${database}"
	fi

elif [[ -z "${search:-}" ]]
then
	pathway="input"
	echo "hier komt een pad naar de files ${input}"
	if [[ -z "${mapping:-}" ]]
	then 
		showHelp 
		echo "no mapping file"
		exit 1 
	else
		echo "mapping=${mapping}"
	fi
fi

module load SAMtools

if [[ ${pathway} == 'search' ]]
then
	if [[ -z "${workDir:-}" ]]
	then
		echo "workdir"
		workDir="$(pwd)/pseudoTmp/bam/"
	else
		workDir="${workDir}/pseudoTmp/bam/"
	fi
		mkdir -p "${workDir}/"{input,tmp,output}
	module load BCFtools
	## read input file
	while read line
	do
		# get dna number from first column
		dnaNumber=$(echo "${line}" | awk '{print $1}')
		# get pseudo id from second column
		pseudo=$(echo "${line}" | awk '{print $2}')
		filePath=$(grep -m 1 "${dnaNumber}" "${database}" | awk '{print $1}')
		vcfFile=$(basename "${filePath}")
		sample=${vcfFile%%.*}
		notFound="false"
		
		if grep -m 1 "${dnaNumber}" "${database}"
		then
			
			if [[ "${filePath}" == *"GAVIN"* ]]
			then
				resultsDir=$(dirname $(dirname "$(dirname ${filePath})"))
			else
				resultsDir=$(dirname "$(dirname ${filePath})")
			fi
			if ssh -n chaperone "test -e ${resultsDir}/alignment/${sample}*bam"
			then
				## copy file from prm to tmp
				echo "chaperone:${resultsDir}/alignment/${sample}*bam" "${workDir}/input/"
				rsync -av "chaperone:${resultsDir}/alignment/${sample}*bam" "${workDir}/input/"
			elif ssh -n chaperone "test -e ${resultsDir}/alignment/${sample}*cram"
			then
				 echo "${sample}.cram is found, please convert it to bam and rerun this again" >> "${workDir}/tmp/notFound.txt"
			else
				notFound="true"
				echo "${resultsDir}/alignment/${sample}*bam not found on prm05/prm06"
				echo "${sample}" >> "${workDir}/tmp/notFound.txt"
				continue
			fi
		else	
			echo "${dnaNumber} cannot be found back in the database (${database})"
			echo "${dnaNumber}" >> "${workDir}/tmp/notFound.txt"
			continue 
		fi
		bam=$(ls ${workDir}/input/*${dnaNumber}*.bam) ## 20000_DNA12345_000_12312.merged.bam
		fileName=$(basename "${bam}") ## 20000000_DNA12345_0000000_1231244
		sampleName=${fileName%%.*} ## 20000000_DNA12345_0000000_1231244

		reheader "${workDir}" "${sampleName}" "${fileName}" "${pseudo}"
	done<"${search}"
	
	echo "pseudo anonimizing done, results can be found here: ${workDir}/output/"
	echo "samples that could not be found back are reported here: ${workDir}/tmp/notFound.txt"
fi

if [[ "${pathway}" == 'input' ]]
then
	if [[ -z "${workDir:-}" ]]
	then
		workDir="${input}/pseudoTmp/bam/"
	else
		workDir="${workDir}/pseudoTmp/bam/"
	fi

	mkdir -p "${workDir}"/{input,tmp,output}
	mv "${input}/"*".bams" "${workDir}/input/"
	module load SAMtools

	while read line
	do
		oldKey=$(echo "${line}" | awk '{print $1}')	 ## DNA12345
		pseudo=$(echo "${line}" | awk '{print $2}')	 ##sample1
		bam=$(ls "${workDir}/input/"*"${oldKey}"*".bam") ## 20000_DNA12345_000_12312.merged.bam
		fileName=$(basename "${bam}") ## 20000000_DNA12345_0000000_1231244
		sampleName=${fileName%%.*} ## 20000000_DNA12345_0000000_1231244
		echo "BAM: ${bam} ${oldKey} ${pseudo} ${sampleName}"

		reheader "${workDir}" "${sampleName}" "${fileName}" "${pseudo}"

	done<"${mapping}"
fi
