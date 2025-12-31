#!/bin/bash
set -eu


#path1="${1:-/groups/umcg-atd/tmp07/concordance/ngs/OA99/}"
#path2="${2:-/groups/umcg-atd/tmp07/concordance/ngs/OA99/compareWith/}"
path1="${1:-0}"
path2="${2:-0}"


if [[ "${path1}" == "0" || "${path2}" == "0" ]]
then
	echo "FATAL: please use 2 arguments for the input folders. \$1 are the the vcf files that you want to compare. \$2 is the folder the vcf files you want to COMPARE WITH"
fi

rm -f 'samplesheet.txt'

count=1
for i in $(ls "${path1}/"*oarray.txt)
do
	f1filename=$(basename "${i}")
	f1dnaNumber=$(echo "${f1filename}" | awk 'BEGIN {FS="_"}{print $3}')
	f1projectName=$(echo "${f1filename}" | awk 'BEGIN {FS="_"}{print $5}')
	for j in $(ls "${path2}/"*oarray.txt)
	do
		f2filename=$(basename "${j}")
		f2dnaNumber=$(echo "${f2filename}" | awk 'BEGIN {FS="_"}{print $3}')
		f2projectName=$(echo "${f2filename}" | awk 'BEGIN {FS="_"}{print $5}')
		echo -e "${f1projectName}\t${f2projectName}\t${f1dnaNumber}\t${f2dnaNumber}\tspId${count}" >> 'samplesheet.txt' 
		count=$((count+1))

	done


done

