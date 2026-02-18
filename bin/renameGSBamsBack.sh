#!/bin/bash
set -eu

gsBatch="${1:-0}"
projectName="${2:-0}"

if [[ "${gsBatch}" == "0" || "${projectName}" == "0" ]]
then
	echo "FATAL: please use 2 arguments for gsBatch and projectName. \$1 is gsBatch (e.g. 103601-002) \$2 is the projectname (e.g. GS_315-WGS_v1"
	exit 1
fi

oldPath="/groups/umcg-genomescan/tmp06/${gsBatch}/Analysis/"
bamLocation="/groups/umcg-genomescan/tmp06/projects/NGS_DNA/${projectName}/run01/results/alignment/"
count=0

while read line
do
	if [[ "${count}" == 0 ]]
	then
		echo "first line"
		count=1
	else 
		GSID=$(echo "${line}" | awk 'BEGIN {FS=","}{print $1}')
		identifier=$(echo "${line}" | awk 'BEGIN {FS=","}{print $2}')
		sampleprocessStepid=$(echo ${identifier} | awk 'BEGIN {FS="-"}{print $3}')
		oldID="${GSID}-${identifier}"
	
		for i in bam bam.bai bam.md5 bam.md5sum
		do
			filename=$(find ${bamLocation} -name *${sampleprocessStepid}.${i} 2>/dev/null | wc -l)
			if [[ "${filename}" -eq 0 ]]
			then
				echo "there is no such file as: ${bamLocation}/*${sampleprocessStepid}.${i}"
			else
				mv -vf "${bamLocation}/"*"${sampleprocessStepid}.${i}" "${oldPath}/${oldID}/${oldID}.${i}"
			fi
		done
	fi
done<"UMCG_CSV_${gsBatch}.csv"