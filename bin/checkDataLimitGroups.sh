#!/bin/bash

set -eu

if [[ -z "${1:-}" ]]
then
	echo "missing file as an argument \$1"
	echo "USAGE: bash checkDataLimitGroups.sh PATHTOFILE/myfile.txt"
	echo "file should have this format: directory DEFAULT    /ifs/rekencluster/umcgst10/groups/umcg-bios/tmp01             No    -     -     20.00T 11.4589T"
	exit 1
fi
## read  Quota vs. disk usage for file system umcgst10. mail in the UMCG HPC helpdesk --> write to mail.txt
## do not copy header in the file
## FORMAT looks like below:
## directory DEFAULT    /ifs/rekencluster/umcgst10/groups/umcg-bios/tmp01             No    -     -     20.00T 11.4589T 

# replacing all white characters with one tab
perl -pi -e 's| +|\t|g' "${1}"

module load cluster-utils

while read -r line
do
	groupname=$(echo "${line}" | awk '{print $3}' | awk 'BEGIN {FS="/"}{print $6}')
	limit=$(echo "${line}" | awk '{print $7}')
	used=$(echo "${line}" | awk '{print $8}')

	lastChar="${used: -1}"
	usedWithoutLastChar="${used::-1}"
	limitWithoutLastChar="${limit::-1}"
	limitInMb=$(echo "${limitWithoutLastChar}" | awk '{print $1 * 1000000}')

	if [[ "${lastChar}" == [0-9] || "${lastChar}" == 'k' ]]
	then
		#echo "only a couple of (kilo)bytes"
		usedInMb=0
	elif [[ "${lastChar}" == 'T' ]]
	then
		##used is in Tb
		usedInMb=$(echo "${usedWithoutLastChar}" | awk '{print $1 * 1000000}')
	elif [[ "${lastChar}" == 'G' ]]
	then
		##used is in Gb
		usedInMb=$(echo "${usedWithoutLastChar}" | awk '{print $1 * 1000}')
	else
		##echo "is already in Mb"
		usedInMb="${limitWithoutLastChar}"
	fi
	mapfile -t owners < <(colleagues -g "${groupname}" | sed -n -e '/owner/,/=/ p' | head -n -1 | awk '{if (NR>2 && $1 != "MIA"){print $0}}' | awk '{if($NF ~ />/){print substr($NF, 2, length($NF)-2)}else {print substr($(NF-1), 2, length($(NF-1))-2)}}')
	mapfile -t managers < <(colleagues -g "${groupname}" | sed -n -e '/manager/,/=/ p' | head -n -1 | awk '{if (NR>2 && $1 != "MIA"){print $0}}'| awk '{if($NF ~ />/){print substr($NF, 2, length($NF)-2)}else {print substr($(NF-1), 2, length($(NF-1))-2)}}')

	limitInTb=$(echo "${limitInMb}" | awk '{print $1 / 1000000}')
	usedInTb=$(echo "${usedInMb}" | awk '{print $1 / 1000000}')
	echo "################"
	if [[ "${usedInMb}" -gt "${limitInMb}" ]]
	then
		## adding managers to the owners array
		if [[ "${#managers[@]:-0}" -ne '0' ]]
		then
			for manager in "${managers[@]}"
			do
				owners+=("${manager}")
			done
		fi
		myMailText="Dear group owners/datamanagers of the group ${groupname},\n\nThe data on tmp01 on Gearshift has exceeded the quota for group ${groupname}.\nThe diskspace used is ${usedInTb}TB, the limit is ${limitInTb}TB.\nPlease cleanup on the tmp01 storage of the Gearshift HPC cluster. This storage system is almost full, when it is full nobody can do anything anymore.\nWe really need your help.\nIf it is really not possible to clean up, please contact the HPC helpdesk to increase the quota limit.\n\nCheers,\nHPC Helpdesk\n"
		echo "SUBJECT: Diskspace limits exceeded for group ${groupname}"
		echo -e "BODY: ${myMailText}"
		printf -v joined '%s;' "${owners[@]}"
		echo "MAILTO: ${joined%;}"
	else
		echo -e "${groupname}\tused:${usedInTb}TB\tlimit:${limitInTb}TB"
	fi
done<"${1}"
