#!/bin/bash
set -e
set -u


#declare -a searchBasePathsNGS=('/groups/umcg-gd/prm05' '/groups/umcg-gd/prm06')

function showHelp() {
	#
	# Display commandline help on STDOUT.
	#
	cat <<EOH
======================================================================================================================
Script to collect all available data from prm

Usage:
	$(basename "${0}") OPTIONS
Options:
	-h	Show this help.
	-c	recovering all the coverageFiles from prm
	-b	recovering all the bam and cram files from prm
	-v	recovering all the vcf files from prm
	-f	recovering all the fastq files from prm
	-a	recover all files from prm
	-o	location of the output file
	-p	prmlocations from where the script runs

===============================================================================================================
EOH

}

function _findCoverageReports() {
	local _searchBasePath
	## add writing header of the file + name
	printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' 'id' 'umcgID' 'familyID' 'dnaID' 'testID' 'filename' 'filepath' 'filetype' 'md5' 'dateCreated' \
		> "${outputLocation}/Coverage_reports_stored_prm.txt"
	
	for _searchBasePath in "${searchBasePathsNGS[@]}"
	do 
		readarray -t _filesFound< <(find "${_searchBasePath}/projects/" -wholename '*/results/coverage/*/*coverage*.txt')
		local _fileFound
		for _fileFound in "${_filesFound[@]}"
		do
			local _filePath=$(dirname "${_fileFound}")
			local _fileName=$(basename "${_fileFound}")
			_regex='^([0-9]+)_([0-9]+)_(D?N?A?[0-9]+)_.*'
			local _familyID='N/A'
			local _umcgID='N/A'
			local _dnaID='N/A'
			if [[ "${_fileName}" =~ ${_regex} ]]
			then
				_familyID="${BASH_REMATCH[1]}"
				_umcgID="${BASH_REMATCH[2]}"
				_dnaID="${BASH_REMATCH[3]}"
			fi
			#local _fileType=$(echo "${_fileName}" | awk 'BEGIN { FS = "." } ; { print $NF }')
			local _fileType="${_fileName##*.}"
			local _dateCreated=$(stat ${_fileFound} | grep 'Change' | awk 'BEGIN { FS = " " } ; { print $2 }')
			local _md5='n/a'
			if [[ -e "${_fileFound/\/results*/}.md5" ]]
			then
				_md5=$(grep "${_fileName}$" "${_fileFound/\/results*/}.md5" | cut -d ' ' -f 1)
			fi
			## Adding the information to a file
			printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' '' "${_umcgID}" "${_familyID}" "${_dnaID}" '' "${_fileName}" "${_filePath}" "${_fileType}" "${_md5}" "${_dateCreated}" \
				>> "${outputLocation}/Coverage_reports_stored_prm.txt"
		done
	done
}

function _findBamsAndCramfiles() {
	local _searchBasePath
	## add writing header of the file + name
	printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' 'id' 'umcgID' 'familyID' 'dnaID' 'testID' 'filename' 'filepath' 'filetype' 'md5' 'dateCreated' \
		> "${outputLocation}/Bams_Crams_files_stored_prm.txt"
	
	for _searchBasePath in "${searchBasePathsNGS[@]}"
	do 
		readarray -t _filesFound< <(find "${_searchBasePath}/projects/" -wholename '*/results/alignment/*.bam' -or -wholename '*/results/alignment/*.cram')
		local _fileFound
		for _fileFound in "${_filesFound[@]}"
		do
			local _filePath=$(dirname "${_fileFound}")
			local _fileName=$(basename "${_fileFound}")
			_regex='^([0-9]+)_([0-9]+)_(D?N?A?[0-9]+)_.*'
			local _familyID='N/A'
			local _umcgID='N/A'
			local _dnaID='N/A'
			if [[ "${_fileName}" =~ ${_regex} ]]
			then
				_familyID="${BASH_REMATCH[1]}"
				_umcgID="${BASH_REMATCH[2]}"
				_dnaID="${BASH_REMATCH[3]}"
			fi
			#local _fileType=$(echo "${_fileName}" | awk 'BEGIN { FS = "." } ; { print $NF }')
			local _fileType="${_fileName##*.}"
			local _dateCreated=$(stat ${_fileFound} | grep 'Change' | awk 'BEGIN { FS = " " } ; { print $2 }')
			local _md5='n/a'
			if [[ -e "${_fileFound/\/results*/}.md5" ]]
			then
				_md5=$(grep "${_fileName}$" "${_fileFound/\/results*/}.md5" | cut -d ' ' -f 1)
			fi
			## Adding the information to a file
			printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' '' "${_umcgID}" "${_familyID}" "${_dnaID}" '' "${_fileName}" "${_filePath}" "${_fileType}" "${_md5}" "${_dateCreated}" \
				>> "${outputLocation}/Bams_Crams_files_stored_prm.txt"
		done
	done
}

function _findVCFFiles() {
	local _searchBasePath
	## add writing header of the file + name
	printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' 'id' 'umcgID' 'familyID' 'dnaID' 'testID' 'filename' 'filepath' 'filetype' 'md5' 'dateCreated' \
		> "${outputLocation}/VCF_files_stored_prm.txt"
	
	for _searchBasePath in "${searchBasePathsNGS[@]}"
	do 
		readarray -t _filesFound< <(find "${_searchBasePath}/projects/" -wholename '*/results/variants/*.vcf.gz' -or -wholename '*/results/variants/*.vcf')
		local _fileFound
		for _fileFound in "${_filesFound[@]}"
		do
			local _filePath=$(dirname "${_fileFound}")
			local _fileName=$(basename "${_fileFound}")
			_regex='^([0-9]+)_([0-9]+)_(D?N?A?[0-9]+)_.*'
			local _familyID='N/A'
			local _umcgID='N/A'
			local _dnaID='N/A'
			if [[ "${_fileName}" =~ ${_regex} ]]
			then
				_familyID="${BASH_REMATCH[1]}"
				_umcgID="${BASH_REMATCH[2]}"
				_dnaID="${BASH_REMATCH[3]}"
			fi
			#local _fileType=$(echo "${_fileName}" | awk 'BEGIN { FS = "." } ; { print $NF }')
			local _fileType="${_fileName##*.}"
			local _dateCreated=$(stat ${_fileFound} | grep 'Change' | awk 'BEGIN { FS = " " } ; { print $2 }')
			local _md5='n/a'
			if [[ -e "${_fileFound/\/results*/}.md5" ]]
			then
				_md5=$(grep "${_fileName}$" "${_fileFound/\/results*/}.md5" | cut -d ' ' -f 1)
			fi
			## Adding the information to a file
			printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' '' "${_umcgID}" "${_familyID}" "${_dnaID}" '' "${_fileName}" "${_filePath}" "${_fileType}" "${_md5}" "${_dateCreated}" \
				>> "${outputLocation}/VCF_files_stored_prm.txt"
		done
	done
}


function _findFastQFiles() {
	local _searchBasePath
	## add writing header of the file + name
	printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' 'id' 'umcgID' 'familyID' 'dnaID' 'testID' 'Project' 'fastQname1' 'fastQname2' 'md5_1' 'md5_2' 'filepath' 'filetype' 'dateCreated1' 'dateCreated2' \
	> "${outputLocation}/FastQ_files_stored_prm.txt"

	for _searchBasePath in "${searchBasePathsNGS[@]}"
	do
		readarray -t _filesFound< <(find "${_searchBasePath}/rawdata/ngs" -regextype posix-basic -regex "${_searchBasePath}"'/rawdata/ngs/[0-9]\{6\}_[A-Z0-9_-]*/[0-9]\{6\}_[A-Z0-9_-]*\.csv')
		local _fileFound
		for _fileFound in "${_filesFound[@]}"
		do
			declare -a _sampleSheetColumnNames=()
			declare -A _sampleSheetColumnOffsets=()
			local _externalSampleIDFieldIndex
			declare -a _externalSampleID=()
			IFS="," read -r -a _sampleSheetColumnNames <<< "$(head -1 "${_fileFound}")"
			for (( _offset = 0 ; _offset < ${#_sampleSheetColumnNames[@]:-0} ; _offset++ ))
			do
				if [[ ! "${_sampleSheetColumnNames[${_offset}]}" == "" ]]
				then
					_sampleSheetColumnOffsets["${_sampleSheetColumnNames[${_offset}]}"]="${_offset}"
				fi
			done

			_externalSampleIDFieldIndex=$((${_sampleSheetColumnOffsets['externalSampleID']} + 1))
			_laneFieldIndex=$((${_sampleSheetColumnOffsets['lane']} + 1))
			_barcodeFieldIndex=$((${_sampleSheetColumnOffsets['barcode']} + 1))
			_projectFieldIndex=$((${_sampleSheetColumnOffsets['project']} + 1))

			readarray -t _externalSampleID < <(tail -n +2 "${_fileFound}" | cut -d "," -f "${_externalSampleIDFieldIndex}","${_laneFieldIndex}","${_barcodeFieldIndex}","${_projectFieldIndex}" --output-delimiter='|')
			
			if [[ "${#_externalSampleID[@]:-0}" -lt '1' ]]
			then
				echo "No externalSampleIDs found for ${_fileFound}"
				continue
			fi
			local _externalSample
			for _externalSample in "${_externalSampleID[@]}"
			do
				echo "${_externalSample}"
				local _sampleName
				_sampleName=$(echo "${_externalSample}" | awk '{split($0,a,"|"); print a[1]}')
				echo "${_sampleName}"
				_regex='^([0-9]+)_([0-9]+)_(D?N?A?[0-9]+)_.*'
				local _familyID='N/A'
				local _umcgID='N/A'
				local _dnaID='N/A'
				if [[ "${_sampleName}" =~ ${_regex} ]]
				then
					_familyID="${BASH_REMATCH[1]}"
					_umcgID="${BASH_REMATCH[2]}"
					_dnaID="${BASH_REMATCH[3]}"
				fi
				local _project=$(echo "${_externalSample}" | awk '{split($0,a,"|"); print a[2]}')
				local _lane=$(echo "${_externalSample}" | awk '{split($0,a,"|"); print a[3]}')
				local _barcode=$(echo "${_externalSample}" | awk '{split($0,a,"|"); print a[4]}')
				local _filePath=$(dirname "${_fileFound}")
				local _ngsRawdataName=$(basename "${_fileFound}" | cut -f1 -d'.')
				
				#Getting the names and creation dates of the fastQfiles
				if [[ -e "${_filePath}/${_ngsRawdataName}_L${_lane}_${_barcode}_1.fq.gz" ]]
				then
					local _fastQfile1="${_ngsRawdataName}_L${_lane}_${_barcode}_1.fq.gz"
					local _dateCreated_1=$(stat ${_filePath}/${_fastQfile1} | grep 'Change' | awk 'BEGIN { FS = " " } ; { print $2 }')
				else
					local _fastQfile1='N/A'
					local _dateCreated_1='N/A'
				fi
				if [[ -e "${_filePath}/${_ngsRawdataName}_L${_lane}_${_barcode}_2.fq.gz" ]]
				then
					local _fastQfile2="${_ngsRawdataName}_L${_lane}_${_barcode}_2.fq.gz"
					local _dateCreated_2=$(stat ${_filePath}/${_fastQfile2} | grep 'Change' | awk 'BEGIN { FS = " " } ; { print $2 }')
				else
					local _fastQfile2='N/A'
					local _dateCreated_2='N/A'
				fi
				#Getting the md5s from both fastq files if they are there
				if [[ -e "${_filePath}/${_ngsRawdataName}_L${_lane}_${_barcode}_1.fq.gz.md5" ]]
				then
					#local _md5_1=$(cat "${_filePath}/${_ngsRawdataName}_L${_lane}_${_barcode}"_1.fq.gz.md5 | cut -d ' ' -f 1)
					local _md5_1=$(cut -d ' ' -f 1 "${_filePath}/${_ngsRawdataName}_L${_lane}_${_barcode}"_1.fq.gz.md5)
				else
					local _md5_1='N/A'
				fi
				
				if [[ -e "${_filePath}/${_ngsRawdataName}_L${_lane}_${_barcode}_2.fq.gz.md5" ]]
				then
					#local _md5_2=$(cat "${_filePath}/${_ngsRawdataName}_L${_lane}_${_barcode}_2.fq.gz.md5" | cut -d ' ' -f 1)
					local _md5_2=$(cut -d ' ' -f 1 "${_filePath}/${_ngsRawdataName}_L${_lane}_${_barcode}"_2.fq.gz.md5)
				else
					local _md5_2='N/A'
				fi
				local _fileType="${_fastQfile1##*.}"

				## Adding the information to a file
				printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' '' "${_umcgID}" "${_familyID}" "${_dnaID}" '' "${_project}" "${_fastQfile1}" "${_fastQfile2}" "${_md5_1}" "${_md5_2}" "${_filePath}" "${_fileType}"  "${_dateCreated_1}" "${_dateCreated_2}" \
					>> "${outputLocation}/FastQ_files_stored_prm.txt"
			done
		done
	done
}

#
# Get commandline arguments.
#

coverageFiles='false'
bamCramFiles='false'
vcfFiles='false'
fastqFiles='false'
recoverAllFiles='false'


while getopts "o:p:fcbvah" opt
do
	case "${opt}" in
		h)
			showHelp
			;;
		c)
			coverageFiles='true'
			;;
		b)
			bamCramFiles='true'
			;;
		v)
			vcfFiles='true'
			;;
		f)
			fastqFiles='true'
			;;
		a)
			recoverAllFiles='true'
			;;
		o)
			outputLocation="${OPTARG}"
			;;
		p)
			searchBasePathsNGS="${OPTARG}"
			;;
		\?)
				echo "FATAL: Invalid option -${OPTARG}. Try $(basename "${0}") -h for help."
				;;
		:)
				echo "FATAL: Option -${OPTARG} requires an argument. Try $(basename "${0}") -h for help."
				;;
		*)
				echo "FATAL: Unhandled option. Try $(basename "${0}") -h for help."
				;;
	esac
done

#
# Check commandline options.
#
if [[ -z "${outputLocation:-}" ]]
then
	showHelp
	echo -e 'ERROR: Must specify an outputLocation with -o'
	exit 1
fi

if [[ -z "${searchBasePathsNGS:-}" ]]
then
	showHelp
	echo -e 'ERROR: Must specify the prm locations from where this script must run example :  /groups/umcg-gd/prm06'
	exit 1
fi

if [[ "${coverageFiles}" == 'true' ]]
then
	_findCoverageReports
fi

if [[ "${bamCramFiles}" == 'true' ]]
then
	_findBamsAndCramfiles
fi

if [[ "${vcfFiles}" == 'true' ]]
then
	_findVCFFiles
fi

if [[ "${fastqFiles}" == 'true' ]]
then
	_findFastQFiles
fi

if [[ "${recoverAllFiles}" == 'true' ]]
then 
	_findCoverageReports
	_findBamsAndCramfiles
	_findVCFFiles
	_findFastQFiles
fi
if [[ "${coverageFiles}" == 'false' && "${bamCramFiles}" == 'false' && "${vcfFiles}" == 'false' && "${recoverAllFiles}" == 'false' && "${fastqFiles}" == 'false' ]]
then
	showHelp
	echo -e 'ERROR: Must specify one of the options -c ; -b ; -v ; -f ; -a'
	exit 1
fi

