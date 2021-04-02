#!/bin/bash

#
##
### Environment and bash sanity.
##
#
set -u
set -e
umask 0027
SCRIPT_NAME=$(basename $0)
if [[ "${BASH_VERSINFO}" -lt 4 || "${BASH_VERSINFO[0]}" -lt 4 ]]
then
	echo "Sorry, you need at least bash 4.x to use ${SCRIPT_NAME}." >&2
	exit 1
fi

#
# Make sure dots are used as decimal separator.
#
LANG='en_US.UTF-8'
LC_NUMERIC="${LANG}"

#
##
### Functions.
##
#

function _Usage() {
	echo
	echo 'Script to rename all FastQs matching a specified pattern in the current directory'
	echo 'from GenomeScan FastQ naming conventions to UMCG Genetics dept. FastQ naming conventions.'
	echo
	echo "    ${SCRIPT_NAME} [options]"
	echo
	echo " Options:"
	echo "    -v                         Enable verbose logging (optional)."
	echo "    -s YYMMDD                  sequencingStartDate"
	echo "    -f 'FastQ_filename_regex'  Regex pattern to find the FastQ files that must be renamed."
	echo "                               Note: the pattern must be single quoted to prevent expansion by the shell."
	echo "                               E.g. -f 'HTVTYBBXX_103373-003*'"
	echo "    -n                         Allow Ns in the DNA barcodes (optional)."
	echo "                               By default FastQ file with Ns in the DNA barcodes will be skipped and not renamed."
	echo
}

#
# Custom signal trapping functions (one for each signal) required to format log lines depending on signal.
#
function _trapSig() {
	local _trap_function="${1}"
	local _line="${2}"
	local _function="${3}"
	local _status="${4}"
	shift 4
	for _sig; do
		trap "${_trap_function} ${_sig} ${_line} ${_function} ${_status}" "${_sig}"
	done
}

function _trapHandler() {
	local _signal="${1}"
	local _problematicLine="${2}"
	local _function="${3}"
	local _exitStatus="${4}"
	local _errorMessage="Unknown error in function ${_function}."
	_reportFatalError "${_problematicLine}" "${_exitStatus}" "${_errorMessage}"
}

function _reportFatalError() {
	local _problematicLine="${1}"
	local _exitStatus="${2:-$?}"
	local _errorMessage="Unknown error."
	_errorMessage="${3:-${_errorMessage}}"
	#
	# Notify on STDOUT.
	#
	echo "
$(hostname) - ${SCRIPT_NAME}:${_problematicLine}: FATAL: exit code = ${_exitStatus}
$(hostname) - ${SCRIPT_NAME}:${_problematicLine}:        error message = ${_errorMessage}
"
	#
	# Reset trap and exit.
	#
	trap - EXIT
	exit 1
}


#
# Trap all exit signals: HUP(1), INT(2), QUIT(3), TERM(15), ERR
#
_trapSig '_trapHandler' '${LINENO}' '${FUNCNAME:-main}' '$?' HUP INT QUIT TERM EXIT ERR

function _RenameFastQ() {
	local _fastqPath="${1}"
	local _sequencingStartDate="${2}"
	echo "INFO: Processing FastQ ${_fastqPath}..."
	#
	# Example for FastQs from external lab:
	#  1. GenomeScan FastQ file name format:
	#         Flowcell_Customer-BatchNumber-SampleNumber_Barcode1-Barcode2_Lane_Read.fastq.gz
	#     E.g.:
	#         HTVTYBBXX_103373-003-001_AGGCAGAA-AGGCTTAG_L001_R1.fastq.gz
	#
	# Example of sequence read header line as found in Illumina FastQ files produced with Casava / bcl2fastq >= 1.8:
	#         @sequencer:run:flowcell:lane:tile:xCoordinateOfClusterInTile:yCoordinateOfClusterInTile sequenceReadOfPair:filtered:barcode[+barcode2]
	#     E.g.:
	#         @K00296:315:HVLN7BBXX:7:1101:2940:1173 1:Y:0:GTAGAGGA+TCTACTAT
	#
	#
	local _fastqDir="$(dirname "${_fastqPath}")"
	local _fastqFile="$(basename "${_fastqPath}")"
	#
	#  Get essential meta-data from the sequence read IDs in the FastQ file.
	#  (Do NOT rely on parsing FastQ filenames!)
	#
	local _firstReadID=$(zcat "${_fastqPath}" | head -1)
	if [[ "${enableVerboseLogging}" -eq 1 ]]
	then
		echo "DEBUG:    Found _firstReadID ............ = ${_firstReadID}"
	fi
	local _regex='^@([A-Z0-9][A-Z0-9]*):([0-9][0-9]*):([A-Z0-9][A-Z0-9]*):([1-8]):[0-9]*:[0-9]*:[0-9]* ([1-2]):[YN]:[0-9][0-9]*:[ATCGN][ATCGN+]*'
	if [[ "${_firstReadID}" =~ ${_regex} ]]
	then
		local _sequencer="${BASH_REMATCH[1]}"
		local _run="${BASH_REMATCH[2]}"
		local _flowcell="${BASH_REMATCH[3]}"
		local _lane="${BASH_REMATCH[4]}"
		local _sequenceReadOfPair="${BASH_REMATCH[5]}"
		#
		# Note: we do require the ID of the first read to contain a DNA barcode at the end,
		# but we don't parse barcodes here. See below for why...
		#
		#
		# Sanity check for run number and add leading zero when run number < 4.
		#
		if [[ "${#_run}" -lt 1 ]]
		then
			_reportFatalError ${LINENO} '1' 'Run number detected in ID of first read is too short (< 1): '"${_run:-}."
		elif [[ "${#_run}" -eq 1 ]]
		then
			_run="000${_run}"
		elif [[ "${#_run}" -eq 2 ]]
		then
			_run="00${_run}"
		elif [[ "${#_run}" -eq 3 ]]
		then
		_run="0${_run}"
		fi
		#
		if [[ "${enableVerboseLogging}" -eq 1 ]]
		then
			echo "DEBUG:    Found _sequencer .............. = ${_sequencer}"
			echo "DEBUG:    Found _run .................... = ${_run}"
			echo "DEBUG:    Found _flowcell ............... = ${_flowcell}"
			echo "DEBUG:    Found _lane ................... = ${_lane}"
			echo "DEBUG:    Found _sequenceReadOfPair ..... = ${_sequenceReadOfPair}"
		fi
	else
		_reportFatalError ${LINENO} '1' "Failed to parse required meta-data values from ID of first read of ${_fastqPath}".
	fi
	#
	# Due to sequencing errors the barcode may deviate slightly, so looking only at the first one won't fly.
	# In addition the start and end of a FastQ file tends to be enriched for sequencing errors / low quality.
	# Therefore we:
	#   A. use a combination of head and tail commands to skip the first 100,000 reads (400,000 lines)
	#      and get data from in the middle of the FastQ file.
	#      (assuming the yield was high enough; otherwise you get the last 1000 reads).
	#   B. Next we parse the barcodes from the read ID lines, sort, count and determine the most abundant barcode.
	#
	local _mostAbundandBarcode=$(zcat "${_fastqPath}" | head -n 440000 | tail -n 4000 | awk 'NR % 4 == 1' | awk -F ':' '{print $10}' | sort | uniq -c | sort -k 1,1nr | head -1 | awk '{print $2}' | tr -d '\n')
	local _barcodeRegex='^([ATCG][ATCG+]*)$'
	if [[ "${_mostAbundandBarcode}" =~ ${_barcodeRegex} ]]
	then
		local _barcodes="${BASH_REMATCH[1]}"
		#
		# In case the sample was double barcoded change the separator from '+' to '-'.
		#
		_barcodes="${_barcodes/+/-}"
		if [[ "${enableVerboseLogging}" -eq 1 ]]
		then
			echo "DEBUG:    Found _barcode(s) ............. = ${_barcodes}"
		fi
	elif [[ "${_mostAbundandBarcode}" =~ N ]]
	then
		if [[ "${allowN}" -eq '1' ]]
		then
			echo "WARN: Most abundant barcode(s) in max 1000 reads from middle of FastQ contains Ns: ${_mostAbundandBarcode}."
			echo "WARN: Will continue processing FastQ ${_fastqFile} despite poor sequencing quality of barcode(s), because commandline option -n was specified."
		else
			qualityControl='failed'
			echo "ERROR: Most abundant barcode(s) in max 1000 reads from middle of FastQ contains Ns: ${_mostAbundandBarcode}."
			echo "ERROR: Skipping discarded FastQ ${_fastqFile} due to poor sequencing quality of barcode(s)."
			return
		fi
	else
		qualityControl='failed'
		echo "ERROR: Failed to determine the most abundant barcodes from max 1000 reads from middle of FastQ."
		echo "ERROR: Failed to parse barcode(s) from read IDs of FastQ file ${_fastqFile}."
		return
	fi
	#
	# Get checksum
	#
	local _fastqChecksum=$(cat "${_fastqDir}/"*.md5 | grep "${_fastqFile}" | awk '{print $1}')
	#
	# Create sequence run subdir if it did not already exist.
	#
	local _newFastqDir="${_fastqDir}/${_sequencingStartDate}_${_sequencer}_${_run}_${_flowcell}"
	mkdir -p "${_newFastqDir}"
	#
	# Check if new FastQ file path does not already exist to prevent overwriting a sample with another one.
	#
	local _newFastqFile="${_sequencingStartDate}_${_sequencer}_${_run}_${_flowcell}_L${_lane}_${_barcodes}_${_sequenceReadOfPair}.fq.gz"
	if [[ -e "${_newFastqDir}/${_newFastqFile}" ]]
	then
		_reportFatalError ${LINENO} '1' "${_newFastqDir}/${_newFastqFile} already exists; will NOT move ${_fastqPath} -> ${_newFastqDir}/${_newFastqFile}."
	fi
	#
	# Move and rename FastQ + MD5 checksum on the fly.
	#
	echo "INFO:    Moving ${_fastqPath} -> ${_newFastqDir}/${_newFastqFile}"
	printf '%s  %s\n' "${_fastqChecksum}" "${_newFastqFile}" > "${_newFastqDir}/${_newFastqFile}.md5"
	mv "${_fastqPath}" "${_newFastqDir}/${_newFastqFile}"
}

#
##
### Main
##
#

#
# Get commandline arguments.
#
enableVerboseLogging=0 # Disabled by default.
allowN=0 # Do not allow Ns in the DNA barcodes by default.
while getopts "s:f:hnv" opt
do
	case ${opt} in
		h)
			_Usage
			#
			# Reset trap and exit.
			#
			trap - EXIT
			exit 0
			;;
		s)
			sequencingStartDate="${OPTARG}"
			;;
		f)
			fastqFilePattern="${OPTARG}"
			;;
		n)
			allowN=1
			;;
		v)
			enableVerboseLogging=1
			;;
		\?)
			_reportFatalError "${LINENO}" '1' "Invalid option -${OPTARG}. Try \"$(basename $0) -h\" for help."
			;;
		:)
			_reportFatalError "${LINENO}" '1' "Option -${OPTARG} requires an argument. Try \"$(basename $0) -h\" for help."
			;;
	esac
done

#
# Make sure there are no extra arguments we did not expect nor need.
#
shift $(($OPTIND - 1))
if [[ ! -z ${1:-} ]]
then
	_reportFatalError "${LINENO}" '1' "Invalid argument \"$1\". Try \"$(basename $0) -h\" for help."
fi

#
# Check if required args are present.
#
if [[ -z "${sequencingStartDate:-}" || -z "${fastqFilePattern:-}" ]]
then
	_reportFatalError "${LINENO}" '1' "One ore more required arguments is missing. Try \"$(basename $0) -h\" for help."
fi

#
# Check if sequencingStartDate is in the expected format.
#
ssd_regex='^[1-9][0-9][0-1][0-9][0-3][0-9]$'
if [[ "${sequencingStartDate}" =~ ${ssd_regex} ]]
then
	echo "INFO: Using sequencingStartDate ${sequencingStartDate}"
else
	_reportFatalError "${LINENO}" '1' "sequencingStartDate in unsupported format. Must be YYMMDD, but got ${sequencingStartDate}."
fi

#
# Process FastQ files.
#
qualityControl='unknown'
for FastQ in $(ls -1 ${fastqFilePattern})
do
	_RenameFastQ "${FastQ}" "${sequencingStartDate}"
done

#
# Reset trap and exit with 0 or 1 depending on wether QC was Ok or failed.
#
if [[ "${qualityControl}" == 'failed' ]]
then
	_reportFatalError "${LINENO}" '1' "One or more FastQ files failed QC and was not renamed!."
else
	echo "INFO: Finished successfully!"
	trap - EXIT
	exit 0
fi
