#!/bin/bash

set -eu

SCRIPT_NAME="$(basename ${0})"
SCRIPT_NAME="${SCRIPT_NAME%.*sh}"

function showHelp() {
	#
	# Display commandline help on STDOUT.
	#
	cat <<EOH
======================================================================================================================
Script to check whether the samplesheet is good
Usage:
	$(basename $0) OPTIONS
Options:
	-h   Show this help.
	-s   Scratch directory.
======================================================================================================================
EOH
	trap - EXIT
	exit 0
}


#
# Get commandline arguments.
#
declare group=''
while getopts "s:h" opt
do
	case $opt in
		h)
			showHelp
			;;
		s)
			SCR_DIR="${OPTARG}"
			;;
	esac
done

#
# Check commandline options.
#
if [[ -z "${SCR_DIR:-}" ]]
then
	echo 'Must specify a SCR_DIR (e.g. /groups/umcg-gd/scr01/).'
fi

sampleSheetsDir="${SCR_DIR}/Samplesheets"
echo "INFO: Processing samplesheets from ${sampleSheetsDir}/new/..."

declare -a sampleSheets=($(find "${sampleSheetsDir}/"*"/new/" -name '*.csv'))
if [[ "${#sampleSheets[@]:-0}" -eq '0' ]]
then
	echo "WARN: No samplesheets found in ${sampleSheetsDir}/*/new/."
	exit 0
else
	echo "DEBUG: samplesheets found: ${sampleSheets[@]}."
fi

for sampleSheet in "${sampleSheets[@]}"
do
	#
	# Create a copy of the original, so we preserve the original owner 
	# and can check who to notify when something is wrong with the samplesheet.
	#
	cp "${sampleSheet}"{,.converted}
	
	#
	# Make sure
	#  1. The last line ends with a line end character.
	#  2. We have the right line end character: convert any carriage return (\r) to newline (\n).
	#  3. We remove empty lines.
	#
	printf '\n'     >> "${sampleSheet}.converted"
	sed -i 's/\r/\n/g' "${sampleSheet}.converted"
	sed -i '/^\s*$/d'  "${sampleSheet}.converted"
	#
	# Parse content with Python sanity check script.
	#
	check='failed' # default.
	fileName=$(basename "${sampleSheet}")
	if "${EBROOTNGSMINUTILS}"/"${SCRIPT_NAME}".py --input "${sampleSheet}.converted" --log "${sampleSheet}.converted.log"
	then
		check=$(cat "${sampleSheet}.converted.log")
	fi
	
	if [[ "${check}" == 'OK' ]]
	then
		echo "INFO: Samplesheet is OK, moving ${sampleSheet}.converted to ${sampleSheetsDir}/${fileName}..."
		mv "${sampleSheet}.converted" "${sampleSheetsDir}/${fileName}"
		rm -f "${sampleSheet}"* # cleanup.
	else
		echo "ERROR: Samplesheet ${fileName} is not correct, see ${sampleSheet}.converted.log."
		if [[ -e "${sampleSheet}.converted.log.mailed" ]]
		then
			echo "INFO: Notification was already sent."
		else
			echo "INFO: Trying to send email notification ..."
			#
			# Get email addresses for list of users that should always receive mail.
			#
			declare mailAddress=''
			if [[ -e "${SCR_DIR}/logs/${SCRIPT_NAME}.mailinglist" ]]
			then
				mailAddress="$(cat "${SCR_DIR}/logs/${SCRIPT_NAME}.mailinglist" | tr '\n' ' ')"
			else
				printf '%s\n' "ERROR: ${SCR_DIR}/logs/${SCRIPT_NAME}.mailinglist is missing on $(hostname -s)." \
					| mail -s "Samplesheet is wrong, but we cannot send email to the relevant users." 'helpdesk.gcc.groningen@gmail.com'
			fi
			#
			# Get email address for owner of the samplesheet.
			#
			fileOwner=$(stat -c '%U' "${sampleSheet}" | tr -d '\n')
			mailAddressOwner="$(getent passwd "${fileOwner}" | cut -d ':' -s -f 5)"
			if [[ -z "${mailAddressOwner:-}" ]]
			then
				printf '%s\n' "WARN: We do not have an email address for this user: ${fileOwner}." \
					| mail -s "Samplesheet is wrong on $(hostname -s), but we cannot email the owner." "${mailAddress:-}"
			else
				mailAddress="${mailAddress:-} ${mailAddressOwner:-}"
			fi
			#
			# Prepare message content.
			#
			header="Dear ${fileOwner},"
			body="${SCRIPT_NAME} detected an error when parsing ${sampleSheet} on $(hostname -s): $(<"${sampleSheet}.converted.log")"
			footer='Cheers from the GCC.'
			#
			# Send email to notify users.
			#
			printf '%s\n\n%s\n\n%s\n' "${header}" "${body}" "${footer}" \
				| mail -s "Samplesheet is wrong on $(hostname -s)." "${mailAddress:-}"
			touch "${sampleSheet}.converted.log.mailed"
		fi
	fi
done

echo "INFO: finished processing samplesheets."
