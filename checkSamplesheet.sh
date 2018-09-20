#!/bin/bash

set -eu

declare baseDir='/groups/umcg-gd/scr01/'
declare sampleSheetsDir="${baseDir}"'Samplesheets/'

SCRIPT_NAME="$(basename ${0})"
SCRIPT_NAME="${SCRIPT_NAME%.*sh}"

echo "INFO: processing samplesheets from ${sampleSheetsDir}/new/..."
for sampleSheet in $(ls -1 "${sampleSheetsDir}/new/"*'.csv')
do
	#
	# Make sure
	#  1. The last line ends with a line end character.
	#  2. We have the right line end character: convert any carriage return (\r) to newline (\n).
	#  3. We remove empty lines.
	#
	cp "${sampleSheet}"{,.converted} \
		&& printf '\n'     >> "${sampleSheet}.converted" \
		&& sed -i 's/\r/\n/g' "${sampleSheet}.converted" \
		&& sed -i '/^\s*$/d'  "${sampleSheet}.converted" \
		&& mv -f "${sampleSheet}"{.converted,}
	#
	# Parse content with Python sanity check script.
	#
	"${sampleSheetsDir}"/checkSampleSheet_v2.py --input "${sampleSheet}" --logfile "${sampleSheet}.logger"
	filename=$(basename "${sampleSheet}")
	check=$(cat "${sampleSheet}.logger")
	if [[ "${check}" == "OK" ]]
	then
		echo "INFO: Samplesheet is OK, moving ${sampleSheet} to ${sampleSheetsDir}..."
		mv "${sampleSheet}" "${sampleSheetsDir}"
		rm -f "${sampleSheet}.logger.mailed"
		rm -f "${sampleSheet}.logger"
		rm -f "${sampleSheet}.logger.text"
	else
		echo "ERROR: Samplesheet ${filename} is not correct, see logger."
		if [[ -e "${sampleSheet}.logger.mailed" ]]
		then
			echo "INFO: Notification was already sent."
		else
			echo "INFO: Trying to send email notification ..."
			#
			# Get list of users that should always receive mail
			#
			declare mailAddress=''
			if [[ -e "${baseDir}/logs/${SCRIPT_NAME}.mailinglist" ]]
			then
				mailAddress="$(cat "${baseDir}/logs/${SCRIPT_NAME}.mailinglist" | tr '\n' ' ')"
			else
				echo -e "ERROR: ${baseDir}/logs/${SCRIPT_NAME}.mailinglist is missing on $(hostname -s)\n" \
					| mail -s "Samplesheet is wrong, but we cannot send email to the relevant users."
			fi
			fileOwner=$(stat -c "%U" "${sampleSheet}" | tr -d '\n')
			mailAddressOwner="$(getent passwd "${fileOwner}" | cut -d ':' -s -f 5)"
			if [[ -z "${mailAddressOwner:-}" ]]
			then
				echo -e "WARN: We do not have an email address of this user: ${fileOwner}\n" \
					| mail -s "Samplesheet is wrong on $(hostname -s), but we cannot email the owner." "${mailAddress:-}"
			else
				
			fi
			mailAddress="${mailAddress:-} ${mailAddressOwner:-}
			
			echo -e "Dear ${fileOwner},\n" > "${sampleSheet}.logger.text"
			cat "${sampleSheet}.logger"   >> "${sampleSheet}.logger.text"
			echo -e "\n\nCheers from GCC" >> "${sampleSheet}.logger.text"
			
			cat "${sampleSheet}.logger.text" \
				| mail -s "Samplesheet is wrong on $(hostname -s)" "${mailAddress:-}"
			touch "${sampleSheet}.logger.mailed"
			rm -f "${sampleSheet}.logger.text"
		fi
	fi
done
