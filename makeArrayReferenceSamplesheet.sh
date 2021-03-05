#!/bin/bash


set -e
set -u

# nadenken over dubbele samplesheets op prm05/prm06 waardoor samples dubbel in de referentie samplesheet komen
# voor het maken van de referentie samplesheet gaan we de archief mappen door op prm06 misschien dat dit voor de meest recente samples niet optimaal is
# outputten hoeveel samples er in de samplesheet zitten

function showHelp() {
        #
        # Display commandline help on STDOUT.
        #
        cat <<EOH
======================================================================================================================
Script generates a big samplesheet and a temporary directory where all idats of those samples are stored so diagnostics can use these for making the array reference in genomestudio.
Usage:
        $(basename "${0}") OPTIONS
Options:
        -h      Show this help.
        -s      list of samples which needs to be added to one big samplesheet in csv format.
        -d      The directory where the samplesheet and the idats are stored for use in genomestudio.
===============================================================================================================
EOH
        trap - EXIT
        exit 0
}

#
##
### Main.
##
#

#
# Get commandline arguments.
#

while getopts "d:s:h" opt
do
        case "${opt}" in
                h)
                        showHelp
                        ;;
                s)
                        sampleList="${OPTARG}"
                        ;;
                d)
                        workDirectory="${OPTARG}"
                        ;;
                *)
                        "Unhandled option. Try $(basename "${0}") -h for help."
                        ;;
        esac
done

#
# Check commandline options.
#
if [[ -z "${sampleList:-}" ]]
then
        echo 'Must specify a sample list with -s.'
        exit 0
fi

if [[ -z "${workDirectory:-}" ]]
then
        echo 'Must specify a directory where the idat files and samplesheet should be stored -d.'
        exit 0
fi

mkdir -p ${workDirectory}

datum=$(date '+%Y%m%d')
echo "${datum}"
outputFile="${workDirectory}/reference_samplesheet_${datum}.csv"
rm -rf "${outputFile}"
touch "${outputFile}"
input="${sampleList}"


# Making the reference samplesheet

while read line
 do
        echo "Making reference samplesheet for these samples:"
        echo "Warning: may contain duplicate samples, because samplesheets are sometimes stored on both prm05 and prm06 so please check the samplesheet before using it in GenomeStudio"
        echo "${line}"
        if grep "${line}" /groups/umcg-gap/prm0*/Samplesheets/archive/*.csv
        then
                grep "${line}" /groups/umcg-gap/prm0*/Samplesheets/archive/*.csv >> "${outputFile}"
        else
                echo "cannot find ${line}"
        fi
done < "${input}"

# Getting the barcode and position and adding them to a file

echo "test + ${outputFile}"
awk 'BEGIN {FS=","}{print $7 "_" $8}' "${outputFile}" > "${workDirectory}/barcodes.txt"

# Looping through the list of barcodes and positions

while read line
do
        echo "${line}"
        col1=$(echo "${line}" | awk 'BEGIN {FS="_"}{print $1}')
        col2=$(echo "${line}" | awk 'BEGIN {FS="_"}{print $2}')
        echo "Copying data to location for making reference in GenomeStudio:"
        echo "location is:"
        echo "slide is : ${col1}"
        echo "Copying data /groups/umcg-gap/prm0*/rawdata/array/IDAT/${col1}/${col1}_${col2}*.idat"
        ls -latrh /groups/umcg-gap/prm0*/rawdata/array/IDAT/"${col1}/${col1}_${col2}"*.idat
        rsync -av /groups/umcg-gap/prm0*/rawdata/array/IDAT/"${col1}/${col1}_${col2}"*.idat "${workDirectory}/scanData"
done <"${workDirectory}/barcodes.txt"
