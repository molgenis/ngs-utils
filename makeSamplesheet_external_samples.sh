#!/bin/bash

set -e
set -u


function showHelp() {
    #
    # Display commandline help on STDOUT.
    #
    cat <<EOH
===============================================================================================================
Usage:
        $(basename $0) OPTIONS
        
Options:    
        --h  Show this help.

        Required:
        -p  projectName
        -d  sequencingStartDate (YYMMDD)
        -o  externalSampleDirectory
        -f  FastQ_filename_regex  Regex pattern to find the FastQ files that must be renamed.
            Note: the pattern must be single quoted to prevent expansion by the shell.
            E.g. -f 'HTVTYBBXX_103373-003*'

        Optional:
        -w  workDir (default is this directory)
        -b  barcodeType (default: RPI. AGI,LEX, NEB enz.)
        -c  capturingKit (default:None, name of capturingKit (e.g. Agilent\/ONCO_v3)
        -k  prepKit (default:Exceptional Prep kit. Important for RNA: lexogen, TruSeq, RiboZero)
        -g  species (default:homo_sapiens|GRCh37, if other species, supply single quoted 'homo_sapiens|GRCh37')
        -s  sampleType (default:RNA, DNA/RNA)
        -t  seqType (default:PE. otherwise specify SR)
        -m  pairedMateNameOne (default= _R1_, difference between mateNames, e.g. R1 for mate 1 and R2 for mate 2)
                N.B. please also specify the character before and after the 1, e.g. _1.)
        
        
Output will be written in workDir with the name: {projectName}.csv
===============================================================================================================
EOH
    trap - EXIT
    exit 0
}
sequencingStartdate=
while getopts "w:s:k:t:g:p:c:d:f:b:m:o:h" opt;
do

    case $opt in h)showHelp;; 
                 w)workDir="${OPTARG}";; 
                 s)sampleType="${OPTARG}";; 
                 k)prepKit="${OPTARG}";; 
                 t)seqType="${OPTARG}";; 
                 g)species="${OPTARG}";; 
                 p)projectName="${OPTARG}";; 
                 c)capturingKit="${OPTARG}";; 
                 d)sequencingStartDate="${OPTARG}";; 
                 f)FastQ_filename_regex="${OPTARG}";; 
                 b)barcodeType="${OPTARG}";; 
                 m)pairedMateNameOne="${OPTARG}";;
                 o)externalSampleDirectory="${OPTARG}";
    esac
done



if [[ -z "${sampleType:-}" ]]
then
    sampleType="RNA"
fi

if [[ -z "${capturingKit:-}" ]]
then
    capturingKit="None"
fi

if [[ -z "${projectName:-}" ]]
then
    echo -e '\nERROR: Must specify a projectName\n'
    exit 1
fi

if [[ -z "${seqType:-}" ]]
then
    seqType="PE"
fi

if [[ -z "${prepKit:-}" ]]
then
    prepKit="'Exceptional Prep kit'"
fi

if [[ -z "${workDir:-}" ]]
then
    workDir=$(pwd)
fi

if [[ -z "${species:-}" ]]
then
    species="homo_sapiens|GRCh37"
fi

if [[ -z "${barcodeType:-}" ]]
then
    barcodeType="RPI"
fi

if [[ -z "${pairedMateNameOne:-}" ]]
then
    pairedMateNameOne="_R1_"
fi

if [[ -z "${externalSampleDirectory:-}" ]]
then
    echo -e '\n ERROR: must specify externalSampleDirectory'
    exit 1
else 
    externalSampleDirectoryPath=$(pwd)/"${externalSampleDirectory}"
fi

if [[ -z "${fastqFilePattern:-}" ]]
then
    echo -e "ERROR: must specify fastQfilePattern"
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
##
### Function to retrieve from the supplied fastQ files: Sequencer, flowcell, run, lane, barcode.
##
#


function _makeSampleInfo() {
    local _fastqPath="${1}"
    local _sequencingStartDate="${2}"
    echo "INFO: Processing FastQ ${_fastqPath}"
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
    
    local _newFastqFileInfo=",${_sequencingStartDate},${_sequencer},${_run},${_flowcell},L${_lane},${_barcodes}\n"
    echo "INFO: new FastQ file information: ${_sequencingStartDate},${_sequencer},${_run},${_flowcell},L${_lane},${_barcodes}"
    printf "${_newFastqFileInfo}" >> "${workDir}/${projectName}.csv"
    }


#
# Make header for your sample sheet 
#

printf "externalSampleID,externalFastQ_1,externalFastQ_2,project,capturingKit,sampleType,seqType,prepKit,species,Gender,arrayFile,sequencingStartDate,sequencer,run,flowcell,lane,barcode\n" > "${workDir}/${projectName}.csv"

#
# Process FastQ files, and add other necessary information
#

for FastQ in $(ls "${externalSampleDirectoryPath}/"*"${pairedMateNameOne}"*".gz")
do
    if [[ "${seqType}" == "PE" ]]
    then
        fileName=${FastQ%%.*}
        baseNameFile=$(basename ${fileName})
        withoutExtension=${baseNameFile%%.*}
        printf "${withoutExtension}" >> "${workDir}/${projectName}.csv" ##externalSampleID
        printf ",${FastQ}," >> "${workDir}/${projectName}.csv" ## externalFastQ1 

        pairedMateNameTwo=$(echo "${pairedMateNameOne}" | sed 's|1|2|')
        cleanPairedMateNameOne=$(echo $pairedMateNameOne | sed 's|\.|\\.|')
        cleanPairedMateNameTwo=$(echo $pairedMateNameTwo | sed 's|\.|\\.|')

        printf "${FastQ}" | sed -r "s|${cleanPairedMateNameOne}|${cleanPairedMateNameTwo}|g"  >> "${workDir}/${projectName}.csv" ## externalFastQ2
        printf ",${projectName}" >> "${workDir}/${projectName}.csv" ## project
        printf ",${capturingKit}" >> "${workDir}/${projectName}.csv" ## capturingKit
        printf ",${sampleType}" >> "${workDir}/${projectName}.csv" ## sampleType
        printf ",${seqType}" >> "${workDir}/${projectName}.csv" ## seqType
        printf ",${prepKit}" >> "${workDir}/${projectName}.csv" ## prepKit
        printf ",${species}" >> "${workDir}/${projectName}.csv" ## species
        printf "," >> "${workDir}/${projectName}.csv"  ## Gender
        printf "," >> "${workDir}/${projectName}.csv" ## arrayFile
        _makeSampleInfo "${FastQ}" "${sequencingStartDate}" ## SequencingStartDate,sequencer,run,flowcell,lane,barcode
    elif [[ "${seqType}" == "SR" ]]
    then
        fileName=${FastQ%%.*}
        baseNameFile=$(basename ${fileName})
        withoutExtension=${baseNameFile%%.*}
        printf "${withoutExtension}" >> "${workDir}/${projectName}.csv" ##externalSampleID
        printf ",${FastQ}" >> "${workDir}/${projectName}.csv" ## externalFastQ1 
        printf ","  >> "${workDir}/${projectName}.csv" ## externalFastQ2
        printf ",${projectName}" >> "${workDir}/${projectName}.csv" ## project
        printf ",${capturingKit}" >> "${workDir}/${projectName}.csv" ## capturingKit
        printf ",${sampleType}" >> "${workDir}/${projectName}.csv" ## sampleType
        printf ",${seqType}" >> "${workDir}/${projectName}.csv" ## seqType
        printf ",${prepKit}" >> "${workDir}/${projectName}.csv" ## prepKit
        printf ",${species}" >> "${workDir}/${projectName}.csv" ## species
        printf "," >> "${workDir}/${projectName}.csv"  ## Gender
        printf "," >> "${workDir}/${projectName}.csv" ## arrayFile
        _makeSampleInfo "${FastQ}" "${sequencingStartDate}" ## SequencingStartDate,sequencer,run,flowcell,lane,barcode
    fi
done

echo "Samplesheet is finished!"






