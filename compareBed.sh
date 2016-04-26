#!/bin/sh

set -u
set -e
underline=`tput smul`
normal=`tput sgr0`
bold=`tput bold`
space="           " # format output

function usage () {
echo "
${bold}This script prints the mismatches between your target panel and two given (RefSeq and CCDS) bed files.

Arguments${normal}
        Required:
        -p|--panel		Your gene panel (BED file)

        Optional:
        -r|--refseq		Path to RefSeq file (default: /apps/data/RefSeq/refseq.gz)
        -c|--ccdsq		Path to CCDS file (default: /apps/data/CCDS/ccds.gz)
        -o|--output		Output folder, will be created if non-existent (default: ./compare/)"
}

PARSED_OPTIONS=$(getopt -n "$0"  -o p:r:c:o: --long "panel:,refseq:ccds:output" -- "$@")

#
# Bad arguments, something has gone wrong with the getopt command.
#
if [ $? -ne 0 ]; then
        usage

        echo "FATAL: Wrong arguments."
        exit 1
fi

eval set -- "$PARSED_OPTIONS"

# Iterate through all the options with a case and using shift to analyse 1 argument at a time.
while true; do
  case "$1" in
        -p|--panel)
                case "$2" in
                "") shift 2 ;;
                *) PANEL=$2 ; shift 2 ;;
            esac ;;
        -r|--refseq)
                case "$2" in
                *) REFSEQ=$2 ; shift 2 ;;
            esac ;;
        -c|--ccds)
                case "$2" in
                *) CCDS=$2 ; shift 2 ;;
            esac ;;
        -o|--output)
                case "$2" in
                *) OUTPUT=$2 ; shift 2 ;;
            esac ;;
        --) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
  esac
done

# Check required options were provided.
if [[ -z "${PANEL-}" ]]; then
	usage
	echo ""
        echo "FATAL: missing required parameter."
        echo ""
	exit 1
fi

# Fill in defaults for parameters that were not required
if [[ -z "${REFSEQ-}" ]]; then
        REFSEQ="/apps/data/RefSeq/refseq.gz"
fi
if [[ -z "${CCDS-}" ]]; then
        CCDS="/apps/data/CCDS/ccds.gz"
fi
if [[ -z "${OUTPUT-}" ]]; then
        OUTPUT="compare/"
fi

# exit if they don't exist
if [ ! -f $PANEL ]; then
        echo "ERROR file not found: $PANEL"
        exit 1
fi

if [ ! -f $REFSEQ ]; then
	echo "ERROR file not found: $REFSEQ"
	exit 1
fi

if [ ! -f $CCDS ]; then
        echo "ERROR file not found: $CCDS"
        exit 1
fi

# Print values of parameters
echo
echo "We use these parameters:"
echo "${space}Your panel:       $PANEL"
echo "${space}Refseq:           $REFSEQ"
echo "${space}CCDS:             $CCDS"
echo "${space}Output path:      $OUTPUT"
echo

echo "Start working..."
echo "${space}Create $OUTPUT if non-existent"
mkdir -p ${OUTPUT}

echo "${space}Copy (and unzip) your files to $OUTPUT..."
PANEL_FILE_NAME=$(basename $PANEL)
REFSEQ_FILE_NAME=refseq
CCDS_FILE_NAME=ccds
cp   $PANEL    ${OUTPUT}/${PANEL_FILE_NAME}.original
zcat $REFSEQ > ${OUTPUT}/${REFSEQ_FILE_NAME}.original
zcat $CCDS   > ${OUTPUT}/${CCDS_FILE_NAME}.original

cd $OUTPUT

echo "${space}Remove prefix 'chr' from first column..."
	declare -a files=("$PANEL_FILE_NAME" "$REFSEQ_FILE_NAME" "$CCDS_FILE_NAME")
	for file in ${files[@]}; do
		awk -F '\t' 'BEGIN {OFS=FS} { print $1,$2,$3,$4 }' ${file}.original | sed '/^chr/s/^...//' > $file
	done

echo "${space}Load bedtools..."
#	module load bedtools/2.22.0
	module load BEDTools/2.23.0-goolf-1.7.20
	

	module list

function create_coverage_files {
	bed1=$1
	bed2=$2
	bed1vs2=$1_vs_$2
	bed2vs1=$2_vs_$1
	
	bigspace=${space}${space}

	echo "${bigspace}Create ${bed1vs2}..."
	bedtools coverage -b $bed1 -a $bed2 > ${bed1vs2}
	awk '{ if ($8 < 1) print $0 }' $bed1vs2 > ${bed1vs2}.notFullyCovered
	sort -V -k 1,2 ${bed1vs2}.notFullyCovered > ${bed1vs2}.notFullyCovered.sorted
	bedtools merge -i ${bed1vs2}.notFullyCovered.sorted > ${bed1vs2}.notFullyCovered.sorted.merged

	echo "${bigspace}Create ${bed2vs1}..."
	bedtools coverage -a $bed1 -b $bed2 > ${bed2vs1}

	echo "${bigspace}Get N?_* names of regions covered at least one base..."
	awk '{ if ($8 > 0) {split($4,a,"_exon_"); print a[1]} }' ${bed2vs1} | sort -u > ${bed2vs1}.regionNames
	while read line; do
		grep ${line}_ ${bed2vs1} >> ${bed2vs1}.touchedRegions
	done < ${bed2vs1}.regionNames

	awk '{ if ($8 < 1) print $0 }' ${bed2vs1}.touchedRegions > ${bed2vs1}.touchedRegions.notFullyCovered
	sort -V -k 1,2 ${bed2vs1}.touchedRegions.notFullyCovered > ${bed2vs1}.touchedRegions.notFullyCovered.sorted
	bedtools merge -i ${bed2vs1}.touchedRegions.notFullyCovered.sorted > ${bed2vs1}.touchedRegions.notFullyCovered.sorted.merged
}

echo "${space}Determine overlap your panel with REFSEQ..."
	create_coverage_files $PANEL_FILE_NAME $REFSEQ_FILE_NAME

echo "${space}Determine overlap your panel with CCDS..."
	create_coverage_files $PANEL_FILE_NAME $CCDS_FILE_NAME
echo
echo Done! Please find your results in $OUTPUT
echo
