#!/bin/bash

set -u
set -e

function usage () {
        echo "
        Required:
        -i|--input              inputfile (bed format)
	
	Optional:
	-o|--outputfolder	(default: /gcc/resources/b37/intervals)
        "
}

PARSED_OPTIONS=$(getopt -n "$0"  -o i:o: --long "input:,outputfolder:"  -- "$@")

#
# Bad arguments, something has gone wrong with the getopt command.
#
if [ $? -ne 0 ]; then
        usage
        echo "FATAL: Wrong arguments."
        exit 1
fi

eval set -- "$PARSED_OPTIONS"

#
# Now goes through all the options with a case and using shift to analyse 1 argument at a time.
# $1 identifies the first argument, and when we use shift we discard the first argument, so $2 becomes $1 and goes again through the case.
#

while true; do
  case "$1" in
        -i|--input)
                case "$2" in
                "") shift 2 ;;
                *) INPUT=$2 ; shift 2 ;;
            esac ;;
	-o|--outputfolder)
                case "$2" in
                *) OUTPUTFOLDER=$2 ; shift 2 ;;
            esac ;;
        --) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
  esac
done

#
# Check required options were provided.
#
if [[ -z "${INPUT-}" ]]; then
        usage
        echo "FATAL: missing required parameter."
        exit 1
fi

if [[ -z "${OUTPUTFOLDER-}" ]]; then
	OUTPUTFOLDER=/gcc/resources/b37/intervals/
fi

NAME=$(basename $INPUT)
perl /gcc/tools/scripts/create_per_base_intervals.pl -input $INPUT -output $NAME -outputfolder $OUTPUTFOLDER

echo "INPUTNAME:$INPUT"
echo "OUTPUTFOLDER:$OUTPUTFOLDER"
echo "NAME:$NAME"

sort -V -k1 -k2 -k3 $OUTPUTFOLDER/$NAME.per_base.bed | uniq -u > $OUTPUTFOLDER/$NAME.uniq.per_base.bed
rm $OUTPUTFOLDER/$NAME.per_base.bed
