#!/bin/bash

set -e 
set -u

underline=`tput smul`
normal=`tput sgr0`
bold=`tput bold`

function usage () {
echo "
${bold}
echo 'to run: sh bamout.sh \${region} \${bam}'
echo "e.g. sh bamout 1:1000-2000 /path/to/bam"

${bold}Arguments${normal}

        Required:
        -p|--region             region of interest (e.g. 1:100-120)
	-i|--bam		bam input file 

        Optional:
	-w|--workingdir		default (this dir)
        -r|--reference          Which reference file is used (default: /apps/data/1000G/phase1/human_g1k_v37_phiX.fasta)
	-g|--gender		When the region is chromosome X then a gender should be specified (Male/Female)"
}

module load ngs-utils
PARSED_OPTIONS=$(getopt -n "$0"  -o p:i:r:g:w: --long "region:bam:reference:gender:workingdir:"  -- "$@")

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
	-p|--region)
                case "$2" in
                "") shift 2 ;;
                *) REGION=$2 ; shift 2 ;;
            esac ;;
        -r|--reference)
                case "$2" in
                *) REFERENCE=$2 ; shift 2 ;;
            esac ;;
	-g|--gender)
                case "$2" in
                *) GENDER=$2 ; shift 2 ;;
            esac ;;
	-w|--workingdir)
                case "$2" in
                *) WORKDIR=$2 ; shift 2 ;;
            esac ;;
        -i|--bam)
                case "$2" in
                *) INPUT=$2 ; shift 2 ;;
            esac ;;
        --) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
  esac
done

ploidy="2"

#
# Check required options were provided.
if [[ -z "${REGION-}" ]]; then
        usage
        exit 1
fi
if [[ -z "${INPUT-}" ]]; then
        usage
	exit 1
fi
if [[ -z "${REFERENCE-}" ]]; then
        REFERENCE="/apps/data/1000G/phase1/human_g1k_v37_phiX.fasta"
fi

if [[ -z "${WORKDIR-}" ]]; then
        WORKDIR=$(pwd)
fi


if [[ $REGION == *"X"* || $REGION == *"x"* ]]
then
	if [[ -z "${GENDER-}" ]]; then
		echo "No gender supplied, please fill in a gender (Male or Female) when using chromosome X region"
        	usage
		exit 1
	fi

	if [ $GENDER == "Male" ]
	then
		ploidy="1"
	elif [ $GENDER == "Female" ]
	then
		ploidy="2"
	else
		echo "Gender is not Female or Male, exiting"
		exit 1
	fi	
fi


ml GATK/3.7-Java-1.8.0_74 

myregion=$(echo "${REGION}" | tr : _)

BAS=$(basename ${INPUT})

java -XX:ParallelGCThreads=2 -Xmx12g -jar \
${EBROOTGATK}/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R ${REFERENCE} \
-I ${INPUT} \
-bamout ${WORKDIR}/${BAS}.${myregion}.bam \
-variant_index_type LINEAR \
-variant_index_parameter 128000 \
-o ${WORKDIR}/"${BAS}".vcf \
-L "${REGION}" \
--emitRefConfidence GVCF \
-ploidy ${ploidy}
