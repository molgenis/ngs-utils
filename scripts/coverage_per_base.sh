set -u
set -e

function usage () {
        echo "
        This script is making a coverage_per_base bed file and interval_list

        Required:
        -n|--name              BED name (without extension)

        Optional:
        -o|--outputfolder       path to outputfolder (default: /gcc/resources/b37/intervals)
        -e|--extension          extension to the bed file (default: _b37_human_g1k_v37)
        -r|-reference           Which reference file is used (default: /gcc/resources/b37/indices/human_g1k_v37.dict)
       "
}

PARSED_OPTIONS=$(getopt -n "$0"  -o n:o:e:r: --long "name:,outputfolder:extension:reference:"  -- "$@")

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
        -n|--name)
                case "$2" in
                "") shift 2 ;;
                *) NAME=$2 ; shift 2 ;;
            esac ;;
        -o|--outputfolder)
                case "$2" in
                *) OUTPUTFOLDER=$2 ; shift 2 ;;
            esac ;;
        -e|--extension)
                case "$2" in
                *) EXTENSION=$2 ; shift 2 ;;
            esac ;;
        -r|--reference)
                case "$2" in
                *) REFERENCE=$2 ; shift 2 ;;
            esac ;;
        --) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
  esac
done

#
# Check required options were provided.
if [[ -z "${NAME-}" ]]; then
        usage
        echo "FATAL: missing required parameter."
        exit 1
fi

if [[ -z "${OUTPUTFOLDER-}" ]]; then
        OUTPUTFOLDER=/gcc/resources/b37/intervals/
fi

if [[ -z "${EXTENSION-}" ]]; then
        EXTENSION="_b37_human_g1k_v37"
fi
if [[ -z "${REFERENCE-}" ]]; then
        REFERENCE="/gcc/resources/b37/indices/human_g1k_v37.dict"
fi

BAITS=${OUTPUTFOLDER}/${NAME}_baits${EXTENSION}
EXONS=${OUTPUTFOLDER}/${NAME}_exons${EXTENSION}

TMP="/gcc/groups/gcc/tmp01/tmp"
if [ ! -f ${BAITS}.uniq.per_base.bed ]
then
	perl /gcc/tools/scripts/create_per_base_intervals.pl -input ${BAITS}.bed -output ${NAME}_baits${EXTENSION} -outputfolder $TMP

        sort -V -k1 -k2 -k3 ${TMP}/${NAME}_baits${EXTENSION}.per_base.bed | uniq -u > ${BAITS}.uniq.per_base.bed
        rm ${TMP}/${NAME}_baits${EXTENSION}.per_base.bed

        echo "intervals per base done: ${BAITS}.uniq.per_base.bed"
else
        echo "${BAITS}.uniq.per_base.bed already exists, skipped!"
fi

#make interval_list coverage per base
if [ ! -f ${BAITS}.uniq.per_base.interval_list ] 
	head -85 ${OUTPUTFOLDER}/Agilent_SureSelect_Human_Inherited_Disease_V1_S0684402_UMCG_baits_b37_human_g1k_v37.interval_list > ${BAITS}.uniq.per_base.interval_list
	cat ${BAITS}.uniq.per_base.bed >> ${BAITS}.uniq.per_base.interval_list
else
	echo "coverage per base skipped"
fi
fi
