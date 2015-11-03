set -e 
set -u

function usage () {
echo "
To run this script, please run the following command:
sh extractRuntimes.sh -i INPUTFOLDER -o OUTPUTFOLDER
required arguments:

	-i|--input	Path to inputfolder e.g. /groups/umcg-gaf/tmp04/projects/projectX/run01/jobs/  (file will always be molgenis.bookkeeping.log)	
	-o|--output	Path to outputfolder (also the tmp files will be written to this folder (in a seperate folder called TMP))
"
}
module load ngs-utils
PARSED_OPTIONS=$(getopt -n "$0"  -o i:o: --long "input:output"  -- "$@")

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
	-o|--output)
                case "$2" in
                "") shift 2 ;;
                *) OUTPUT=$2 ; shift 2 ;;
            esac ;;
        --) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
  esac
done

# 
# Check required options were provided.
if [[ -z "${INPUT-}" ]]; then
        usage
        exit 1
fi
if [[ -z "${OUTPUT-}" ]]; then
        usage
        exit 1
fi

tmpFolder=${OUTPUT}/TMP
mkdir -p ${tmpFolder}
awk '{print $8"\t"$5}' ${INPUT}/molgenis.bookkeeping.log > ${tmpFolder}/cuttedRuntimes.txt
echo ${INPUT}/molgenis.bookkeeping.log

awk 'BEGIN { FS="_" }; {print $1"_"$2}' ${tmpFolder}/cuttedRuntimes.txt > ${tmpFolder}/runtimeNames.txt
awk 'BEGIN { FS="\t" }; {print $2}' ${tmpFolder}/cuttedRuntimes.txt > ${tmpFolder}/runtimeTimes.txt

paste -d "\t" ${tmpFolder}/runtimeNames.txt ${tmpFolder}/runtimeTimes.txt | sort -V -k1 > ${tmpFolder}/combinedSortedRuntimeNameAndTime.txt

python count.py --input ${tmpFolder}/combinedSortedRuntimeNameAndTime.txt --output ${tmpFolder}/calculationRuntimes.txt

cat ${tmpFolder}/calculationRuntimes.txt | sort -V -k1 | sed '1istep\tmax\tmean'
echo -e "\nOutput is also written to: ${OUTPUT}/runtimes.txt"

cat combinedSortedRuntimeNameAndTime.txt | sort -V -k1 | sed '1istep\tmax\tmean' > ${OUTPUT}/runtimes.txt
