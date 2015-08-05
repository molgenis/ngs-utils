set -e 
set -u

function usage () {
echo "
Arguments
        Required:
        -n|--name	name of FastQ (SequencingStartData_Sequencer_Run_Flowcell, e.g. 150803_SN163_0661_AHKYM5ADXX)

	Optional:
	-t|--tmp	where to write intermediate files (default: /gcc/groups/gaf/tmp03/tmp/)
	-p|--prm	location of the rawdata/ngs directory on the permanent storage (default:/gcc/groups/gaf/prm02/rawdata/ngs)
	-o|--output	outputfolder (default: /gcc/groups/gaf/tmp03/rawdata/ngs/)
"
}


PARSED_OPTIONS=$(getopt -n "$0"  -o n:t:p:o: --long "name:tmp:prm:output"  -- "$@")

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
	 -t|--tmp)
                case "$2" in
                "") shift 2 ;;
                *) TMP=$2 ; shift 2 ;;
            esac ;;
	 -p|--prm)
                case "$2" in
                "") shift 2 ;;
                *) PRM=$2 ; shift 2 ;;
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
if [[ -z "${NAME-}" ]]; then
        usage
        echo "FATAL: missing required parameter."
        exit 1
fi
if [[ -z "${TMP-}" ]]; then
        TMP="/gcc/resources/b37/intervals/"
fi
if [[ -z "${PRM-}" ]]; then
        PRM="/gcc/groups/gaf/prm02/rawdata/ngs"
fi
if [[ -z "${OUTPUT-}" ]]; then
        OUTPUT="/gcc/groups/gaf/tmp03/rawdata/ngs/"
fi

FASTQ=${NAME}
FASTQDIR=${PRM}/${FASTQ}
RAWDATATMP=${OUTPUT}/${FASTQ}_combined/

if [ ! -d ${RAWDATATMP} ]
then    
        mkdir -p ${RAWDATATMP}
	echo "mkdir -p ${RAWDATATMP}"
fi

TMP="/gcc/groups/gaf/tmp03/tmp/"

if [ -f ${TMP}/allBarcodes.txt ]
then
        rm ${TMP}/allBarcodes.txt
fi 

OLDIFS=$IFS
IFS="_"
for i in $(ls -1 ${FASTQDIR}/*.fq.gz)
do

IN="${i}"
set -- "$IN"
declare -a Array=($*)
echo "${Array[8]}" >> ${TMP}/allBarcodes.txt
done

IFS=$OLDIFS

sort -u ${TMP}/allBarcodes.txt > ${TMP}/allUniqBarcodes.txt

while read line
do 
cat ${FASTQDIR}/${FASTQ}_L1_${line}_1.fq.gz ${FASTQDIR}/${FASTQ}_L2_${line}_1.fq.gz > ${RAWDATATMP}/${FASTQ}_combined_${line}_1.fq.gz
echo "${RAWDATATMP}/${FASTQ}_combined_${line}_1.fq.gz done"
cat ${FASTQDIR}/${FASTQ}_L1_${line}_2.fq.gz ${FASTQDIR}/${FASTQ}_L2_${line}_2.fq.gz > ${RAWDATATMP}/${FASTQ}_combined_${line}_2.fq.gz
echo "${RAWDATATMP}/${FASTQ}_combined_${line}_2.fq.gz done"

done<${TMP}/allUniqBarcodes.txt

rm ${TMP}/allUniqBarcodes.txt
rm ${TMP}/allBarcodes.txt

echo "removed intermediate files"
echo "all combined fastq files can be found in ${RAWDATATMP}"


