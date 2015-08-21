WHO=`whoami`
bold=`tput bold`
normal=`tput sgr0`

function usage {
echo "
${bold}This tool is to split a vcf by sample${normal}

	Required:
	-i|--input		vcf file (complete path) 

        Optional:
	-o|--outputfolder	specify the outputfolder (default:/gcc/groups/gaf/tmp03/${WHO}/SplitVcf --> folder will be created if not exists )
	-r|--reference    	which reference file is used (default:/gcc/resources/b37/indices/Homo_sapiens.GRCh37.GATK.illumina.fasta)
"
}

PARSED_OPTIONS=$(getopt -n "$0"  -o i:r:o --long "inputfile:reference:outputfolder:"  -- "$@")

#
# Bad arguments, something has gone wrong with the getopt command.
#
if [ $? -ne 0 ]; then
        usage
        echo "FATAL: Wrong arguments."
        exit 1
fi

eval set -- "$PARSED_OPTIONS"

while true; do
  case "$1" in
        -i|--input)
                case "$2" in
                "") shift 2 ;;
                *) INPUT=$2 ; shift 2 ;;
            esac ;;
	-r|--reference)
                case "$2" in
                "") shift 2 ;;
                *) REFERENCE=$2 ; shift 2 ;;
            esac ;;
	-o|--outputfolder)
                case "$2" in
                "") shift 2 ;;
                *) OUTPUT=$2 ; shift 2 ;;
            esac ;;
	--) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
  esac
done


if [ ! -d /gcc/groups/gaf/tmp03/${WHO}/SplitVcf ]
then
	mkdir -p /gcc/groups/gaf/tmp03/${WHO}/SplitVcf
fi

if [[ -z "${INPUT-}" ]]; then
        usage
        echo "FATAL: missing required parameter."
        exit 1
fi
if [[ -z "${REFERENCE-}" ]]; then
        REFERENCE="/gcc/resources/b37/indices/Homo_sapiens.GRCh37.GATK.illumina.fasta"
fi
if [[ -z "${OUTPUT-}" ]]; then
        OUTPUT="/gcc/groups/gaf/tmp03/${WHO}/SplitVcf"
fi

module load GATK/3.4-0-g7e26428

IETS=$(awk '{if($1=="#CHROM")print $0}' ${INPUT})
teller=0

allsamples=()
for i in $IETS
do 
	if [ ${teller} -gt 9 ]
	then
		allsamples+=(${i}) 
	fi
	(( teller++ ))

done

echo "length: ${#allsamples[@]}"

for i in ${allsamples[@]}
do
	java -Xmx2g -jar /gcc/tools/GATK-3.4-0-g7e26428/GenomeAnalysisTK.jar \
	-R ${REFERENCE} \
	-T SelectVariants \
	--variant ${INPUT} \
	-o ${OUTPUT}/${i}.splitted \
	-sn ${i}
done
