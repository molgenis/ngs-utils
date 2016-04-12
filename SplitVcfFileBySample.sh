WHO=`whoami`
bold=`tput bold`
normal=`tput sgr0`

function usage {
echo "
${bold}This tool is to split a vcf by sample${normal}

	Required:
	-i|--input		vcf file (complete path) 

        Optional:
	-o|--outputfolder	specify the outputfolder (default:/home/${WHO}/SplitVcf --> folder will be created if not exists )
	-r|--reference    	which reference file is used (default:/apps/data/GRC/GRCh37/Homo_sapiens.GRCh37.GATK.illumina.fasta)
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


if [ ! -d /home/${WHO}/SplitVcf ]
then
	mkdir -p /home//${WHO}/SplitVcf
fi

if [[ -z "${INPUT-}" ]]; then
        usage
        echo "FATAL: missing required parameter."
        exit 1
fi
if [[ -z "${REFERENCE-}" ]]; then
        REFERENCE="/apps/data/GRC/GRCh37/Homo_sapiens.GRCh37.GATK.illumina.fasta"
fi
if [[ -z "${OUTPUT-}" ]]; then
        OUTPUT="/home/${WHO}/SplitVcf"
fi

module load GATK

IETS=$(awk '{if($1=="#CHROM")print $0}' ${INPUT})
teller=0

allsamples=()
for i in $IETS
do 
	if [ ${teller} -gt 8 ]
	then
		allsamples+=(${i}) 
	fi
	(( teller=$((teller+1)) ))

done

echo "length: ${#allsamples[@]}"

for i in ${allsamples[@]}
do
	java -Xmx2g -jar /apps/software/GATK/3.4-46-Java-1.7.0_80/GenomeAnalysisTK.jar \
	-R ${REFERENCE} \
	-T SelectVariants \
	--variant ${INPUT} \
	-o ${OUTPUT}/${i}.splitted.vcf \
	-sn ${i}
done
