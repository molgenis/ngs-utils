set -eu
count=0

function showHelp() {
        #
        # Display commandline help on STDOUT.
        #
        cat <<EOH
===============================================================================================================
===============================================================================================================

Usage:
	$(basename $0) OPTIONS
Options:
	-h   Show this help.

   required:
	-o   filetype of file one: OPENARRAY (default) or VCF
	-t   filetype of file two: OPENARRAY (default) or VCF
	-s   samplesheet; format: projectname1\tprojectname2\tdnanumber1\tdnanumber2\tsampleprocessstepid (default: ./samplesheet.txt)
	-p   path to files (file 1)
	-c   path to files (file 2)
===============================================================================================================
EOH
        trap - EXIT
        exit 0
}

echo -e"example Usage:\n\tbash makeSamplesheet.sh -o OPENARRAY -t OPENARRAY -p /groups/umcg-atd/tmp07/concordance/ngs/OA99/ -c /groups/umcg-atd/tmp07/concordance/ngs/OA99/compareWith/"
sleep 2
while getopts "ho:t:s:c:p:" opt;
do
	case $opt in h)showHelp;; o)filetypeone="${OPTARG}";; t)filetypetwo="${OPTARG}";; s)samplesheet="${OPTARG}";; p)f1concordancePath="${OPTARG}";; c)f2concordancePath="${OPTARG}";; 
esac
done

if [[ -z "${filetypeone:-}" &&  -z "${filetypetwo:-}" ]]
then
	filetypeone="OPENARRAY"
	filetypetwo="VCF"
fi

if [[ -z "${f1concordancePath:-}" &&  -z "${f2concordancePath:-}" ]]
then
	f1concordancePath='/groups/umcg-atd/tmp07/concordance/ngs/'
	f1concordancePath="${f2concordancePath}"
	filetypeone="OPENARRAY"
	filetypetwo="VCF"
fi



if [[ -z "${samplesheet:-}" ]]
then
	samplesheet="samplesheet.txt"
fi


f1refGenome="GRCh37"
f2refGenome="GRCh37"
f1Extension=''
f2Extension=''

echo "filetypeone=${filetypeone}"
echo "filetypetwo=${filetypetwo}"


if [[ "${filetypeone}" == "OPENARRAY" ]]
then
	f1Extension='.oarray.txt'
elif [[ "${filetypeone}" == "VCF" ]]
then
	f1Extension='.concordanceCheckCalls.vcf'
else
	echo "FATAL: unknown filetype: ${filetypeone}"
	exit 1
fi

if [[ "${filetypetwo}" == "OPENARRAY" ]]
then
	f2Extension='.oarray.txt'
elif [[ "${filetypetwo}" == "VCF" ]]
then
	f2Extension='.concordanceCheckCalls.vcf'
else
	echo "FATAL: unknown filetype: ${filetypetwo}"
	exit 1
fi

count=0
while read line 
do
	if [[ "${count}" == 0 ]]
	then
		echo -e "data1Id\tdata2Id\tlocation1\tlocation2\tfileType1\tfileType2\tbuild1\tbuild2\tproject1\tproject2\tfileprefix\tprocessStepId" > newSamplesheetje.txt
		count=1
	fi

	f1project=$(echo "${line}" | awk 'BEGIN {FS="\t"}{print $1}')
	f2project=$(echo "${line}" | awk 'BEGIN {FS="\t"}{print $2}')
	f1sample=$(echo "${line}" | awk 'BEGIN {FS="\t"}{print $3}')
	f2sample=$(echo "${line}" | awk 'BEGIN {FS="\t"}{print $4}')
	sampleProcess=$(echo "${line}" | awk 'BEGIN {FS="\t"}{print $5}')
	f1Sample=$(ls "${f1concordancePath}/"*"${f1sample}"*"${f1Extension}")
	f2Sample=$(ls "${f2concordancePath}/"*"${f2sample}"*"${f2Extension}")
	f1Id=$(basename "${f1Sample}" "${f1Extension}")
	f2Id=$(basename "${f2Sample}" "${f2Extension}")
	filePrefix="${sampleProcess}_${f1project}_${f1sample}_${f2project}_${f2sample}"

	echo -e "${f1Id}\t${f2Id}\t${f1Sample}\t${f2Sample}\t${filetypeone}\t${filetypetwo}\t${f1refGenome}\t${f2refGenome}\t${f1project}\t${f2project}\t${filePrefix}\t${sampleProcess}" >> 'newSamplesheetje.txt'

done<"${samplesheet}"

mkdir -p './output'
input='newSamplesheetje.txt'
header=$(head -n 1 "${input}")

#skip header
tail -n +2 "$input" | while IFS=$'\t' read -r data1Id data2Id location1 location2 fileType1 fileType2 build1 build2 project1 project2 fileprefix processStepId; do
	# maak een bestandsnaam op basis van de prefix
	filename="${fileprefix}.sampleId.txt"

    # schrijf header en de regel naar het bestand
	{
		echo "$header"
		echo -e "${data1Id}\t${data2Id}\t${location1}\t${location2}\t${fileType1}\t${fileType2}\t${build1}\t${build2}\t${project1}\t${project2}\t${fileprefix}\t${processStepId}"
} > "output/$filename"


