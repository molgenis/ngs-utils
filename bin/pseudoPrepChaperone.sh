
## This script works with this input format:
## UMCGnr	DNA	new_id	file
## The projectname should already be splitted in the file!!!
## e.g. originally it says _QXTR888_ but should already been renamed to _QXTR_888_
## e.g.2. original HSR825A --> HSR_825A_
## DUPLICATES ARE COPIED, running the pseudo_vcf script will error on the duplicates, just remove it from the inputs on the diagnostic cluster and rerun script till success


mappingfile=$1 ## file described as above UMCGnr DNA new_id file
mappingDestinationFolder=$2 ## foldername in /groups/umcg-atd/tmp06/pseudo/

if [[ -z "${mappingfile:-}" || -z "${mappingDestinationFolder:-}" ]]
then
	echo "we need a mappingfile and mappingDestinationFolder as inputs"
	exit 1
fi

awk 'BEGIN {FS="_"} {print $7"_"$8}' "${mappingfile}" > "${mappingfile%.*}.project.txt"

paste -d "\t" "${mappingfile}" "${mappingfile%.*}.project.txt" > "${mappingfile%.*}.merged.txt"
ssh cf-porch+copperfist "mkdir -p /groups/umcg-atd/tmp06/pseudo/${mappingDestinationFolder}"
rm -f "${mappingfile%.*}.final.txt"


while read line
do
	project=$(echo "$line" | awk '{print $5}')
	sampleID=$(echo "$line" | awk '{print $2}')
	newID=$(echo "$line" | awk '{print $3}')
	echo -e "${project}\t${sampleID}\t${newID}"
	rsync -v "/groups/umcg-gd/prm0"*"/projects/${project}"*"/run01/results/variants/"*"${sampleID}"*".vcf.gz" "cf-porch+copperfist:/groups/umcg-atd/tmp06/pseudo/${mappingDestinationFolder}"

	echo -e "${sampleID}\t${newID}" >>  "${mappingfile%.*}.final.txt"

done<"${mappingfile%.*}.merged.txt"

