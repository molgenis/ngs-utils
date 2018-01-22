set -e
set -u


function showHelp() {
	#
	# Display commandline help on STDOUT.
	#
	cat <<EOH
===============================================================================================================
Script to copy (sync) data from a succesfully finished analysis project from tmp to prm storage.
Usage:
	$(basename $0) OPTIONS
Options:
	-h   Show this help.
	-i   inputFolder
	-t   inputType (vcf or vcf.gz) (default= vcf.gz)
	-o   outputFolder (default:\${inputFolder}/output/)
	-v   validationFolder, folder where the vcfs are with the SNPs that should be found back (default=/groups/umcg-gd/prm02/projects/NGS_DNA_Verification_test_ZF/validationVcfs/)
===============================================================================================================
EOH
	trap - EXIT
	exit 0
}


while getopts "i:o:v:t:h" opt; 
do
	case $opt in h)showHelp;; i)inputFolder="${OPTARG}";; o)outputFolder="${OPTARG}";; v)validationFolder="${OPTARG}";; t)inputType="${OPTARG}";;
esac 
done

if [[ -z "${inputFolder:-}" ]]; then showHelp ; echo "inputFolder is not specified" ; fi ; echo "inputFolder=${inputFolder}"
if [[ -z "${outputFolder:-}" ]]; then mkdir -p ${inputFolder}/output/ ; outputFolder=${inputFolder}/output/ ; fi ; echo "outputFolder=${outputFolder}"
if [[ -z "${inputType:-}" ]]; then inputType="vcf.gz" ; fi ; echo "inputType=${inputType}"
if [[ -z "${validationFolder:-}" ]]; then validationFolder="/groups/umcg-gd/prm02/projects/NGS_DNA_Verification_test_ZF/validationVcfs/" ; fi ; echo "validationFolder=${validationFolder}"

ml GATK 


if [ $(hostname) != "calculon" ]
then
	mkdir -p ${inputFolder}/validationVcfs/
	echo "copying validationVcfs"
	if [ -f ${inputFolder}/validationVcfs/DNA087244.vcf ]
	then
		echo "already copied, skipped"
	else
		scp calculon.hpc.rug.nl:/groups/umcg-gd/prm02/projects/NGS_DNA_Verification_test_ZF/validationVcfs/* ${inputFolder}/validationVcfs/
	fi
	validationFolder=${inputFolder}/validationVcfs/
fi
if [ 1 == 0 ]
then
for i in $(ls ${validationFolder}/*.vcf)
do
	name=$(basename $i ".vcf")

	java -jar ${EBROOTGATK}/GenomeAnalysisTK.jar \
	-T VariantEval \
	-R /apps/data/1000G/phase1/human_g1k_v37_phiX.fasta \
	-o ${outputFolder}/output.${name}.eval.grp \
	--eval ${inputFolder}/*${name}*.${inputType} \
	--comp $i

done
fi
echo "Sample  Chr  Pos  Ref  Alt  Found?"
for i in ${validationFolder}/*.vcf
do
	name=$(basename $i ".vcf")

	referenceCall=$(grep "0/0" ${i} | wc -l)
	referenceCallMale=$(grep -P "\t0:" ${i} | wc -l)

	check=$(awk '{if (NR==5){if ($11 == "100.00"){print "correct"}}}' ${outputFolder}/output.${name}.eval.grp)

	if [ "${check}" == "correct" ]
	then
		awk -v sample=${name} 'BEGIN {OFS="  "}{if ($1 !~ /^#/){print sample,$1,$2,$4,$5,"FOUND BACK"}}' ${i}
	elif [[ ${referenceCall} -ne 0 || ${referenceCallMale} -ne 0 ]]
	then
		awk -v sample=${name} 'BEGIN {OFS="  "}{if ($1 !~ /^#/){print sample,$1,$2,$4,$5,"FOUND BACK,REF CALL"}}' ${i}
	else
		awk -v sample=${name} 'BEGIN {OFS="  "}{if ($1 !~ /^#/){print sample,$1,$2,$4,$5,"Not 100% concordant!"}}' ${i}
		exit 1
	fi
done
