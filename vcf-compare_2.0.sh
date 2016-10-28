#!/usr/bin/env bash

#
# Bash sanity
#
set -e  # Exit if any subcommand or pipeline returns a non-zero status.
set -u  # Exit if any uninitialised variable is used.

usage() {
    echo '##################################################################################################'
    echo ' This tool compares 2 VCF files (or GZIPPED vcf) and shows only the different rows based on the following columns:'
    echo '   CHROM, POS, REF, ALT and SampleID:GT'
    echo
    echo ' Usage:'
    echo '   bash vcf-compare.sh -1 <input_vcf_file> -2 <input_vcf_file>  -o <output_folder>'
    echo
    echo ' Example:'
    echo '   bash vcf-compare.sh -1 old_analysis.snps.final.vcf \'
    echo '                       -2 new_analysis.snps.final.vcf \'
    echo '                       -o ./results/'
    echo '##################################################################################################'
}

declare VCF1=""
declare VCF2=""
declare OUT=""
#
# Take the parameters given on the commandline and assign them to variables.
#
while getopts ":1:2:d:o:h" option; do
    case "${option}" in
        1)  VCF1=${OPTARG};;
        2)  VCF2=${OPTARG};;
        o)   OUT=${OPTARG};;
        h)
            usage
            exit 0
            ;;
        \?)
            echo "ERROR: Invalid option -${OPTARG}. Try \"$(basename $0) -h\" for help."
            exit 1
            ;;
        :)
            echo "ERROR: Option -${OPTARG} requires an argument. Try \"$(basename $0) -h\" for help."
            ;;
    esac
done

#
# Check if all parameters are set.
#
if [[ ${VCF1} && ${VCF2} && ${OUT} ]]; then
    echo
    echo "Comparing ${VCF1} to ${VCF2}..."
    echo
else
    usage
    echo
    echo "ERROR: missing required argument. Try \"$(basename $0) -h\" for help."
    echo
    exit 1
fi

#
# Create tmp folders.
#
SCRATCH="${OUT}/TMP/"
rm -Rf   ${SCRATCH}
mkdir -p ${SCRATCH}

printf "vcf1:${VCF1}\nvcf2:${VCF2}"> "${OUT}/runparameters.txt"
#
# select filenames from given path and check if they are unique.
#

checkformatVCF1=${VCF1##*.}
checkformatVCF2=${VCF2##*.}

if [ $checkformatVCF1 != $checkformatVCF2 ]
then
	echo "formats are not the same ($checkformatVCF1 vs $checkformatVCF2)"
	exit 1
fi

if [ "${checkformatVCF1}" == "gz" ]
then
	inputVcf1=${VCF1%.*}
	gzip -c -d $VCF1 > ${inputVcf1}

	inputVcf2=${VCF2%.*}
	gzip -c -d $VCF2 > ${inputVcf2}
elif [ "${checkformatVCF1}" == "vcf" ]
	inputVcf1=$VCF1
	inputVcf2=$VCF2
else
	echo "not a correct format, please use vcf or vcf.gz"
	exit 1
fi

vcf01="$(basename ${inputVcf1})"
vcf02="$(basename ${inputVcf2})"
if [ ${vcf01} == ${vcf02} ]; then
	usage
	echo
	echo "ERROR: ${vcf01} is equal to ${vcf02}."
	echo "       Make sure the filenames are unique!"
	echo
	exit 1
fi

##Remove header
grep -v '^#' ${inputVcf1} > ${SCRATCH}/${vcf01}.removedHeader.txt
## Get only necessary columns
awk '{OFS="\t"} {print $1,$2,$4,$5,$10}' ${SCRATCH}/${vcf01}.removedHeader.txt > ${SCRATCH}/${vcf01}.stripColumns_removedHeader.txt
##get only genotype call
awk -F '[\t:]' '{OFS="-"}{print $1,$2,$3,$4,$5}' ${SCRATCH}/${vcf01}.stripColumns_removedHeader.txt > ${SCRATCH}/${vcf01}.stripped.txt


awk '{OFS="\t"} {print $1,$2,$4,$5,$10}' ${inputVcf2} > ${SCRATCH}/${vcf02}.stripped.txt
grep -v '^#' ${SCRATCH}/${vcf02}.stripped.txt > ${SCRATCH}/${vcf02}.stripColumns_removedHeader.txt
awk -F '[\t:]' '{OFS="-"}{print $1,$2,$3,$4,$5}' ${SCRATCH}/${vcf02}.stripColumns_removedHeader.txt > ${SCRATCH}/${vcf02}.stripped.txt

declare -A arrVCF1
declare -A arrVCF2


while read line 
do
	OLDIFS=$IFS
        IFS="-"
        set $line
        chr=$1
        pos=$2
        ref=$3
        alt=$4
        gen=$5
        IFS=$OLDIFS
	myvalue="${ref}-${alt}-${gen}"
	mykey="${chr}-${pos}"
	arrVCF1["${mykey}"]="${myvalue}"

done<${SCRATCH}/${vcf01}.stripped.txt

while read line 
do
	OLDIFS=$IFS
        IFS="-"
        set $line
        chr=$1
        pos=$2
        ref=$3
        alt=$4
        gen=$5
        IFS=$OLDIFS
        arrVCF2["${chr}-${pos}"]="${ref}-${alt}-${gen}"

done<${SCRATCH}/${vcf02}.stripped.txt
printf "" > ${SCRATCH}/differences.txt
printf "" > ${SCRATCH}/diff.txt

for i in "${!arrVCF1[@]}"
do
	if [ "${arrVCF2[$i]+abc}" ] && printf "$i ${arrVCF1[$i]}\n" >> ${SCRATCH}/truePos.txt
	then
		if [ ${arrVCF1[$i]} != ${arrVCF2[$i]} ]
		then
			printf "$i \t| ${arrVCF1[$i]}\t|\t${arrVCF2[$i]} \n" >> ${SCRATCH}/diff.txt
			printf "$i \t| ${arrVCF1[$i]}\t|\t${arrVCF2[$i]} \n" >> ${SCRATCH}/inconsistent.txt
		fi
	else
		printf "$i\t${arrVCF1[$i]}\n" >> ${SCRATCH}/diff.txt
		printf "$i\t${arrVCF1[$i]}\n" >> ${SCRATCH}/notInVcf2.txt
	fi
done
for i in "${!arrVCF2[@]}"
do
        if [ "${arrVCF1[$i]+abc}" ]
        then
                if [ ${arrVCF1[$i]} != ${arrVCF2[$i]} ] 
                then
                        echo "already done"
                fi
        else
                printf "$i\t${arrVCF2[$i]}\n" >> ${SCRATCH}/diff.txt
                printf "$i\t${arrVCF2[$i]}\n" >> ${SCRATCH}/notInVcf1.txt
        fi
done

bold=`tput bold`
normal=`tput sgr0`
underline=`tput smul`

sort -n -k1 ${SCRATCH}/diff.txt > ${SCRATCH}/differences.txt
perl -pi -e 's|-|\t|' ${SCRATCH}/differences.txt
printf "" > ${OUT}/vcfStats.txt

alarm=0
falseNegative=0
falsePositive=0

##NOT IN VCF2
if [ -f ${SCRATCH}/notInVcf2.txt ]
then
	falseNegative=$(cat ${SCRATCH}/notInVcf2.txt | wc -l)
	sort -V -k1 ${SCRATCH}/notInVcf2.txt > ${SCRATCH}/notInVcf2.txt.sorted
	printf "chr\tpos\tref\talt\tgen\n" > ${OUT}/notInVcf2.txt
	cat ${SCRATCH}/notInVcf2.txt.sorted >> ${OUT}/notInVcf2.txt
	perl -pi -e 's|-|\t|g' ${OUT}/notInVcf2.txt
fi

##NOT IN VCF1
if [ -f ${SCRATCH}/notInVcf1.txt ]
then
	falsePositive=$(cat ${SCRATCH}/notInVcf1.txt | wc -l)
	sort -V -k1 ${SCRATCH}/notInVcf1.txt > ${SCRATCH}/notInVcf1.txt.sorted
	printf "chr\tpos\tref\talt\tgen\n" > ${OUT}/notInVcf1.txt
	cat ${SCRATCH}/notInVcf1.txt.sorted >> ${OUT}/notInVcf1.txt
	perl -pi -e 's|-|\t|g' ${OUT}/notInVcf1.txt
fi

## INCONSISTENT
if [ -f ${SCRATCH}/inconsistent.txt ]
then
	alarm=$(cat ${SCRATCH}/inconsistent.txt | wc -l)
	sort -V -k1 ${SCRATCH}/inconsistent.txt > ${SCRATCH}/inconsistent.txt.sorted
	printf "${bold}$underline\t\t\t|\tvcf1\t\t|\t\tvcf2\t\t\n$normal" > ${OUT}/inconsistent.txt
	printf "${bold}chr\tposition\t| ref\talt\tgen\t|\tref\talt\tgen\n${normal}" >> ${OUT}/inconsistent.txt
	cat ${SCRATCH}/inconsistent.txt.sorted >> ${OUT}/inconsistent.txt
	perl -pi -e 's|-|\t|g' ${OUT}/inconsistent.txt
	
fi

##TRUE POS
truePos=$(cat ${SCRATCH}/truePos.txt | wc -l)
sort -V -k1 ${SCRATCH}/truePos.txt > ${SCRATCH}/truePos.txt.sorted
printf "chr\tpos\tref\talt\tgen\n" > ${OUT}/truePos.txt
cat ${SCRATCH}/truePos.txt.sorted >> ${OUT}/truePos.txt

perl -pi -e 's|-|\t|g' ${OUT}/truePos.txt


total=$((truePos + alarm + falsePositive + falseNegative))
printf "TotalNumber:${total}\n" >> ${OUT}/vcfStats.txt
##TRUE POS
tpr=$(awk "BEGIN {printf \"%.2f\n\", ((${truePos}/$total)*100)}")
printf "TP:$truePos, TP rate: ${tpr}%%\n" >> ${OUT}/vcfStats.txt

if [ ${falseNegative} -ne 0 ]
then
	fnr=$(awk "BEGIN {printf \"%.2f\n\", ((${falseNegative}/$total)*100)}")
	printf "FN:$falseNegative, FN rate: ${fnr}%%\n" >> ${OUT}/vcfStats.txt	
else
	printf "FP:$falseNegative, FN rate: 0%%\n" >> ${OUT}/vcfStats.txt

fi
if [ ${falsePositive} -ne 0 ]
then	
	fpr=$(awk "BEGIN {printf \"%.2f\n\", ((${falsePositive}/$total)*100)}")
	printf "FP:$falsePositive, FP rate: ${fpr}%%\n" >> ${OUT}/vcfStats.txt
else
	printf "FP:$falsePositive, FP rate: 0%%\n" >> ${OUT}/vcfStats.txt
fi
if [ ${alarm} -ne 0 ]
then
	alarmRate=$(awk "BEGIN {printf \"%.2f\n\", ((${alarm}/$total)*100)}")
	printf "Inconsistent:$alarm, InconsistentRate: ${alarmRate}%%\n" >> ${OUT}/vcfStats.txt
else
	printf "Inconsistent:$alarm, InconsistentRate: 0%%\n" >> ${OUT}/vcfStats.txt
fi

printf "done..\nComparison can be found: $OUT \n"


exit 0
