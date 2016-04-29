#!/usr/bin/env bash

#
# Bash sanity
#
set -e  # Exit if any subcommand or pipeline returns a non-zero status.
set -u  # Exit if any uninitialised variable is used.

usage() {
    echo '##################################################################################################'
    echo ' This tool compares 2 VCF files and shows only the different rows based on the following columns:'
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
vcf01="$(basename ${VCF1})"
vcf02="$(basename ${VCF2})"
if [ ${vcf01} == ${vcf02} ]; then
	usage
	echo
	echo "ERROR: ${vcf01} is equal to ${vcf02}."
	echo "       Make sure the filenames are unique!"
	echo
	exit 1
fi
##Remove header
grep -v '^#' $VCF1 > ${SCRATCH}/${vcf01}.removedHeader.txt
## Get only necessary columns
awk '{OFS="\t"} {print $1,$2,$4,$5,$10}' ${SCRATCH}/${vcf01}.removedHeader.txt > ${SCRATCH}/${vcf01}.stripColumns_removedHeader.txt
##get only genotype call
awk -F '[\t:]' '{OFS="-"}{print $1,$2,$3,$4,$5}' ${SCRATCH}/${vcf01}.stripColumns_removedHeader.txt > ${SCRATCH}/${vcf01}.stripped.txt


awk '{OFS="\t"} {print $1,$2,$4,$5,$10}' $VCF2 > ${SCRATCH}/${vcf02}.stripped.txt
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


sort -n -k1 ${SCRATCH}/diff.txt > ${SCRATCH}/differences.txt
perl -pi -e 's|-|\t|' ${SCRATCH}/differences.txt
printf "" > ${OUT}/vcfStats.txt

falseNegative=$(cat ${SCRATCH}/notInVcf2.txt | wc -l)
falsePositive=$(cat ${SCRATCH}/notInVcf1.txt | wc -l)
alarm=$(cat ${SCRATCH}/inconsistent.txt | wc -l)
truePos=$(cat ${SCRATCH}/truePos.txt | wc -l)

sort -V -k1 ${SCRATCH}/notInVcf2.txt > ${SCRATCH}/notInVcf2.txt.sorted
sort -V -k1 ${SCRATCH}/notInVcf1.txt > ${SCRATCH}/notInVcf1.txt.sorted
sort -V -k1 ${SCRATCH}/inconsistent.txt > ${SCRATCH}/inconsistent.txt.sorted
sort -V -k1 ${SCRATCH}/truePos.txt > ${SCRATCH}/truePos.txt.sorted

bold=`tput bold`
normal=`tput sgr0`
underline=`tput smul`
printf "chr\tpos\tref\talt\tgen\n" > ${OUT}/notInVcf2.txt
printf "chr\tpos\tref\talt\tgen\n" > ${OUT}/notInVcf1.txt
printf "${bold}$underline\t\t\t|\tvcf1\t\t|\t\tvcf2\t\t\n$normal" > ${OUT}/inconsistent.txt

printf "${bold}chr\tposition\t| ref\talt\tgen\t|\tref\talt\tgen\n${normal}" >> ${OUT}/inconsistent.txt
printf "chr\tpos\tref\talt\tgen\n" > ${OUT}/truePos.txt

cat ${SCRATCH}/notInVcf2.txt.sorted >> ${OUT}/notInVcf2.txt
cat ${SCRATCH}/notInVcf1.txt.sorted >> ${OUT}/notInVcf1.txt
cat ${SCRATCH}/inconsistent.txt.sorted >> ${OUT}/inconsistent.txt
cat ${SCRATCH}/truePos.txt.sorted >> ${OUT}/truePos.txt

perl -pi -e 's|-|\t|g' ${OUT}/notInVcf2.txt
perl -pi -e 's|-|\t|g' ${OUT}/notInVcf1.txt
perl -pi -e 's|-|\t|g' ${OUT}/inconsistent.txt
perl -pi -e 's|-|\t|g' ${OUT}/truePos.txt

total=$((truePos + alarm + falsePositive + falseNegative))

tpr=$(awk "BEGIN {printf \"%.2f\n\", ((${truePos}/$total)*100)}")
fpr=$(awk "BEGIN {printf \"%.2f\n\", ((${falsePositive}/$total)*100)}")
fnr=$(awk "BEGIN {printf \"%.2f\n\", ((${falseNegative}/$total)*100)}")
alarmRate=$(awk "BEGIN {printf \"%.2f\n\", ((${alarm}/$total)*100)}")

printf "TotalNumber:${total}\n" >> ${OUT}/vcfStats.txt
printf "TP:$truePos, TP rate: ${tpr}%%\n" >> ${OUT}/vcfStats.txt
printf "FP:$falsePositive, FP rate: ${fpr}%%\n" >> ${OUT}/vcfStats.txt
printf "FN:$falseNegative, FN rate: ${fnr}%%\n" >> ${OUT}/vcfStats.txt
printf "Inconsistent:$alarm, InconsistentRate: ${alarmRate}%%\n" >> ${OUT}/vcfStats.txt

exit 0
