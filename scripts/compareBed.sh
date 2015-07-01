#!/bin/sh

if [ -z "$1" ]; then echo "This script compares overlap between your target panel and ./refseq.bed and ./ccds.bed. Please provide path to your panel as argument!"; fi

# define refseq and ccds bed files
refseq="refseq.bed"
ccds="ccds.bed"

space="           " # format output

# exit if they don't exist
if [ ! -f $refseq ]; then
	echo "ERROR file not found: ./$refseq"
	exit 1
fi

if [ ! -f $ccds ]; then
        echo "ERROR file not found: ./$ccds"
        exit 1
fi

mkdir -p compare

echo "${space}Remove prefix 'chr' from first column..."
	declare -a files=("$1" "$refseq" "$ccds")
	for file in ${files[@]}; do
		awk -F '\t' 'BEGIN {OFS=FS} { print $1,$2,$3,$4 }' $file | sed '/^chr/s/^...//' > compare/$file
	done

cd compare/

echo "${space}Load bedtools..."
	module load bedtools/2.22.0

function create_coverage_files {
	bed1=$1
	bed2=$2
	bed1vs2=$1_vs_$2
	bed2vs1=$2_vs_$1
	
	bigspace=${space}${space}

	echo "${bigspace}Create ${bed1vs2}..."
	bedtools coverage -b $bed1 -a $bed2 > ${bed1vs2}
	awk '{ if ($8 < 1) print $0 }' $bed1vs2 > ${bed1vs2}.notFullyCovered
	sort -V -k 1,2 ${bed1vs2}.notFullyCovered > ${bed1vs2}.notFullyCovered.sorted

	echo "${bigspace}Create ${bed2vs1}..."
	bedtools coverage -a $bed1 -b $bed2 > ${bed2vs1}

	echo "${bigspace}Get N?_* names of regions covered at least one base..."
	awk '{ if ($8 > 0) {split($4,a,"_exon_"); print a[1]} }' ${bed2vs1} | sort -u > ${bed2vs1}.regionNames
	while read line; do
		grep ${line}_ ${bed2vs1} >> ${bed2vs1}.touchedRegions
	done < ${bed2vs1}.regionNames

	awk '{ if ($8 < 1) print $0 }' ${bed2vs1}.touchedRegions > ${bed2vs1}.touchedRegions.notFullyCovered
	sort -V -k 1,2 ${bed2vs1}.touchedRegions.notFullyCovered > ${bed2vs1}.touchedRegions.notFullyCovered.sorted
}

echo "${space}REFSEQ..."
	create_coverage_files $1 $refseq

echo "${space}CCDS..."
	create_coverage_files $1 $ccds
echo
echo Done! Please find your results in ./compare/
