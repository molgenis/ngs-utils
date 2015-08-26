#!/usr/local/bin/bash

##
# MergeDatasets.sh <dataset_1> <dataset_2> <dataset_N> > <output_file>/bgzip <outputfile>
##
# B Miedema
##
# This scripts merges two or more vcf files into one vcf file and strips away all info data
#
##
#  
# 

vcfFiles=($@)

if [ ${#vcfFiles[@]} -eq 0 ]; then 
	echo "This scripts merges two or more vcf files into one vcf file and strips away all info data"
	echo "usage: This scripts merges two or more vcf files into one vcf file and strips away all info data"
fi

lines=()
header=()


Exists (){
  if [ "$2" != in ]; then
    echo "Incorrect usage."
    echo "Correct usage: exists {key} in {array}"
    return
  fi   
  eval '[ ${'$3'[$1]+muahaha} ]'  
}

declare -A variantLines
headerPresent="false"

for item in ${vcfFiles[@]}
do
	while read ln; do 
		if [[ ${ln} == "#"* ]];
			then
			if [ ${headerPresent} == "false" ]; 
			then			
				if [[ ${ln} != "##"*"INFO"* ]];
				then
					echo ${ln[@]}
				fi
			fi
		else
			IFS=$'\t' read -r -a VcfFileLineCollumns <<< "$ln"
			
			IFS=$',' read -r -a altAlleles <<< "${VcfFileLineCollumns[4]}"
			
			for altAllele in ${altAlleles[@]};do
				
			
				
				key="${VcfFileLineCollumns[0]}:${VcfFileLineCollumns[1]}:${VcfFileLineCollumns[3]}:${altAllele}"
				
			
				if Exists key in variantLines;
				then
					variantLines[${key}]=variantLines[${key}]+1
				else
					variantLines[${key}]="${VcfFileLineCollumns[0]}\t${VcfFileLineCollumns[1]}\t${VcfFileLineCollumns[2]}\t${VcfFileLineCollumns[3]}\t${altAllele}\t${VcfFileLineCollumns[5]}\t${VcfFileLineCollumns[6]}\t."
				fi
			done	
		fi
			
	done < <(gzcat ${item})
	
	headerPresent="true"
done

for line in ${variantLines[@]};
do
	echo -e $line
done









