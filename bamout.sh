#!/bin/bash

set -e 
set -u

echo 'to run: sh bamout.sh ${region} ${bam} (optional ${gender})'
echo "e.g. sh bamout 1:1000-2000 /path/to/bam"
region=$1
bam=$2

ml GATK/3.7-Java-1.8.0_74 

gender=""

if [[ $region == *"X"* || $region == *"x"* ]]
then
	if [ -z "$3" ]
  	then
    		echo "No gender supplied, please fill in a gender (Male or Female) when using chromosome X region"
	else

		gender=$3
	fi
fi

myregion=$(echo "${region}" | tr : _)
if [ 1 == 0 ]
then
java -XX:ParallelGCThreads=2 -Xmx12g -jar \
${EBROOTGATK}/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R /apps/data/UCSC/GRCh38/GRCh38_no_alt_plus_hs38d1_analysis_set/reference.fa \
-I ${bam} \
-bamout ${bam}.${myregion}.bam
-o "${bam}".vcf \
-L "${region}" \
--emitRefConfidence GVCF \
-ploidy 2
fi

echo "${region}"

java -XX:ParallelGCThreads=2 -Xmx12g -jar \
${EBROOTGATK}/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R /apps/data/UCSC/GRCh38/GRCh38_no_alt_plus_hs38d1_analysis_set/reference.fa \
-I ${bam} \
-bamout ${bam}.${myregion}.bam \
-variant_index_type LINEAR \
-variant_index_parameter 128000 \
-o "${bam}".vcf \
-L "${region}" \
--emitRefConfidence GVCF \
-ploidy 2
