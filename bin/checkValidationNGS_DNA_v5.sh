#!/bin/bash

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
	-i   inputFile
	-t   inputType (vcf or vcf.gz) (default= vcf.gz)
	-t   workDir working directory (default= CURRENTDIR)
	-b   which build, ucsc (chr1, chr2,chrX etc) or regular (1,2,X etc) default is regular
	-o   outputFolder (default:\${workDir}/output/)
	-v   validationFolder, folder where the vcfs are with the SNPs that should be found back (default=/groups/umcg-gd/prm06/projects/validationVcfs/)
	-l   validationLevel (all|1|2|3|4) default is 4
					1	is old validation (finding back some SNPs in 11 samples)
					2	vkgl standard = giab sample vs hc callset and vice versa
					3	frankenstein
					4	is running both option 1 and 2
					all	is running options 1,2 and 3 
===============================================================================================================
EOH
	trap - EXIT
	exit 0
}

function doVariantEval(){

	folder="${validationFolderTmp}"
	
	mapfile -t validationFiles < <(find "${folder}" -maxdepth 1 -name "*.${inputType}")
	if [[ "${#validationFiles[@]:-0}" -eq '0' ]]
	then
		echo "There are no files found in: ${folder}"
	else
		echo "##OPEN## COMPARING vcf files, all variants should be found back!##" >>  "${outputFolder}/output.txt"

		for i in "${validationFiles[@]}"
		do
			name=$(basename "${i}" ".${inputType}")
			inputFile=$(ls "${inputFolder}/"*"${name}"*".${inputType}")
			bgzippedInput=${inputFile%.*}.bgz

			zcat "${inputFile}" | bgzip -c > "${bgzippedInput}" 
			tabix -p vcf "${bgzippedInput}" 
			inputFile="${bgzippedInput}"

			java -jar "${EBROOTGATK}/GenomeAnalysisTK.jar" \
			-T VariantEval \
			-R '/apps/data/1000G/phase1/human_g1k_v37_phiX.fasta' \
			-o "${outputFolder}/output.${name}.eval.grp" \
			--eval "${inputFile}" \
			--comp "${i}"

		
			check=$(awk '{if (NR==5){if ($11 == "100.00"){print "correct"}}}' "${outputFolder}/output.${name}.eval.grp")
			if [ "${check}" == "correct" ]
			then
				zcat "${i}" | awk -v sample="${name}" 'BEGIN {OFS="  "}{if ($1 !~ /^#/){print sample,$1,$2,$4,$5,"FOUND BACK"}}' >> "${outputFolder}/output.txt"
			else
				zcat "${i}" | awk -v sample="${name}" 'BEGIN {OFS="  "}{if ($1 !~ /^#/){print sample,$1,$2,$4,$5,"Not 100% concordant!"}}' >> "${outputFolder}/output.txt"
			
			fi
		done
		echo "##CLOSE## COMPARING vcf files, all variants should be found back!##" >>  "${outputFolder}/output.txt"

	fi

}

function doComparisonFiltered (){
	refCall="${1}"
	if [ "${refCall}" == "referenceCall" ]
	then
		folder="${validationFolderTmp}/filtered/referenceCall"
		output="${outputFolder}/filtered/referenceCall"
	else
		folder="${validationFolderTmp}/filtered/"
		output="${outputFolder}/filtered/"
	fi
	mkdir -p "${output}"

	mapfile -t inputFiles < <(find "${folder}" -maxdepth 1 -name "*.${inputType}")
	if [[ "${#inputFiles[@]:-0}" -eq '0' ]]
	then
		echo "There are no files found in: ${folder}"
	else
		for i in "${inputFiles[@]}"
		do
			name=$(basename "${i}" ".${inputType}")

			if [[ "${inputType}" == "vcf.gz" ]]
			then
				validationSample=$(zcat ${i} | grep -v '^#' | awk '{print $1"-"$2"-"$4"-"$5"-"$7}')
				inputSample=$(zcat "${inputFolder}/"*"${name}"*".${inputType}" | grep 18598089 | awk '{print $1"-"$2"-"$4"-"$5"-"$7}')

				if [ "${validationSample}" == "${inputSample}" ]
				then
					if [ "${refCall}" == "referenceCall" ]
					then
						zcat "${i}" | awk -v sample="${name}" 'BEGIN {OFS="  "}{if ($1 !~ /^#/){print sample,$1,$2,$4,$5,$7,"FOUND BACK,REF CALL"}}' >> "${outputFolder}/output.txt"
					else
						zcat "${i}" | awk -v sample="${name}" 'BEGIN {OFS="  "}{if ($1 !~ /^#/){print sample,$1,$2,$4,$5,$7,"FOUND BACK"}}' >> "${outputFolder}/output.txt"
					fi
				else
					echo -e "\n##OPEN## ${refCall}: Expected variant####" >>  "${outputFolder}/output.txt"
					zcat "${i}" | awk -v sample="${name}" 'BEGIN {OFS="  "}{if ($1 !~ /^#/){print sample,$1,$2,$4,$5,$7,"Not 100% concordant!"}}' >> "${outputFolder}/output.txt"
					zcat "${inputFolder}/"*"${name}"*".${inputType}" | grep '18598089' >> "${outputFolder}/output.txt"
					zcat "${i}" >> "${outputFolder}/output.txt"
					echo -e "\n##CLOSE## ${refCall}: Expected variant####" >>  "${outputFolder}/output.txt"

				fi
			fi

		done
	fi
}

function checkFrankenstein() {

	knownVariants="${validationFolderTmp}/ValidationSet.annotated.vcf.gz"
	knownMissingVariants="${validationFolderTmp}/knownMissingVariants.txt"

	totalNoOfKnownVariants=$(zcat "${knownVariants}" | grep -v '^#' | wc -l)

	echo "Testing Frankenstein on ${totalNoOfKnownVariants} known variants"

	tmpFolder="${outputFolder}/tmp"
	input="${inputFolder}/Frankenstein_1.final.vcf.gz"

	echo "comparing \"true\" set: /groups/umcg-atd/tmp01/ValidationSet.annotated.vcf.gz"
	echo "with new: ${input}"
	echo "working folder: ${tmpFolder}"

	bedtools intersect -header -b "${input}" -a "${knownVariants}" > "${tmpFolder}/SameAsFrank.vcf"
	bgzip -c "${tmpFolder}/SameAsFrank.vcf" > "${tmpFolder}/SameAsFrank.vcf.gz"

	totalNoOfKnownVariants=$(zcat "${knownVariants}" | grep -v '^#' | wc -l)
	totalSame=$(cat "${tmpFolder}/SameAsFrank.vcf" | grep -v '^#' | wc -l)

	echo "found ${totalSame} of in total ${totalNoOfKnownVariants} back in the ${input}"

	tortilla.sh -v -1  "${tmpFolder}/SameAsFrank.vcf.gz" -2 "${knownVariants}" -o "${outputFolder}"


	if diff -q "${knownMissingVariants}" "${outputFolder}/notInVcf1.txt"
	then
		echo "Everything is OK, the missing variants are the same as the already missing variants"
	else
		echo "diff ${knownMissingVariants} ${outputFolder}/notInVcf1.txt"
		echo "OOPS, there are differences between the known missing variants and the sample!"
		diff -y "${knownMissingVariants}" "${outputFolder}/notInVcf1.txt"
	fi
}

function giabvshc(){
	giabSample="${1}"
	path=$(dirname "${giabSample}")	
	sampleName=$(basename "${giabSample}")
	sampleName=${sampleName%%.*}

	if [ "${buildType}" == "regular" ]
	then
		namingConvention=""
		refGenome="/apps/data/1000G/phase1/human_g1k_v37_phiX.fasta"
		bedfile="/apps/data/Agilent/Exoom_v3/human_g1k_v37/captured.merged.bed"
	elif [ "${buildType}" == "ucsc" ]
	then
		namingConvention="ucsc"
		refGenome="/apps/data/GRC/GCA_000001405.14_GRCh37.p13_full_analysis_set.fna"
		bedfile="/apps/data/Agilent/Exoom_v3-ucsc/human_g1k_v37/captured.merged.bed"
	else
		echo "buildType ${buildType} unknown"
	fi
		

	##capture no union bedfile
	bedtools intersect \
		-header \
		-a "${giabSample}" \
		-b "/apps/data/NIST/${namingConvention}/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.19_2mindatasets_5minYesNoRatio_noMT.bed" \
		> "${outputFolder}/tmp/${sampleName}.union.vcf"

	#split giabFile into INDEL and SNP	
	java -XX:ParallelGCThreads=1 -Xmx5g -jar "${EBROOTGATK}/GenomeAnalysisTK.jar" \
	-R "${refGenome}" \
	-T SelectVariants \
	--variant "${outputFolder}/tmp/${sampleName}.union.vcf" \
	-o "${outputFolder}/tmp/${sampleName}_INDEL.union.vcf" \
	--selectTypeToInclude INDEL \
	-sn "${sampleName}"

	echo "${sampleName}_INDEL.vcf done"

	#Select SNPs and MNPs
	java -XX:ParallelGCThreads=1 -Xmx5g -jar "${EBROOTGATK}/GenomeAnalysisTK.jar" \
	-R "${refGenome}" \
	-T SelectVariants \
	--variant "${outputFolder}/tmp/${sampleName}.union.vcf" \
	-o  "${outputFolder}/tmp/${sampleName}_SNP.union.vcf" \
	--selectTypeToExclude INDEL \
	-sn "${sampleName}"

	echo "${sampleName}_SNP.vcf done"
	
	## capture 20bp bedfile SNP
	bedtools intersect \
		-header \
		-a "${outputFolder}/tmp/${sampleName}_SNP.union.vcf" \
		-b "${bedfile}" \
		> "${outputFolder}/tmp/${sampleName}_SNP.union.20bp.vcf"
	
	## capture 20bp bedfile INDEL
	bedtools intersect \
		-header \
		-a "${outputFolder}/tmp/${sampleName}_INDEL.union.vcf" \
		-b "${bedfile}" \
		> "${outputFolder}/tmp/${sampleName}_INDEL.union.20bp.vcf"

	## Do comparison per type
	firstLine="false"
	for type in SNP INDEL
	do
		grep "^#" "${outputFolder}/tmp/${sampleName}_${type}.union.20bp.vcf" > "${outputFolder}/tmp/${sampleName}_${type}.union.20bp_filtered.vcf"
		grep -v "^#" "${outputFolder}/tmp/${sampleName}_${type}.union.20bp.vcf" |  grep -v "0/0" | grep -v "\./\." | grep "PASS" | grep -v "1/0" >> "${outputFolder}/tmp/${sampleName}_${type}.union.20bp_filtered.vcf"
		echo "comparing ${type}: ${sampleName}vsHC"
		## sample vs HC callset (precision)
		java -jar "${EBROOTGATK}/GenomeAnalysisTK.jar" \
		-T VariantEval \
		-R "${refGenome}" \
		-o "${outputFolder}/tmp/precision-${sampleName}vsHC_union.${type}.vcf" \
		--eval "${outputFolder}/tmp/${sampleName}_${type}.union.20bp_filtered.vcf" \
		--comp "/apps/data/NIST/${namingConvention}/GIAB_HC.${type}_20bp.vcf"

		if [[ ${firstLine} == "false" ]]
		then
			## print header too
			head -4 "${outputFolder}/tmp/precision-${sampleName}vsHC_union.${type}.vcf" | tail -1 | awk '{OFS="\t"}{print "measurement","Type",$6,$7,$8,$9,$10,$11,"FinalConcordance"}' > "${outputFolder}/precision_output.txt"
		fi
		head -5 "${outputFolder}/tmp/precision-${sampleName}vsHC_union.${type}.vcf" | tail -1 | awk -v type=${type} '{OFS="\t"}{print "precision",type,$6,$7,$8,$9,$10,$11,(($10/$6)*100)}' >> "${outputFolder}/precision_output.txt"
	
		echo "comparing ${type}: HCvs${sampleName}"
		##HC callset vs sample (sensitivity)
		java -jar "${EBROOTGATK}/GenomeAnalysisTK.jar" \
		-T VariantEval \
		-R "${refGenome}" \
		-o "${outputFolder}/tmp/sensitivity-HCvs${sampleName}_union.${type}.vcf" \
		--comp "${outputFolder}/tmp/${sampleName}_${type}.union.20bp_filtered.vcf" \
		--eval "/apps/data/NIST/${namingConvention}/GIAB_HC.${type}_20bp.vcf"

		if [[ ${firstLine} == "false" ]]
		then
			## print header too
			head -4 "${outputFolder}/tmp/sensitivity-HCvs${sampleName}_union.${type}.vcf" | tail -1 | awk '{OFS="\t"}{print "measurement","Type",$6,$7,$8,$9,$10,$11,"FinalConcordance"}' > "${outputFolder}/sensitivity_output.txt"
			firstLine="true"
		fi

		head -5 "${outputFolder}/tmp/sensitivity-HCvs${sampleName}_union.${type}.vcf" | tail -1 | awk -v type=${type} '{OFS="\t"}{print "sensitivity",type,$6,$7,$8,$9,$10,$11,(($10/$6)*100)}' >> "${outputFolder}/sensitivity_output.txt"
	done
}

while getopts "i:o:w:v:t:b:l:h" opt; 
do
	case $opt in h)showHelp;; i)inputFile="${OPTARG}";; w)workDir="${OPTARG}";; o)outputFolder="${OPTARG}";; v)validationFolderPrm="${OPTARG}";; b)buildType="${OPTARG}";; t)inputType="${OPTARG}";; l)validationLevel="${OPTARG}";;
esac 
done

if [[ -z "${inputFile:-}" ]]; then showHelp ; echo "inputFile is not specified" ; fi ; echo "inputFile=${inputFile}"
if [[ -z "${workDir:-}" ]]; then workDir="$(pwd)" ; mkdir -p "${workDir}/input/" ; fi ; echo "workDir=${workDir}"
if [[ -z "${outputFolder:-}" ]]; then mkdir -p "${workDir}/output/tmp" ; outputFolder="${workDir}/output/" ; fi ; echo "outputFolder=${outputFolder}"
if [[ -z "${buildType:-}" ]]; then buildType="regular" ; fi ; echo "buildType=${buildType}" 
if [[ -z "${inputType:-}" ]]; then inputType="vcf.gz" ; fi ; echo "inputType=${inputType}"	
if [[ -z "${validationFolderPrm:-}" ]]; then validationFolderPrm="/groups/umcg-gd/prm06/projects/validationVcfs/" ; fi ; echo "validationFolderPrm=${validationFolderPrm}"
if [[ -z "${validationLevel:-}" ]]; then validationLevel="all" ; fi ; echo "validationLevel=${validationLevel}" 

ml GATK 
ml HTSlib
ml BEDTools
ml ngs-utils

rm -rf "${outputFolder}/tmp/"
mkdir -p "${outputFolder}/tmp/"

echo '' > "${outputFolder}/output.txt"

if [[ "${validationLevel}" != "all" && "${validationLevel}" != "1" && "${validationLevel}" != "2" && "${validationLevel}" != "3" && "${validationLevel}" != "4" ]] 
then
	echo "this is an unknown validationLevel [${validationLevel}]"
	echo "bye bye"
	exit 1
fi
if [[ "${validationLevel}" == "all" || "${validationLevel}" == "1" || "${validationLevel}" == "4" ]]
then

	whichHost=$(hostname -s)

	if [[ "${whichHost}" == "leucine-zipper"  || "${whichHost}" == "zinc-finger" ]]
	then
		validationFolderTmp="${workDir}/input/validationVcfs/"
		mkdir -p "${outputFolder}/filtered/"

		echo "copying validationVcfs"
		if [[ -f "${validationFolderTmp}/DNA087244.${inputType}" ]]
		then
			echo "already copied, skipped"
		else
			rsync -av chaperone:${validationFolderPrm}/ "${validationFolderTmp}/"
		fi
	else
		echo "please run on leucine-zipper or zinc-finger"
	fi

	doVariantEval
	doComparisonFiltered "findVariant"
	doComparisonFiltered "referenceCall"
fi

if [[ "${validationLevel}" == "all" || "${validationLevel}" == "2" || "${validationLevel}" == "4" ]]
then
	echo "starting giab hc callset"
	mapfile -t giabSample < <(find ${inputFile})
	if [[ "${#giabSample[@]:-0}" -eq '0' ]]
	then
		echo "${inputFile} not found" 
		echo "${inputFile} not found" >> "${outputFolder}/output.txt"
	else
		echo -e "##OPEN## GIAB vs HC and versa" >> "${outputFolder}/output.txt"
		for i in "${giabSample[@]}"
		do
			echo "processing ${i}.."
			giabvshc "${i}"
		done
	fi

#	echo "giab hc callset done, output can be found here: ${outputFolder}/output.txt" 
#	echo -e "Name type\t% Concordance\t((%variantsOverlap * %concordanceOverlap) / 100)" >> "${outputFolder}/output.txt"
#	lineNumberStart=$(grep -n '##OPEN## GIAB' "${outputFolder}/output.txt" | awk 'BEGIN{FS=":"}{print $1}')
#	lineNumberStop=$(grep -n 'Name type' "${outputFolder}/output.txt" | awk 'BEGIN{FS=":"}{print $1}')
#	echo "${lineNumberStart} && ${lineNumberStop}"
	#echo -e "Measurement\tType\tnVariantsEval\tnVariantsComp\tnDifference\t\tConcordance(in %)"
	
	##precision
	echo -e "Measurement\tType\tTP\tFP\tTP/TP+FP"
	echo -e "Measurement\tType\tTP\tFP\tTP/TP+FP" >> "${outputFolder}/output.txt"
	
	awk '{if (NR>1){print $1,$2"\t"$5"\t"$4"\t"$6}}' "${outputFolder}/precision_output.txt"
	awk '{if (NR>1){print $1,$2"\t"$5"\t"$4"\t"$6}}' "${outputFolder}/precision_output.txt" >> "${outputFolder}/output.txt"
	#awk -v start=${lineNumberStart} -v stop=${lineNumberStop} '{if (NR>(start+1) && NR<stop){if($0!=""){print $1,$2"\t"$5"\t"$4"\t"$6}}}' "${outputFolder}/precision_output.txt" >> "${outputFolder}/output.txt"

	##sensitivity
	echo -e "Measurement\tType\tTP\tFN\tTP/TP+FN"
	echo -e "Measurement\tType\tTP\tFN\tTP/TP+FN" >> "${outputFolder}/output.txt"
	
	awk '{if (NR>1){print $1,$2"\t"$5"\t"$4"\t"$6}}' "${outputFolder}/sensitivity_output.txt"
	awk '{if (NR>1){print $1,$2"\t"$5"\t"$4"\t"$6}}' "${outputFolder}/sensitivity_output.txt" >> "${outputFolder}/output.txt"

fi

if [[ "${validationLevel}" == "all" || "${validationLevel}" == "3" ]]
then
	validationFolderTmp="${workDir}/input/validationVcfs/Frankenstein/"
	checkFrankenstein
fi


