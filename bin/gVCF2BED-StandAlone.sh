#!/bin/bash
set -eu


function showHelp() {
	#
	# Display commandline help on STDOUT.
	#
	cat <<EOH
===============================================================================================================
Script to calculate coverage based on gVCF files
Usage:
	$(basename "${0}") OPTIONS
Options:
	-h   Show this help.

    Required:
	-i   inputFolder
   
    Optional:
	-b   select which bedfile to use for the targets (default="/apps/data/Agilent/Exoom_v3/human_g1k_v37/captured.merged.bed""
	-w   workDir working directory (default= CURRENTDIR)
	-g   gatk version (3 or 4) default="3"
	-o   outputFolder (default:\${workDir}/output/)
===============================================================================================================
EOH
	trap - EXIT
	exit 0
}


while getopts "i:o:w:b:g:h" opt;
do
	# shellcheck disable=SC2249
	# shellcheck disable=SC2220
	case "${opt}" in h)showHelp;; i)inputFolder="${OPTARG}";; o)outputFolder="${OPTARG}";; w)workDir="${OPTARG}";; b)bedfile="${OPTARG}";; g)gatkVersion="${OPTARG}";;
esac
done

if [[ -z "${inputFolder:-}" ]]; then showHelp ; echo "inputFolder is not specified" ; fi ; echo "inputFolder=${inputFolder}"
if [[ -z "${bedfile:-}" ]]; then bedfile="/apps/data/Agilent/Exoom_v3/human_g1k_v37/captured.merged.bed" ; fi ; echo "bedfile=${bedfile}"
if [[ -z "${workDir:-}" ]]; then workDir="$(pwd)" ; mkdir -p "${workDir}/input/" ; fi ; echo "workDir=${workDir}"
if [[ -z "${outputFolder:-}" ]]; then mkdir -p "${workDir}/output/tmp" ; outputFolder="${workDir}/output/" ; fi ; echo "outputFolder=${outputFolder}"
if [[ -z "${gatkVersion:-}" ]]; then gatkVersion="3" ; fi ; echo "gatkVersion=${gatkVersion}"


if [ -d "${inputFolder}/" ]
then
	mapfile -t checkGzFiles < <(find "${inputFolder}/" -name '*.vcf.gz')
	if [[ "${#checkGzFiles[@]:-0}" -eq '0' ]] 
	then
		echo "{inputFolder} existing, but no .gz files detected"
		exit 1
	fi
else
	echo "stopped, ${inputFolder} not existing"
	exit 1
fi
## select all uniq samples from inputfolder
uniqSamples=$(for i in "${inputFolder}/"*"vcf.gz" ; do printf '%s\n' "${i%%.*}";done | sort -V | uniq)

printf "number of unique samples is: " 
IFS=' ' ; echo "${uniqSamples}" | wc -w
IFS=' ' read -r -a uniqSamplesArray <<< "${uniqSamples}"

if [ "${gatkVersion}" == "3" ]
then
	module load GATK/3.7-Java-1.8.0_74
elif [[ "${gatkVersion}" == "4" ]]
then
	module load GATK/4.1.4.1-Java-8-LTS
else
	echo "${gatkVersion} is not a known GATK version, please choose 3 or 4"
fi

module load HTSlib
module load gVCF2BED
module load BCFtools

for i in "${uniqSamplesArray[@]}"
do
	
	sample=$(basename "${i}")
	path=$(dirname "${i}")
	
	## GATK3
	if [ "${gatkVersion}" == "3" ]
	then
		if [ -f "${workDir}/${sample}.merged.g.vcf.gz" ]
		then
			echo "already created ${workDir}/${sample}.merged.g.vcf.gz, skipped"
		else
			xpFile=""
			INPUTS=()

			for gVCF in "${path}/${sample}"*".vcf.gz"
			do 
				if [[ "${gVCF}" == *"batch-Xp.variant.calls.g.vcf.gz"* ]]
				then
					echo "skip the Xp batch for now, this will be added later on"
					xpFile=$(ls "${path}/${sample}"*"batch-Xp.variant.calls.g.vcf.gz")
				else
					INPUTS+=("--variant ${gVCF}")
				fi
			done
			IFS=$'\n' mapfile -t sortedINPUTS < <(sort -V <<<"${INPUTS[*]}")
#			IFS=$'\n' sortedINPUTS=($(sort -V <<<"${INPUTS[*]}")) ; unset IFS

			# shellcheck disable=SC2154
			java -Xmx5g -cp "${EBROOTGATK}/GenomeAnalysisTK.jar" org.broadinstitute.gatk.tools.CatVariants \
			-R /apps/data/1000G/phase1/human_g1k_v37_phiX.fasta \
			"${sortedINPUTS[@]}" \
			-out "${workDir}/${sample}.mergedWithoutXp.g.vcf.gz"
		
			zcat "${workDir}/${sample}.mergedWithoutXp.g.vcf.gz" > "${workDir}/${sample}.mergedWithoutXp.g.vcf"
			zcat "${xpFile}" | grep -v '^#' >> "${workDir}/${sample}.mergedWithoutXp.g.vcf"

			bcftools sort "${workDir}/${sample}.mergedWithoutXp.g.vcf" -o "${workDir}/${sample}.merged.g.vcf"
			bgzip "${workDir}/${sample}.merged.g.vcf"
			tabix -p vcf "${workDir}/${sample}.merged.g.vcf.gz"
		fi
			

	## GATK4
	else
		INPUTS=()
	       
		for gVCF in "${path}/${sample}"*".vcf.gz"
		do 
			INPUTS+=("--INPUT ${gVCF}")
		done

		gatk GatherVcfs \
		"${INPUTS[@]}" \
		--OUTPUT "${workDir}/${sample}.merged.g.vcf.gz"
	fi

	outputFile="${outputFolder}/${sample}_${bedfile}.output.csv"

	echo "starting to do the calculations"

	gvcf2bed2.py \
	-I "${workDir}/${sample}.merged.g.vcf.gz" \
	-O "${outputFile}" \
	-b "${bedfile}"

	awk '{sumDP+=$11;sumTargetSize+=$12;sumCoverageInDpLow+=$13;sumZeroCoverage+=14}END{print "avgCov: "(sumDP/sumTargetSize)"\t%coverageBelow20: "((sumCoverageInDpLow/sumTargetSize)*100)"\t%ZeroCoverage: "((sumZeroCoverage/sumTargetSize)*100)}' "${outputFile}" > "${outputFile%.*}.incl_TotalAvgCoverage_TotalPercentagebelow20x.txt"

	## creeer de coveragePerTarget.txt file
	awk 'BEGIN{OFS="\t"}{if (NR>1){print (NR-1),$1,$2,$3,$8,$4,$12,"CDS","1"}else{print "Index\tChr\tChr Position Start\tChr Position End\tAverage Counts\tDescription\tReference Length\tCDS\tContig"}}' "${outputFile}" > "${outputFile%%*.}.coveragePerTarget.txt"

	echo "Raw output file is here: ${outputFile}"
	echo "final statistics can be found here: ${outputFile}.incl_TotalAvgCoverage_TotalPercentagebelow20x.txt"
	echo "coveragePerTarget file can be found here: ${outputFile%%.*}.coveragePerTarget.txt"


done



