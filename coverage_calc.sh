#!/bin/bash

set -e
set -u

underline=`tput smul`
normal=`tput sgr0`
bold=`tput bold`

function usage () {
echo -e "
${bold}
echo 'to run: sh coverage_calc.sh \${kit} \${bam}'
echo "e.g. sh coverage_calc NEURO_v3 /path/to/bam"

${bold}Arguments${normal}

        Required:
        -k|--kit                Kit for which coverage has to be determined
        -i|--input              bam input file

        Optional:
        -w|--workingdir         default (this dir)
        -r|--reference          Which reference file is used (default: /apps/data/1000G/phase1/human_g1k_v37_phiX.fasta)
        -b|--base               calculation of coverage per base [type false or true] (default:false)
        -t|--target             calculation of coverage per target[type false or true] (default:false)"
}

#
# Now goes through all the options with a case and using shift to analyse 1 argument at a time.
# $1 identifies the first argument, and when we use shift we discard the first argument, so $2 becomes $1 and goes again through the case.
#

while getopts "k:i:w:r:b:t" opt;
do
	case $opt in k)KIT="${OPTARG}";; i)INPUT="${OPTARG}";; w)WORKDIR="${OPTARG}";; r)REF="${OPTARG}";; b)BASE=${OPTARG};; t)TARGET="${OPTARG}";;
        esac
done

# Check required options were provided.

if [[ -z "${KIT-}" ]]; then
        usage
        exit 1
fi
if [[ -z "${INPUT-}" ]]; then
        usage
        exit 1
fi
if [[ -z "${BASE-}" ]]; then
        BASE="false"
else
	BASE='true'
fi
if [[ -z "${TARGET-}" ]]; then
        TARGET="false"
fi
if [[ -z "${WORKDIR-}" ]]; then
        WORKDIR=$(pwd)
fi
if [[ -z "${REF-}" ]]; then
        REF="/apps/data/1000G/phase1/human_g1k_v37_phiX.fasta"
fi

module load ngs-utils
module load GATK/3.7-Java-1.8.0_74
module list

echo ${BASE}
coveragePerBaseDir="/apps/data/UMCG/Diagnostics/CoveragePerBase/"
coveragePerTargetDir="/apps/data/UMCG/Diagnostics/CoveragePerTarget/"

if [[ "${BASE}" == 'true' ]]
then
	### Per base bed files
	bedfile="${KIT}"
	if [ -d "${coveragePerBaseDir}/${bedfile}" ]
	then
		for i in $(ls -d "${coveragePerBaseDir}/${bedfile}"/*)
		do
			perBase=$(basename "${i}")
			perBaseDir=$(echo $(dirname "${i}")"/${perBase}/human_g1k_v37/")
			echo "perBaseDir: ${perBaseDir}"
			java -Xmx10g -XX:ParallelGCThreads=4 -jar "${EBROOTGATK}/GenomeAnalysisTK.jar" \
			-R "${REF}" \
			-T DepthOfCoverage \
			-o "${INPUT}.${bedfile}.coveragePerBase" \
			--omitLocusTable \
			-I "${INPUT}" \
			-L "${perBaseDir}/${perBase}.interval_list"

			sed '1d' "${INPUT}.${bedfile}.coveragePerBase" > "${INPUT}.${bedfile}.coveragePerBase_withoutHeader"
			sort -V "${INPUT}.${bedfile}.coveragePerBase_withoutHeader" > "${INPUT}.${bedfile}.coveragePerBase_withoutHeader.sorted"
			paste "${perBaseDir}/${perBase}.uniq.per_base.bed" "${INPUT}.${bedfile}.coveragePerBase_withoutHeader.sorted" > "${INPUT}.${bedfile}.combined_bedfile_and_samtoolsoutput.txt"

			##Paste command produces ^M character
			perl -p -i -e "s/\r//g" "${INPUT}.${bedfile}.combined_bedfile_and_samtoolsoutput.txt"

			echo -e "Index\tChr\tChr Position Start\tDescription\tMin Counts\tCDS\tContig" > "${INPUT}.${bedfile}.coveragePerBase.txt"

			awk -v OFS='\t' '{print NR,$1,$2,$4,$6,"CDS","1"}' "${INPUT}.${bedfile}.combined_bedfile_and_samtoolsoutput.txt" >> "${INPUT}.${bedfile}.coveragePerBase.txt"

			#remove phiX
			grep -v "NC_001422.1" "${INPUT}.${bedfile}.coveragePerBase.txt" > "${INPUT}.${bedfile}.coveragePerBase.txt.tmp"
			mv "${INPUT}.${bedfile}.coveragePerBase.txt.tmp" "${INPUT}.${bedfile}.coveragePerBase.txt"
			echo "phiX is removed for ${INPUT} perBase"

		done
	else
		echo "There are no CoveragePerBase calculations for this bedfile: ${bedfile}"
	fi
fi

if [[ "${TARGET}" == 'true' ]]
then

	## Per target bed files
	if [ -d "${coveragePerTargetDir}/${bedfile}" ]
	then
		for i in $(ls -d "${coveragePerTargetDir}/${bedfile}"/*)
		do
			perTarget=$(basename "${i}")
			perTargetDir=$(echo $(dirname "${i}")"/${perTarget}/human_g1k_v37/")
			java -Xmx10g -XX:ParallelGCThreads=4 -jar "${EBROOTGATK}/GenomeAnalysisTK.jar" \
			-R "${REF}" \
			-T DepthOfCoverage \
			-o "${INPUT}.${bedfile}.coveragePerTarget" \
			-I "${INPUT}" \
			--omitDepthOutputAtEachBase \
			-L "${perTargetDir}/${perTarget}.interval_list"

			awk -v OFS='\t' '{print $1,$3}' "${INPUT}.${bedfile}.coveragePerTarget.sample_interval_summary" | sed '1d' > "${INPUT}.${bedfile}.coveragePerTarget.coveragePerTarget.txt.tmp.tmp"
			sort -V "${INPUT}.${bedfile}.coveragePerTarget.coveragePerTarget.txt.tmp.tmp" > "${INPUT}.${bedfile}.coveragePerTarget.coveragePerTarget.txt.tmp"
			paste "${INPUT}.${bedfile}.coveragePerTarget.coveragePerTarget.txt.tmp" "${INPUT}.${bedfile}.genesOnly" > "${INPUT}.${bedfile}.coveragePerTarget_inclGenes.txt"

			##Paste command produces ^M character
			perl -p -i -e "s/\r//g" "${INPUT}.${bedfile}.coveragePerTarget_inclGenes.txt"

			awk 'BEGIN { OFS = "\t" } ; {split($1,a,":"); print a[1],a[2],$2,$3}' "${INPUT}.${bedfile}.coveragePerTarget_inclGenes.txt" | awk 'BEGIN { OFS = "\t" } ; {split($0,a,"-"); print a[1],a[2]}' > "${INPUT}.${bedfile}.coveragePerTarget_inclGenes_splitted.txt"

			if [ -d "${INPUT}.${bedfile}.coveragePerTarget.txt" ]
			then
				rm "${INPUT}.${bedfile}.coveragePerTarget.txt"
			fi

			echo -e "Index\tChr\tChr Position Start\tChr Position End\tAverage Counts\tDescription\tReference Length\tCDS\tContig" > "${INPUT}.${bedfile}.coveragePerTarget.txt"
			awk '{OFS="\t"} {len=$3-$2} {print NR,$0,len,"CDS","1"}' "${INPUT}.${bedfile}.coveragePerTarget_inclGenes_splitted.txt" >> "${INPUT}.${bedfile}.coveragePerTarget.txt"

			#Remove phiX
			grep -v "NC_001422.1" "${INPUT}.${bedfile}.coveragePerTarget.txt" > "${INPUT}.${bedfile}.coveragePerTarget.txt.tmp"
			mv "${INPUT}.${bedfile}.coveragePerTarget.txt.tmp" "${INPUT}.${bedfile}.coveragePerTarget.txt"
			echo "phiX is removed for ${INPUT} perTarget"
		done
	else
		echo "There are no CoveragePerTarget calculations for this bedfile: ${bedfile}"
	fi
fi
