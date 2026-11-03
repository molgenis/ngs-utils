set -eu

function showHelp() {
	#
	# Display commandline help on STDOUT.
	#
	cat <<EOH
===============================================================================================================
Usage:
	$(basename $0) -n NEWSAMPLENAME -o OLDSAMPLENAME
	Step 1) make a copy of the folder containing the files (except for huge bam file, this file can be copied to prm directly (and rename it while copying)).
	N.B. Bam file will be used by creating concordancecheck vcf file, should be in original file location .e.g /groups/umcg-genomescan/\${thisTmpDir}/\${gsBatch}/Analysis/\${oldsample}/\${oldsample}.bam
	Step 2)
	This is script should be in the folder containing all the files. This directory is also the workdir (tmp and output dir will be created automatically).
	Step 3) run script; All files will be renamed to the new name and all important files are processed and copied to the output dir)

Options:
	-h	Show this help.

	Required:
	-n	new sample name (e.g. 20212312_0000000_DNA111111_111222_GS000G_WGS_123456_678999
	-o	old sample name (100000-123-001-GS_000G-WGS_v1-678999)

	Optional:
	-b 	alternative location for the bam file (full path, including bamfile itself)
===============================================================================================================
EOH
	trap - EXIT
	exit 0
}

while getopts "n:o:b:h" opt; 
do
	case $opt in h)showHelp;; n)newsample="${OPTARG}";; o)oldsample="${OPTARG}";; b)bamFileLocation="${OPTARG}";;
	esac
done

if [[ -z "${newsample:-}" ]]
then
	echo 'Must specify newsample with -n.'
	exit 1
fi

if [[ -z "${oldsample:-}" ]]
then
	echo 'Must specify oldsample with -o.'
	exit 1
fi

if [[ "$(hostname)" == "copperfist" ]]
then
	thisTmpDir='tmp16'
elif [[ "$(hostname)" == "betabarrel" ]]
then
        thisTmpDir='tmp15'
elif [[ "$(hostname)" == "wingedhelix" ]]
then
        thisTmpDir='tmp07'
else
	echo "unkown host"
	exit 1
fi

echo "NEWSAMPLE:${newsample}"
echo "OLDSAMPLE:${oldsample}"

# extract gsBatch out of old sample name
gsBatch=$(echo "${oldsample}" | awk 'BEGIN {FS="-"}{print $1"-"$2}')
echo "GSBATCH= ${gsBatch}"

if [[ -z "${bamFileLocation:-}" ]]
then
	bamFileLocation="/groups/umcg-genomescan/${thisTmpDir}/${gsBatch}/Analysis/${oldsample}/${oldsample}.bam"
fi

bamFileName=$(basename "${bamFileLocation}")
bamLocation=$(dirname "${bamFileLocation}")

echo "BAMFILENAME: ${bamFileName}"
echo "BAMLOCATION: ${bamLocation}"

# check if necessary files are here
mapfile -t findGVCF < <(find . -maxdepth 1 -type f -name "*.gvcf.gz")
mapfile -t findCNV < <(find . -maxdepth 1 -type f -name "*.cnv.vcf.gz")
mapfile -t findBam < <(find "${bamLocation}" -name "${bamFileName}")
mapfile -t findROH < <(find . -maxdepth 1 -type f -name "*.roh.bed")
mapfile -t findBafBW < <(find . -maxdepth 1 -type f -name "*.baf.bw")
mapfile -t findTnBW < <(find . -maxdepth 1 -type f -name "*.tn.bw")

if [[ "${#findGVCF[@]}" -eq '0' || "${#findCNV[@]}" -eq '0' || "${#findBam[@]}" -eq '0' || "${#findROH[@]}" -eq '0' || "${#findBafBW[@]}" -eq '0' || "${#findTnBW[@]}" -eq '0' ]]
then
	echo "There are files missing from the list below, cannot continue:"
	echo "gVCF: ${#findGVCF[@]} (should be in this folder)"
	echo "CNV: ${#findCNV[@]} (should be in this folder)"
	echo "Bam: ${#findBam[@]} should be /groups/umcg-genomescan/${thisTmpDir}/${gsBatch}/Analysis/${oldsample}/"
	echo "ROH: ${#findROH[@]} (should be in this folder)"
	echo "BafBW: ${#findBafBW[@]} (should be in this folder)"
	echo "tnBW: ${#findTnBW[@]} (should be in this folder)"
exit 1
else
        echo "found gVCF: ${#findGVCF[@]}"
        echo "found CNV: ${#findCNV[@]}"
        echo "found Bam: ${#findBam[@]}"
        echo "found ROH: ${#findROH[@]}"
        echo "found BafBW: ${#findBafBW[@]}"
        echo "found tnBW: ${#findTnBW[@]}"
fi
exit 1
# make directories and output structure for easier copying to prm
mkdir -p -m 2775 'tmp' 'output/'


echo "${oldsample} ${newsample}" > 'newHeader.txt'
ml BCFtools


mapfile -t findsamplesWithOldNames < <(find . -type f -name "${oldsample}*")

if [[ "${#findsamplesWithOldNames[@]}" -eq '0' ]]
then
	echo "Files are already renamed, skipping renaming"
else
	# renaming all the files
	rename "${oldsample}" "${newsample}" "${oldsample}"*
fi

echo "start reheadering gVCF"
#rehader gvcf
bcftools reheader -s "newHeader.txt" "${newsample}.hard-filtered.gvcf.gz" -o "tmp/${newsample}.hard-filtered.g.vcf.gz"
tabix -p vcf "tmp/${newsample}.hard-filtered.g.vcf.gz"

echo "start reheadering CNV"
#reheader cnv
bcftools reheader -s "newHeader.txt" "${newsample}.cnv.vcf.gz" -o "tmp/${newsample}.cnv.reheadered.vcf.gz"
tabix -p vcf "tmp/${newsample}.cnv.reheadered.vcf.gz"

echo "start making concordancec check vcf"
# make concordanceCheck vcf

fasta="/apps/data//GRC/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa"



bcftools mpileup \
-Ou -f "${fasta}" "${bamFileLocation}" -R "/apps/data//UMCG/concordanceCheckSnps_GRCh38.bed" \
| bcftools call -m -Ob -o "tmp/${newsample}.concordanceCheckCalls.tmp.vcf"

bcftools reheader -s "newHeader.txt" "tmp/${newsample}.concordanceCheckCalls.tmp.vcf" \
> "tmp/${newsample}.concordanceCheckCalls.vcf"

echo "start calculating coverage per target"

#calculate coverage per Target
ml HTSlib

sample="${newsample}"

bedfile="/apps/data/UMCG/ncbiRefSeq_hg38_2022-10-28_exons_slop_50bp_MANE_COMPLETED.merged.bed"

bcftools view --regions-file "${bedfile}"  -O z -o "tmp/${sample}.hard-filtered.PseudoExome.g.vcf.gz"  "tmp/${sample}.hard-filtered.g.vcf.gz"
tabix -p vcf "tmp/${sample}.hard-filtered.PseudoExome.g.vcf.gz"

outputFile="tmp/${sample}.PseudoExome.CoverageOutput.csv"

ml purge
ml gVCF2BED/1.2.0-GCCcore-11.3.0

gvcf2bed2.py \
-I "tmp/${sample}.hard-filtered.PseudoExome.g.vcf.gz" \
-O "${outputFile}" \
-b "${bedfile}"

awk '{sumDP+=$11;sumTargetSize+=$12;sumCoverageInDpLow+=$13;sumZeroCoverage+=14}END{print "avgCov: "(sumDP/sumTargetSize)"\t%coverageBelow10: "((sumCoverageInDpLow/sumTargetSize)*100)"\t%ZeroCoverage: "((sumZeroCoverage/sumTargetSize)*100)}' "${outputFile}" > "tmp/${sample}.incl_TotalAvgCoverage_TotalPercentagebelow10x.txt"

awk 'BEGIN{OFS="\t"}{if (NR>1){print (NR-1),$1,$2+1,$3,$8,$4,$12,"CDS","1"}else{print "Index\tChr\tChr Position Start\tChr Position End\tAverage Counts\tDescription\tReference Length\tCDS\tContig"}}' "${outputFile}" > "tmp/${sample}.pseudoExome.coveragePerTarget.txt"

grep -v "NC_001422.1" "tmp/${sample}.pseudoExome.coveragePerTarget.txt" > "tmp/${sample}.pseudoExome.coveragePerTarget.txt.tmp"
mv -v "tmp/${sample}.pseudoExome.coveragePerTarget.txt.tmp" "tmp/${sample}.pseudoExome.coveragePerTarget.txt"


echo "copying data to $(pwd)/output/"

#copyToOutputDir
#copy gVCF
rsync -v "tmp/${newsample}.hard-filtered.g.vcf.gz"{,.tbi} 'output/'
rsync -v "tmp/${sample}.hard-filtered.PseudoExome.g.vcf.gz"{.tbi} 'output/'
# copy cnv
rsync -v "tmp/${newsample}.cnv.reheadered.vcf.gz" "output/${newsample}.cnv.vcf.gz"
rsync -v "tmp/${newsample}.cnv.reheadered.vcf.gz.tbi" "output/${newsample}.cnv.vcf.gz"
# copy coverage
rsync -v "tmp/${sample}.pseudoExome.coveragePerTarget.txt" 'output/'
# copy roh
rsync -v "${sample}.roh.bed" 'output/'
# copy baf.bw
rsync -v "${sample}.hard-filtered.baf.bw" 'output/'
#copy tn.bw
rsync -v "${sample}.tn.bw" 'output/'