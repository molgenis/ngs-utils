set -eu

function parseVCF() {
	local _vcf _type
	_vcf="${1}"
	_type="${2}"
	if [[ "${extension}" == "vcf" ]]
	then
		if [[ "${bedfile}" != "empty" ]]
		then
			echo " file will be captured on bedfile "
			bedtools intersect -header -a "${_vcf}" -b "${bedfile}" > "${_vcf}.tmp"
			mv -v "${_vcf}" "${_vcf}.original"
			mv -v "${_vcf}.tmp" "${_vcf}"
		fi
			bgzip "${_vcf}"
			tabix -p vcf "${_vcf}.gz"

	else
		echo "already gz extension"
		if [[ "${bedfile}" != "empty" ]]
		then

			bedtools intersect -header -a "${_vcf}" -b "${bedfile}" > "${_vcf}.tmp"
			mv -v "${_vcf}" "${_vcf}.original"
			bgzip "${_vcf}.tmp"
			mv "${_vcf}.tmp.gz" ${_vcf}
			tabix -p vcf "${_vcf}"
		else
			echo "just compare as is."
		fi
	fi
	
}


function showHelp() {
	#
	# Display commandline help on STDOUT.
	#
	cat <<EOH
===============================================================================================================
Script to do ConcordanceChecks for array vs array manually. (see instructions:
https://github.com/molgenis/analysis-team-documents/blob/master/sops/GH-09-ConcordanceCheckArrayVsArray.md)
Usage:
	$(basename $0) OPTIONS
Options:
	-h   Show this help.
	-i   folder with index samples
	-c   folder with the compareWith samples
	-d   dataType. DRAGEN or nonDRAGEN (default) DRAGEN outputdata contains only familyNumber+umcgnumber as the samplename in the header of the vcf
	-l   limit comparison to solely these SNPS (bedfile), original file will be moved to .original  default=no bedfile
	-w   working directory (default=current dir)
	-e   extension (vcf or gz) default=gz
===============================================================================================================
EOH
	trap - EXIT
	exit 0
}

while getopts "w:i:c:d:l:e:h" opt;
do
	case $opt in h)showHelp;; i)indexFolder="${OPTARG}";; c)compareWithFolder="${OPTARG}";; d)data="${OPTARG}";; w)workDir="${OPTARG}";; l)bedfile="${OPTARG}";; e)extension="${OPTARG}";;
	esac
done

if [[ -z "${indexFolder:-}" ]]
then
	echo -e '\nERROR: Must specify an indexFolder!\n'

	showHelp
	exit 1
fi

if [[ -z "${compareWithFolder:-}" ]]
then
	echo -e '\nERROR: Must specify an compareWithFolder!\n'

	showHelp
	exit 1
fi

if [[ -z "${data:-}" ]]
then
	data="nonDRAGEN"
fi
if [[ -z "${workDir:-}" ]]
then
	workDir="$(pwd)"
fi

if [[ -z "${bedfile:-}" ]]
then
	bedfile="empty"
fi

if [[ -z "${extension:-}" ]]
then
	extension="gz"
fi

echo "extension:${extension}"
echo "bedfile:${bedfile}"

ml HTSlib
ml BEDTools/2.30.0-GCCcore-11.3.0
ml CompareGenotypeCalls
module list

tmpDir="${workDir}/tmp"
mkdir -p "${tmpDir}"
mkdir -p "${workDir}/samplesheet/"
mkdir -p "${workDir}/output"
for compare in "${compareWithFolder}/"*".${extension}"
do
	parseVCF "${compare}" 'compare'
	
	if [[ "${extension}" == "vcf" ]]
	then
		compare="${compare}.gz"
	fi
	compareWithSampleName=$(basename "${compare}")
	if [[ "${data}" == 'DRAGEN' ]]
	then
		compareWithSampleName=$(echo "${compareWithSampleName}" | awk 'BEGIN {FS="_"}{print $1"_"$2}')
	else
		compareWithSampleName="${compareWithSampleName%%.*}"
	fi

	for index in "${indexFolder}/"*".${extension}"
	do
		parseVCF "${index}" 'index'
		if [[ "${extension}" == "vcf" ]]
		then
			index="${index}.gz"
		fi
		indexSampleName=$(basename "${index}")
		if [[ "${data}" == 'DRAGEN' ]]
		then
			indexSampleName=$(echo "${indexSampleName}" | awk 'BEGIN {FS="_"}{print $1"_"$2}')
		else
			indexSampleName="${indexSampleName%%.*}"
		fi

		## create samplesheet 
		sampleSheet="${workDir}/samplesheet/${indexSampleName}_${compareWithSampleName}.sampleId.txt"
		echo -e "data1Id\tdata2Id\tlocation1\tlocation2" > "${sampleSheet}"
		echo -e "${indexSampleName}\t${compareWithSampleName}\t${index}\t${compare}" >> "${sampleSheet}"

		java -XX:ParallelGCThreads=1 -Djava.io.tmpdir="${tmpDir}" -Xmx9g -jar "${EBROOTCOMPAREGENOTYPECALLS}/CompareGenotypeCalls.jar" \
		-d1 "${index}" \
		-D1 VCF \
		-d2 "${compare}" \
		-D2 VCF \
		-ac \
		--sampleMap "${sampleSheet}" \
		-o "${workDir}/output/${indexSampleName}_${compareWithSampleName}" \
		-sva

		echo "${indexSampleName}_${compareWithSampleName} done: output/${indexSampleName}_${compareWithSampleName}"
	done
done

