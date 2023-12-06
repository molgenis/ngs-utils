set -eu

function showHelp() {
	#
	# Display commandline help on STDOUT.
	#
	cat <<EOH
===============================================================================================================
Script to do (many-to-many) ConcordanceChecks manually.
Usage:
	$(basename $0) OPTIONS
Options:
	-h   Show this help.
	-i   folder with index samples
	-c   folder with the compareWith samples
	-d   dataType. DRAGEN or nonDRAGEN (default) DRAGEN outputdata contains only familyNumber+umcgnumber as the samplename in the header of the vcf
	-w   working directory (default=current dir)
===============================================================================================================
EOH
	trap - EXIT
	exit 0
}

while getopts "w:i:c:d:h" opt;
do
	case $opt in h)showHelp;; i)indexFolder="${OPTARG}";; c)compareWithFolder="${OPTARG}";; d)data="${OPTARG}";; w)workDir="${OPTARG}";;
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

ml CompareGenotypeCalls


tmpDir="${workDir}/tmp"
mkdir -p "${tmpDir}"
mkdir -p "${workDir}/samplesheet/"
mkdir -p "${workDir}/output"
for compare in "${compareWithFolder}/"*".gz"
do
	compareWithSampleName=$(basename "${compare}")
	if [[ "${data}" == 'DRAGEN' ]]
	then
		compareWithSampleName=$(echo "${compareWithSampleName}" | awk 'BEGIN {FS="_"}{print $1"_"$2}')
	else
		compareWithSampleName="${compareWithSampleName%%.*}"
	fi

	for index in "${indexFolder}/"*".gz"
	do
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
