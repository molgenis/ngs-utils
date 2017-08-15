set -e
set -u

function showHelp() {
	#
	# Display commandline help on STDOUT.
	#
	cat <<EOH
===============================================================================================================
Usage:
	$(basename $0) OPTIONS
Options:
	-h	Show this help.

	Required:
	-i	inputFolder (where all the inputfiles are)
	-c	capturingKit (name of capturingKit (e.g. Agilent\/ONCO_v3)
	-p	projectName

	Optional:
	-w	workDir (default is this directory)
	-m	mateOneName (difference between mateNames, e.g. R1 for mate 1 and R2 for mate 2) (default= _R1_)
		N.B. please also specify the character before and after the 1 (e.g. _1.)

Output will be written in workDir with the name: {projectName}.csv
===============================================================================================================
EOH
	trap - EXIT
	exit 0
}

while getopts "i:m:p:c:h" opt; 
do
	case $opt in h)showHelp;; w)workDir="${OPTARG}";; i)inputFolder="${OPTARG}";; m)pairedMateNameOne="${OPTARG}";; p)projectName="${OPTARG}";; c)capturingKit="${OPTARG}";; 
	esac
done


if [[ -z "${inputFolder:-}" ]]; then
	echo -e '\nERROR: Must specify an inputFolder\n'

	showHelp
	exit 1
fi

if [[ -z "${workDir:-}" ]]; then
	workDir=$(pwd)
fi
if [[ -z "${pairedMateNameOne:-}" ]]; then
	pairedMateNameOne="_R1_"
fi
if [[ -z "${capturingKit:-}" ]]; then
	echo -e '\nERROR: Must specify a capturingKit\n'
	exit 1

fi
if [[ -z "${projectName:-}" ]]; then
	echo -e '\nERROR: Must specify a projectName\n'
	exit 1

fi
count=1
totalCount=$(ls ${inputFolder}/*${pairedMateNameOne}* | wc -l)
printf "externalSampleID,externalFastQ_1,externalFastQ_2,barcode,project,capturingKit,seqType,Gender,arrayFile,lane,sequencingStartDate,sequencer,run,flowcell\n" > "${workDir}/${projectName}.csv"
for i in $(ls ${inputFolder}/*${pairedMateNameOne}*)
do
	fileName=${i%%.*}
	baseNameFile=$(basename ${fileName})
	withoutExtension=${baseNameFile%%_*}
	printf "${withoutExtension}" >> "${workDir}/${projectName}.csv" ##externalSampleID
	printf ",${i}," >> "${workDir}/${projectName}.csv" ## externalFastQ1 
	pairedMateNameTwo=$(echo "${pairedMateNameOne}" | sed 's|1|2|')

	printf "${i}" | sed -r "s|${pairedMateNameOne}|${pairedMateNameTwo}|g"  >> "${workDir}/${projectName}.csv" ## externalFastQ1
	printf ",${withoutExtension}" >> "${workDir}/${projectName}.csv" ## barcode
	printf ",${projectName}" >> "${workDir}/${projectName}.csv" ## project
	printf ",$capturingKit" >> "${workDir}/${projectName}.csv" ## capturingKit
	printf ",PE" >> "${workDir}/${projectName}.csv" ## seqType
	printf "," >> "${workDir}/${projectName}.csv"  ## Gender
	printf "," >> "${workDir}/${projectName}.csv" ## arrayFile
	printf ",1" >> "${workDir}/${projectName}.csv" ##lane
	printf ",DummyDate" >> "${workDir}/${projectName}.csv" ##sequencingstartdate
	printf ",Dummy" >> "${workDir}/${projectName}.csv" ##sequencer
	printf ",DummyRun" >> "${workDir}/${projectName}.csv" ##run
	printf ",DummyFlowcell\n" >> "${workDir}/${projectName}.csv" ## flowcell

	echo "${count} of ${totalCount} done"

	if [ "${count}" == "${totalCount}" ]
	then
		echo "Finished"
	fi
	count=$((count+1))

done
