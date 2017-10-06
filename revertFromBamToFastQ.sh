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
        -i   inputBam (for 1 bam file only) 
        -f   inputFolder (location of the inputfolder containing the bam files)
        -o   outputFolder (location where to store the output, including scripts)

===============================================================================================================
EOH
	trap - EXIT
        exit 0
}


while getopts "i:o:f:h" opt;
do
	case $opt in h)showHelp;; i)inputBam="${OPTARG}";; o)outputFolder="${OPTARG}";; f)inputFolder="${OPTARG}";;
esac
done

if [[ -z "${inputBam:-}" && -z "${inputFolder:-}" ]]; then showHelp ; echo "inputBam/inputFolder is not specified" ; fi
if [[ ! -z "${inputBam:-}" && ! -z "${inputFolder:-}" ]]; then showHelp ; echo "inputBam and inputFolder cannot be both specified" ; fi

if [[ -z "${inputBam:-}" ]]; then inputBam="empty"; fi
if [[ -z "${inputFolder:-}" ]]; then inputFolder=$(dirname ${inputBam}) ; fi

if [[ -z "${outputFolder:-}" ]]; then showHelp ; echo "outputFolder is not specified"; fi

mkdir -p ${outputFolder}/scripts

echo "inputBam=${inputBam}"
echo "inputFolder=${inputFolder}"
echo "outputFolder=${outputFolder}"

arrayBams=()
##if inputBam is empty, than inputFolder is set and vice versa
if [ "${inputBam}" == "empty" ]
then
	for i in $(ls "${inputFolder}/"*.bam)
	do
		arrayBams+=(${i})
	done
else
	arrayBams+=(${inputBam})
fi

for inputBam in ${arrayBams[@]}
do
	fileN=$(basename ${inputBam})
	name=${fileN%%.*}
echo -e "#!/bin/bash
#SBATCH --job-name=ConvertToFastQ_${name}
#SBATCH --output=${outputFolder}/scripts/ConvertToFastQ_${name}.out
#SBATCH --error=${outputFolder}scripts/ConvertToFastQ_${name}.err
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task 2
#SBATCH --mem 9gb
#SBATCH --nodes 1
#SBATCH --open-mode=append

set -e
set -u

touch ${outputFolder}/scripts/ConvertToFastQ_${name}.started

module load picard

java -jar -XX:ParallelGCThreads=2 -Djava.io.tmpdir=/tmp -Xmx8g \${EBROOTPICARD}/picard.jar RevertSam \\
    I=${inputBam} \\
    O=${outputFolder}/${name}.revertsam.bam \\
    SANITIZE=true \\
    MAX_DISCARD_FRACTION=0.005 \\
    ATTRIBUTE_TO_CLEAR=XT \\
    ATTRIBUTE_TO_CLEAR=XN \\
    ATTRIBUTE_TO_CLEAR=AS \\
    ATTRIBUTE_TO_CLEAR=OC \\
    ATTRIBUTE_TO_CLEAR=OP \\
    SORT_ORDER=queryname \\
    RESTORE_ORIGINAL_QUALITIES=true \\
    REMOVE_DUPLICATE_INFORMATION=true \\
    REMOVE_ALIGNMENT_INFORMATION=true


java -jar -XX:ParallelGCThreads=2 -Djava.io.tmpdir=/tmp -Xmx8g \${EBROOTPICARD}/picard.jar SamToFastq I=${outputFolder}/${name}.revertsam.bam F=${outputFolder}/${name}_reads_1.fq F2=${outputFolder}/${name}_reads_2.fq


gzip -c ${outputFolder}/${name}_reads_1.fq > ${outputFolder}/${name}_reads_1.fq.gz
gzip -c ${outputFolder}/${name}_reads_2.fq > ${outputFolder}/${name}_reads_2.fq.gz
echo \"gzipping done, removing fq and reverted bam\"
rm -f ${outputFolder}/${name}_reads_[12].fq
rm -f ${outputFolder}/${name}.revertsam.bam
mv ${outputFolder}/scripts/ConvertToFastQ_${name}.{started,finished}" >> ${outputFolder}/scripts/ConvertToFastQ_${name}.sh
done

for i in $(ls ${outputFolder}/scripts/ConvertToFastQ_*.sh)
do
	sbatch $i
done

echo "follow the progress in: ${outputFolder}/scripts/"
