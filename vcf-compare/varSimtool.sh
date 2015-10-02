set -u
set -e

underline=`tput smul`
normal=`tput sgr0`
bold=`tput bold`

function usage () {
echo "
${bold}
### Make intersect file  + runVarsim tool
Step 1: make a WORKDIR, this will be the working directory (e.g. /my/path/to/workdir/)
Step 2: copy the bed file that you need to your WORKDIR and call it captured.bed
Step 3: make 2 directories in the WORKDIR and put the vcf files in the different directories (e.g. /my/path/to/workdir/3.1.2)

Arguments${normal}
        Required:
        -w|--workdir		This will be your working directory
        -c|--compare		What is the name of the project which is to be compared

        Optional:
        -o|--outputfolder     	Path to outputfolder (default: WORKDIR/output
        -t|--truth    		what is the name of the project which is considered as the truth (default: NextGene) "
}

PARSED_OPTIONS=$(getopt -n "$0"  -o w:o:c:t: --long "workdir:,outputfolder:compare:truth:"  -- "$@")

#
# Bad arguments, something has gone wrong with the getopt command.
#
if [ $? -ne 0 ]; then
        usage
        echo "FATAL: Wrong arguments."
        exit 1
fi

eval set -- "$PARSED_OPTIONS"

#
# Now goes through all the options with a case and using shift to analyse 1 argument at a time.
# $1 identifies the first argument, and when we use shift we discard the first argument, so $2 becomes $1 and goes again through the case.
#

while true; do
  case "$1" in
        -w|--workdir)
                case "$2" in
                "") shift 2 ;;
                *) WORKDIR=$2 ; shift 2 ;;
            esac ;;
        -c|--compare)
                case "$2" in
                *) COMPARE=$2 ; shift 2 ;;
            esac ;;
        -o|--output)
                case "$2" in
                *) OUTPUT=$2 ; shift 2 ;;
            esac ;;
        -t|--truth)
                case "$2" in
                *) TRUTH=$2 ; shift 2 ;;
            esac ;;
        --) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
  esac
done

#
# Check required options were provided.
if [[ -z "${WORKDIR-}" ]]; then
        usage
        echo "FATAL: missing required parameter ${bold} workdir."${normal}
        exit 1
fi
if [[ -z "${COMPARE-}" ]]; then
        usage
        echo "FATAL: missing required parameter ${bold} compare."${normal}
        exit 1
fi
if [[ -z "${OUTPUT-}" ]]; then
	OUTPUT="${WORKDIR}/output"
fi
if [[ -z "${TRUTH-}" ]]; then
        TRUTH="NextGene"
fi

module load varsim
module load R
module load ngs-utils
module load BEDTools
module list

### MAKE INTERSECT FILES

if [ ! -d $WORKDIR ]
then
	exit 1
fi
if [ ! -d $WORKDIR/$TRUTH/intersect/ ]
then
	if [ ! -d $WORKDIR/$TRUTH/ ]
	then
		echo "There is no $TRUTH directory in your WORKDIR"
		exit 1
	else
	    	mkdir $WORKDIR/$TRUTH/intersect/
	fi
fi

if [ ! -d $WORKDIR/$COMPARE/intersect/ ]
then
	if [ ! -d $WORKDIR/$COMPARE/ ]
        then
            	echo "There is no $COMPARE directory in your WORKDIR "
                exit 1
        else
                mkdir $WORKDIR/$COMPARE/intersect/
        fi
fi

echo "starting with intersect"

cd ${WORKDIR}/${TRUTH}
for i in $(ls *.vcf)
do
bedtools intersect -a $i -b ${WORKDIR}/captured.bed > $WORKDIR/$TRUTH/intersect/${i}.intersect.vcf
echo "$i done"
done
echo ""
echo "TRUTH done"
echo ""

cd ${WORKDIR}/${COMPARE}
for i in $(ls *.vcf)
do
	bedtools intersect -a ${i} -b ${WORKDIR}/captured.bed > $WORKDIR/$COMPARE/intersect/$i.intersect.vcf
	echo "$i done"
done

echo ""
echo "COMPARE done"
echo ""

### Varsim part ###

if [ ! -d $OUTPUT ] 
then
	mkdir $OUTPUT
fi

THISDIR=`pwd`

cd  ${WORKDIR}/${COMPARE}/
LISTOFSAMPLES=($(for i in $(ls -1 *vcf); do echo $i | sed -e 's/.final.vcf//';done))

cd $THISDIR

for SAMPLE in ${LISTOFSAMPLES[@]}
do
	echo "S1:${TRUTH}/intersect/${SAMPLE}*.intersect.vcf"
	echo "S2:${COMPARE}/${SAMPLE}*.intersect.vcf"
	mkdir -p ${OUTPUT}/${SAMPLE}/

	java -jar ${EBROOTVARSIM}/VarSim.jar vcfcompare -true_vcf ${WORKDIR}/${TRUTH}/intersect/${SAMPLE}*intersect.vcf \
	-prefix ${OUTPUT}/${SAMPLE}/${SAMPLE}_${TRUTH}_VS_${COMPARE} \
	${WORKDIR}/${COMPARE}/intersect/${SAMPLE}*intersect.vcf > ${OUTPUT}/${SAMPLE}/${SAMPLE}_${TRUTH}_VS_${COMPARE}.stats 2>&1

	echo "${OUTPUT}/${SAMPLE}/${SAMPLE}_${TRUTH}_VS_${COMPARE}.stats"
	perl -pi -e 's|Num new variants read.*|Num new variants read \nALLVARIANTS|' ${OUTPUT}/${SAMPLE}/${SAMPLE}_${TRUTH}_VS_${COMPARE}.stats

	Rscript ${EBROOTNGSMINUTILS}/vcf-compare/make_vcfcompare_plots.R 'ALLVARIANTS' ${OUTPUT}/${SAMPLE}/${SAMPLE}_${TRUTH}_VS_${COMPARE}.stats "${TRUTH} vs ${COMPARE}"
	Rscript ${EBROOTNGSMINUTILS}/vcf-compare/make_vcfcompare_plots.R 'SNP' ${OUTPUT}/${SAMPLE}/${SAMPLE}_${TRUTH}_VS_${COMPARE}.stats "${TRUTH} vs ${COMPARE}"
	Rscript ${EBROOTNGSMINUTILS}/vcf-compare/make_vcfcompare_plots.R 'Insertion' ${OUTPUT}/${SAMPLE}/${SAMPLE}_${TRUTH}_VS_${COMPARE}.stats "${TRUTH} vs ${COMPARE}"
	Rscript ${EBROOTNGSMINUTILS}/vcf-compare/make_vcfcompare_plots.R 'Deletion' ${OUTPUT}/${SAMPLE}/${SAMPLE}_${TRUTH}_VS_${COMPARE}.stats "${TRUTH} vs ${COMPARE}"

done

