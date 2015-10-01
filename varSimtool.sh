module load varsim/0.5.1
module load jdk/1.8.0_25
module load R
module load ngs-utils
module list

### The best way to run this script is to the following:
### Step 1: run intersect.sh first if the comparison is between Nextgene and a pipeline (NextGene is calling whole genome)
### Step 2: make folder somewhere that will be HOME
### Step 3: Make 2 directories with vcfs of TRUTH and COMPARE


HOME=""
TRUTH="NextGene"
COMPARE="3.1.2"
OUTPUT=""

if [ ! -d $OUTPUT ] 
then
	mkdir $OUTPUT
fi

THISDIR=`pwd`

##LISTOFSAMPLES should be a list of samples seperated by a space 
for SAMPLE in $LISTOFSAMPLES
do
	echo "S1:${TRUTH}/intersect/${SAMPLE}*.intersect.vcf"  
	echo "S2:${COMPARE}/${SAMPLE}*.intersect.vcf"
	mkdir -p ${OUTPUT}/${SAMPLE}/

	cd /gcc/tools/varsim_0.5.1
	java -jar VarSim.jar vcfcompare -true_vcf ${HOME}/${TRUTH}/intersect/${SAMPLE}*intersect.vcf \
	-prefix ${OUTPUT}/${SAMPLE}/${SAMPLE}_${TRUTH}_VS_${COMPARE} \
	${HOME}/${COMPARE}/intersect/${SAMPLE}*intersect.vcf > ${OUTPUT}/${SAMPLE}/${SAMPLE}_${TRUTH}_VS_${COMPARE}.stats 2>&1
	
	echo "${OUTPUT}/${SAMPLE}/${SAMPLE}_${TRUTH}_VS_${COMPARE}.stats"
	perl -pi -e 's|Num new variants read.*|Num new variants read \nALLVARIANTS|' ${OUTPUT}/${SAMPLE}/${SAMPLE}_${TRUTH}_VS_${COMPARE}.stats 	
	
	Rscript ${EBROOTNGSMINUTILS}/make_plots_v2.R 'ALLVARIANTS' ${OUTPUT}/${SAMPLE}/${SAMPLE}_${TRUTH}_VS_${COMPARE}.stats "${TRUTH} vs ${COMPARE}"
	Rscript ${EBROOTNGSMINUTILS}/make_plots_v2.R 'SNP' ${OUTPUT}/${SAMPLE}/${SAMPLE}_${TRUTH}_VS_${COMPARE}.stats "${TRUTH} vs ${COMPARE}"
	Rscript ${EBROOTNGSMINUTILS}/make_plots_v2.R 'Insertion' ${OUTPUT}/${SAMPLE}/${SAMPLE}_${TRUTH}_VS_${COMPARE}.stats "${TRUTH} vs ${COMPARE}"
	Rscript ${EBROOTNGSMINUTILS}/make_plots_v2.R 'Deletion' ${OUTPUT}/${SAMPLE}/${SAMPLE}_${TRUTH}_VS_${COMPARE}.stats "${TRUTH} vs ${COMPARE}"

done

