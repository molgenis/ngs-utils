#This script should be changed to be able to use it. Fill in the appropriate input and output dir and give the name of the bedfile without the extension

INPUTDIR=/gcc/groups/gcc/home/mdijkstra/development/InSilicoData/data/humanPhiX/intervals
#INPUTDIR=/gcc/groups/gcc/home/gvdvries/5GPM_bedfiles
NAME=None_exons_b37_human_g1k_v37_phiX
OUTPUTDIR=/gcc/groups/gcc/home/mdijkstra/development/InSilicoData/data/humanPhiX/intervals
#OUTPUTDIR=/gcc/groups/gcc/home/gvdvries/5GPM_bedfiles

if [ -f ${OUTPUTDIR}/${NAME}.chr1.bed ]
	then
	echo "Is this bed file already splitted before? If so, please remove the old ones or do not run this script ;)"	
else
	echo "splitting bed file"
	awk '{print $0 >> "'${OUTPUTDIR}'/'${NAME}'.chr" $1".bed" }' ${INPUTDIR}/${NAME}.bed 
	awk '{
	if ($1 != "X"){
		print $0 >> "'${OUTPUTDIR}'/'${NAME}'.withoutChrX.bed"
	}' ${INPUTDIR}/${NAME}.bed
fi
