#This script should be changed to be able to use it. Fill in the appropriate input and output dir and give the name of the bedfile without the extension

INPUTDIR=/gcc/resources/b37/intervals
NAME=bedfilename_without.bed
OUTPUTDIR=/gcc/resources/b37/intervals

awk '{print $0 >> "'${OUTPUTDIR}'/'${NAME}'." $1".bed" }' ${INPUTDIR}/${NAME}.bed 
