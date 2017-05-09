#!/bin/bash
set -u
set -e
underline=`tput smul`
normal=`tput sgr0`
bold=`tput bold`

function usage () {
echo "
${bold}This script is calculating coverage calculations per target and per Gene based on coverage_per_target files created by DepthOfCoverage.
AvgCoverage, median, percentage <10x,<20x,<50x,<100x coverage


${bold}Arguments${normal}

	example sh countCoverage.sh -p Exoom_v1 OPTIONS

	Required:
	-p|--panel		Name of the panel/capturingkit, e.g. Exoom_v1, CARDIO_v3
	Optional:
	-w|--workdir		working directory (default: /groups/umcg-gaf/tmp04/coverage/{panel})
	-s|--structure		relative path from permanentDir that contains coveragepertarget files(default: run01/results/coverage/CoveragePerTarget/)
	-d|--permanentdir	location of the permanantDir (default: /groups/umcg-gd/prm02/projects/)
	-t|--tmp		Give tmpfolder location (default: \${workdir}/tmp)"

It will automatically get all files that are in the structure, the following command will be executed:
e.g.
command: ls \${permanentDir}/*\${panel}/\${structure}/*\${panel}*.coveragePerTarget.txt
weaved command: ls /groups/umcg-gd/prm02/projects/*-Exoom_v1/run01/results/coverage/CoveragePerTarget/*Exoom_v1*.coveragePerTarget.txt


}

module load ngs-utils
PARSED_OPTIONS=$(getopt -n "$0"  -o p:w:s:d:t: --long "name:panel:workdir:structure:permanentdir:tmp:"  -- "$@")

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
        -p|--panel)
                case "$2" in
                "") shift 2 ;;
                *) PANEL=$2 ; shift 2 ;;
            esac ;;
	-w|--workdir)
                case "$2" in
                *) WORDKIR=$2 ; shift 2 ;;
            esac ;;
	-s|--structure)
                case "$2" in
                *) STRUCTURE=$2 ; shift 2 ;;
            esac ;;
	-d|--permanentdir)
                case "$2" in
                *) PRMDIR=$2 ; shift 2 ;;
            esac ;;
	-t|--tmp)
                case "$2" in
                *) TMP=$2 ; shift 2 ;;
            esac ;;
        --) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
  esac
done

empty=""
#
# Check required options were provided.
if [[ -z "${PANEL-}" ]]; then
	usage
	exit 1
fi
if [[ -z "${WORKDIR-}" ]]; then
	WORKDIR="/groups/umcg-gaf/tmp04/coverage/${PANEL}/"
fi
if [[ -z "${STRUCTURE-}" ]]; then
        STRUCTURE="run01/results/coverage/CoveragePerTarget/"
fi
if [[ -z "${PRMDIR-}" ]]; then
        PRMDIR="/groups/umcg-gd/prm02/projects/"
fi
if [[ -z "${TMP-}" ]]; then
	TMP=${WORKDIR}/tmp/
	mkdir -p ${TMP}	
	echo "makedir ${TMP}"
fi


rm -rf ${WORKDIR}
mkdir -p ${WORKDIR}/coverage/
mkdir -p ${WORKDIR}/tmp/

echo "starting to extract all the ${PANEL} projects and retreive the coverage for each sample"
echo "ls ${PRMDIR}/*${PANEL}/${STRUCTURE}/*${PANEL}*.coveragePerTarget.txt"
count=0

for i in $(ls ${PRMDIR}/*${PANEL}/${STRUCTURE}/*${PANEL}*.coveragePerTarget.txt)
do
	if [ $count == 0 ]
	then
		awk '{print $2,$3,$4,$6}' $i > ${TMP}/firstcolumns.txt
		count=$((count+1))
	fi
	 
	NAMEFILE=$(basename $i)
	DIRNAME=$(dirname $i)
	SAMPLENAME=${NAMEFILE%%.*}

	#awk '{sum+=$5}END {print sum/210414}' $i > $DIRNAME/$SAMPLENAME.averageCoverage.txt
	awk '{print $5}' $i > ${WORKDIR}/coverage/${SAMPLENAME}.coverage
done

echo "${WORKDIR}/coverage/"

for i in $(ls -d ${WORKDIR})
do
	find $i/coverage/ -type f -name "*.coverage" | xargs paste > ${TMP}/BIG.pasta
done

echo "## calculate MEDIAN ##"
## calculate MEDIAN ##
awk -v max=1 '
        function median(c,v,  j) { 
           asort(v,j); 
           if (c % 2) return j[(c+1)/2]; 
           else return (j[c/2+1]+j[c/2])/2.0; 
        } 
{ 
         count++;values[count]=$NF;  
         if (count >= max) { 
           print  median(count,values); count=0; 
         } 
} 
END { 
         print  median(count,values); 
}' ${TMP}/BIG.pasta > ${TMP}/BIG.pasta.median

echo "## calculate AVG ##"
## CALCULATE AVG
awk '{ for(i = 1; i <= NF; i++) sum+=$i;print sum/NF;sum=0 }' ${TMP}/BIG.pasta > ${TMP}/BIG.pasta.avg

echo "## Calculate percentage under 10,20,50 and 100x ##"
## Calculate percentage under 10,20,50 and 100x ##
awk '{ for(i = 1; i <= NF; i++) if ($i < 10 )counter+=1;print 100-((counter/NF))*100;counter=0 }'  ${TMP}/BIG.pasta > ${TMP}/BIG.pasta.percentageU10
awk '{ for(i = 1; i <= NF; i++) if ($i < 20 )counter+=1;print 100-((counter/NF))*100;counter=0 }'  ${TMP}/BIG.pasta > ${TMP}/BIG.pasta.percentageU20
awk '{ for(i = 1; i <= NF; i++) if ($i < 50 )counter+=1;print 100-((counter/NF))*100;counter=0 }'  ${TMP}/BIG.pasta > ${TMP}/BIG.pasta.percentageU50
awk '{ for(i = 1; i <= NF; i++) if ($i < 100 )counter+=1;print 100-((counter/NF))*100;counter=0 }'  ${TMP}/BIG.pasta > ${TMP}/BIG.pasta.percentageU100

awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' ${TMP}/firstcolumns.txt > ${TMP}/first4columns.txt

echo "## update column gene with only	one annotation possible"
## update column gene with only one annotation possible
awk '{print $4}' ${TMP}/first4columns.txt | awk '{FS=",|:"}{print $1}' > ${TMP}/UpdatedGenes.txt
awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' ${TMP}/first4columns.txt > ${TMP}/chromstartstop.txt
paste ${TMP}/chromstartstop.txt  ${TMP}/UpdatedGenes.txt > ${TMP}/BigupdatedFile.txt

echo "pasting median,avg,10x,20x,50x and 100x"
rm -f ${WORKDIR}/CoverageOverview.txt

paste ${TMP}/BigupdatedFile.txt ${TMP}/BIG.pasta.median > ${TMP}/BigupdatedFileMedian.txt
paste ${TMP}/BigupdatedFileMedian.txt ${TMP}/BIG.pasta.avg > ${TMP}/BigupdatedFileAvg.txt 

paste ${TMP}/BigupdatedFileAvg.txt ${TMP}/BIG.pasta.percentageU10 > ${TMP}/BigupdatedFile10.txt  
paste ${TMP}/BigupdatedFile10.txt ${TMP}/BIG.pasta.percentageU20 > ${TMP}/BigupdatedFile20.txt 
paste ${TMP}/BigupdatedFile20.txt ${TMP}/BIG.pasta.percentageU50 > ${TMP}/BigupdatedFile50.txt 
paste ${TMP}/BigupdatedFile50.txt ${TMP}/BIG.pasta.percentageU100 > ${TMP}/BigupdatedFile100.txt 

echo -e "Chr\tStart\tStop\tGene\tMedian\tAvgCoverage\tu10\tu20\tu50\tu100" > ${WORKDIR}/CoverageOverview.txt
tail -n+2 ${TMP}/BigupdatedFile100.txt >> ${WORKDIR}/CoverageOverview.txt 
head -n -1 ${WORKDIR}/CoverageOverview.txt > ${TMP}/CoverageOverview.txt.tmp
cp ${TMP}/CoverageOverview.txt.tmp ${WORKDIR}/CoverageOverview.txt
echo "DONE, final file is ${WORKDIR}/CoverageOverview.txt"

echo "now creating per Gene calculations"
python countCoveragePerGene.py > ${WORKDIR}/CoverageOverview_PerGene.txt

echo -e "Gene\tAvgCoverage\tNo of Targets\tMedian\tu10x\tu20x\tu50x\tu100x" > ${WORKDIR}/CoverageOverview_PerGene.sorted.txt

sort -V -k1 ${WORKDIR}/CoverageOverview_PerGene.txt >> ${WORKDIR}/CoverageOverview_PerGene.sorted.txt
echo "done, gene file can be found here: ${WORKDIR}/CoverageOverview_PerGene.sorted.txt"
