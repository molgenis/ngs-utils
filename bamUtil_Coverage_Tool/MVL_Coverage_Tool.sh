#!/bin/bash
set -u
set -e
underline=`tput smul`
normal=`tput sgr0`
bold=`tput bold`

function usage () {
echo "
${bold}This script is calculating coverage per position of pathogenic and likely pathogenic varaints from the MVL.
AvgCoverage, SD, avarage cov <1x,<10x,<30x,
${bold}Arguments${normal}
	example sh MVL_Coverage_Tool.sh -s /groups/umcg-gd/prm02/projects/QXTR_90-Exoom_v1/run01/results/alignment/ OPTIONS
	Required:
	-s|--sampledir		Directory of of bam files of the used samples
	-t|--tmpdir         Give tmpfolder location for inbetween and final files
	Optional:
	-b|--bamstat		path to the ban stats (default: $EBROOTNGSMINUTILS/bamUtil_Coverage_Tool/bamUtil-1.0.14/bin/bam stats)
	-m|--mvl	        location of MVL_per_bp_unique.txt (default: $EBROOTNGSMINUTILS/bamUtil_Coverage_Tool/MVL_per_bp_unique.txt)
	-x|--mvltabix       location of MVL_per_bp_unique_tabix.txt (default: $EBROOTNGSMINUTILS/bamUtil_Coverage_Tool/MVL_per_bp_unique_tabix.txt)
	-q|--qcmvldir          location of MerdgeQCofMVL_file.py necessary for final txt generation (default:$EBROOTNGSMINUTILS/bamUtil_Coverage_Tool/)
"
}

module load HTSlib
module load Python

PARSED_OPTIONS=$(getopt -n "$0"  -o s:b:o:m:x:t:f:z: --long "sampledir:bamstat:outputdir:mvl:mvltabix:tmpdir:outputfilter:finaloutput:"  -- "$@")

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
    -s|--sampledir)
                case "$2" in
                "") shift 2 ;;
                *) SAMPLEDIR=$2 ; shift 2 ;;
            esac ;;
	-b|--bamstat)
                case "$2" in
                *) BAMSTAT=$2 ; shift 2 ;;
            esac ;;
	-m|--mvl)
                case "$2" in
                *) MVL=$2 ; shift 2 ;;
            esac ;;
	-x|--mvltabix)
                case "$2" in
                *) MVLTABIX=$2 ; shift 2 ;;
            esac ;;
    -t|--tmpdir)
                case "$2" in
                *) TMPDIR=$2 ; shift 2 ;;
            esac ;; 
    -q|--qcmvl)
                case "$2" in
                *) QCMVLDIR=$2 ; shift 2 ;;
            esac ;; 
        --) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
  esac
done



empty=""
#
#Check required options were provided.
#
if [[ -z "${SAMPLEDIR-}" ]]; then
	usage
	exit 1
fi
if [[ -z "${BAMSTAT-}" ]]; then
	BAMSTAT="$EBROOTNGSMINUTILS/bamUtil_Coverage_Tool/bamUtil-1.0.14/bin/bam"
fi
if [[ -z "${MVL-}" ]]; then
        MVL="$EBROOTNGSMINUTILS/bamUtil_Coverage_Tool/MVL_per_bp_unique.txt"
fi
if [[ -z "${MVLTABIX-}" ]]; then
        MVLTABIX="$EBROOTNGSMINUTILS/bamUtil_Coverage_Tool/MVL_per_bp_unique_tabix.txt"
fi
if [[ -z "${TMPDIR-}" ]]; then
	usage
	exit 1
fi
if [[ -z "${QCMVLDIR-}" ]]; then
        QCMVLDIR="$EBROOTNGSMINUTILS/bamUtil_Coverage_Tool/"
fi
echo "SAMPLEDIR: ${SAMPLEDIR}"
echo "BAMSTAT: ${BAMSTAT}"
echo "MVL: ${MVL}"
echo "MVLTABIX: ${MVLTABIX}"
echo "TMPDIR: ${TMPDIR}"
echo "QCMVLDIR: ${QCMVLDIR}"
echo "starting"


rm -rf ${TMPDIR}
mkdir -p ${TMPDIR}/input/
mkdir -p ${TMPDIR}/outputfilter/
mkdir -p ${TMPDIR}/output_final/
mkdir -p ${TMPDIR}/output_final_txt_files/


#runBamStatWithRegion3.sh run the BAM stat tool to generate QC report per position
echo "start Bam stats"

count=0

for i in $(ls ${SAMPLEDIR}/*.bam)
do
	echo "processingBamstat ${i}"
	${BAMSTAT} stats \
	--in ${i} \
	--cBaseQC ${TMPDIR}/input/outputBamUtil_${count}.txt \
	--regionlist ${MVL} \
	--withinRegion \
	--nophonehome

	    count=$((count+1))

done
echo "finisched Bam stats"


# test.sh, make QC report complete by adding also position with 0x coverage

echo "start tabix 0x coverage"
for i in $(ls ${TMPDIR}/input/output*.txt)
do 
    b=$(basename ${i})
    bedfile=${b%.*}.bed.txt
    awk 'BEGIN{OFS="\t"}{print $1,$2,($3-1),$4}' ${i} > ${TMPDIR}/${bedfile}

sed -i '1d' ${TMPDIR}/${bedfile}    
bgzip -c ${TMPDIR}/${bedfile} > ${TMPDIR}/${bedfile}.gz
    tabix -0 -p bed ${TMPDIR}/${bedfile}.gz
rm -f ${TMPDIR}/outputfilter/${b%.*}_filter.txt    
    while read line
    do 
        if tabix -0 ${TMPDIR}/${bedfile}.gz ${line} | grep -q "^"
        then
            tabix -0 ${TMPDIR}/${bedfile}.gz ${line}
        else
            echo $line | awk -F":" '{print $1"\t"$2}' | awk -F"-" '{print $1"\t"$2"\t0"}' 
        fi
    done < ${MVLTABIX} >> ${TMPDIR}/outputfilter/${b%.*}_filter.txt
echo "${i} done tabix 0x coverage"
done

# test_chr_pos_cov.sh, make files usable for python script including 0x 
echo "starting making final txt files per patient"
for i in $(ls ${TMPDIR}/outputfilter/output*.txt)
do
        b=$(basename ${i})
        filename=${b%.*}

        awk '{
                        if ($4==0) { 
                                print $1"\t"($2)"\t"$4
                        }
                        else{
                                print $1"\t"$2"\t"$4
                        }
                }' $i > "${TMPDIR}/output_final/${filename}_final.txt" 
echo "${i} done"
done
echo "finisched making final txt per patient"

#add gene name to coverage file. 
echo "add gene name"
sort -V ${MVL} > ${MVL}.sorted
awk '{print $4}' ${MVL} > ${MVL}.genes

for i in $(ls ${TMPDIR}/output_final/outputBamUtil_*_filter_final.txt)
do

    sortedFile=${i%.*}.sorted.txt
    sort -V ${i} > ${sortedFile}
    a=$(diff <(awk '{print $2}' "${MVL}.sorted") <(awk '{print $2}' "${sortedFile}"))
    if [ "${a}" == "" ]
    then
        paste -d"\t" $i ${MVL}.genes > ${i%.*}_gene.txt
    else 
        echo 'positions in MVL and BAMutil output do not match'
    fi
echo "${i} done adding gene name"    
done
echo "start python script"

python ${QCMVLDIR}/MerdgeQCofMVL_file.py --in ${TMPDIR}/output_final/ --ff ${TMPDIR}/output_final_txt_files/

