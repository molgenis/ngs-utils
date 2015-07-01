#!/bin/bash

set -u
set -e
underline=`tput smul`
normal=`tput sgr0`
bold=`tput bold`

function usage () {
echo "
${bold}This script is a small pipeline of 3 seperate scripts creating:
- an interval_list 
- per chromosome bed file 
- optionally a coverage_per_base bed file and interval_list ${normal}
	
These scripts are also available to run seperately, please run the scripts below from the /gcc/tools/scripts folder:
To make interval lists          --> create_interval_listV4.pl
To split bed into batches	--> scatter_gather.sh
Make coverage per base file	--> coverage_per_base.sh	

${bold}Arguments${normal}
        Required:
        -n|--name              BED name (without extension)
	
	Optional:
	-c|--coverageperbase	true or false (default: false)
	-d|--data		What kind of data. If exome: batchsize is automatically 100 , wgs=200 (default: targeted = batchsize of 25)
				choose between exome (92 autosomal, 7X, 1Y), targeted (22 autosomal, 2X, 1Y) or wgs(185 autosomal, 14X, 1Y) 
	-e|--extension		extension to the bed file (default: human_g1k_v37)
	-o|--intervalfolder	path to intervalfolder (default: /gcc/resources/b37/intervals)
	-r|--reference 		Which reference file is used (default: /gcc/resources/b37/indices/human_g1k_v37.dict)"
}

PARSED_OPTIONS=$(getopt -n "$0"  -o n:o:e:r:c:d: --long "name:,intervalfolder:extension:reference:coverageperbase:data"  -- "$@")

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
        -n|--name)
                case "$2" in
                "") shift 2 ;;
                *) NAME=$2 ; shift 2 ;;
            esac ;;
	-c|--coverageperbase)
                case "$2" in
                *) COVPERBASE=$2 ; shift 2 ;;
            esac ;;
	-d|--data)
                case "$2" in
                *) DATA=$2 ; shift 2 ;;
            esac ;;
	-o|--INTERVALFOLDER)
                case "$2" in
                *) INTERVALFOLDER=$2 ; shift 2 ;;
            esac ;;
	-e|--extension)
                case "$2" in
                *) EXTENSION=$2 ; shift 2 ;;
            esac ;;
	-r|--reference)
                case "$2" in
                *) REFERENCE=$2 ; shift 2 ;;
            esac ;;
        --) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
  esac
done

#
# Check required options were provided.
if [[ -z "${NAME-}" ]]; then
        usage
        echo "FATAL: missing required parameter."
        exit 1
fi
if [[ -z "${INTERVALFOLDER-}" ]]; then
	INTERVALFOLDER="/gcc/resources/b37/intervals/"
fi
if [[ -z "${EXTENSION-}" ]]; then
        EXTENSION="human_g1k_v37"
fi
if [[ -z "${REFERENCE-}" ]]; then
        REFERENCE="/gcc/resources/b37/indices/${EXTENSION}"
fi
if [[ -z "${COVPERBASE-}" ]]; then
        COVPERBASE="false"
fi
if [[ -z "${DATA-}" ]]; then
        DATA="targeted"
fi

BATCHCOUNT=22
phiXRef=${REFERENCE}_phiX.dict
phiXExt=${EXTENSION}_phiX
batchCount_X=2


#check which data
if [ $DATA == "wgs" ]
then
	BATCHCOUNT=185
	batchCount_X=14
elif [ $DATA == "exome" ]
then 
	BATCHCOUNT=92
	batchCount_X=7		
fi

echo "NAME: $NAME"
echo "INTERVALSFOLDER: $INTERVALFOLDER"
echo "EXTENSION: $EXTENSION"
echo "REFERENCE: $REFERENCE"
echo "COVPERBASE: $COVPERBASE"
echo "BATCHCOUNT: $((BATCHCOUNT + 1 + batchCount_X))"

MAP="${INTERVALFOLDER}/${NAME}/"

if [[ ${NAME} == *"baits"*  ]] || [[ ${NAME} == *"v37"* ]] || [[ ${NAME} == *"exons"* ]] || [[ ${NAME} == *".bed"* ]]
then
        echo "No need to put extension (baits, exons, v37 or .bed) in the name, only the name of the bed"
        exit 0
fi

if [ ! -d ${MAP} ]
then
        mkdir -p ${MAP}/
	echo "created ${MAP}"
else
	rm -rf ${MAP}
	mkdir -p ${MAP}/
	echo "removed and created ${MAP}"
fi

cp /gcc/resources/b37/intervals/1000G_phase1.indels_Mills_and_1000G_gold_standard.indels.b37.human_g1k_v37.* ${MAP}



baits=${MAP}/${NAME}_baits_${phiXExt}
exons=${MAP}/${NAME}_exons_${phiXExt}

#make copy of the existing BAIT bed file
cp ${INTERVALFOLDER}/${NAME}_baits_${EXTENSION}.bed ${baits}.bed.tmp
if [ ! -f ${INTERVALFOLDER}/${NAME}_exons_${EXTENSION}.bed ]
then
	cp ${INTERVALFOLDER}/${NAME}_baits_${EXTENSION}.bed ${exons}.bed.tmp
else
	cp ${INTERVALFOLDER}/${NAME}_exons_${EXTENSION}.bed ${exons}.bed.tmp
fi	

## If there are 4 columns, it adds an extra column (this is necessary for the GATK batch tool
colcount=`awk '{print NF}' ${INTERVALFOLDER}/${NAME}_baits_${EXTENSION}.bed | sort | tail -n 1`

#check for the presence of phiX region
if [ -f ${INTERVALFOLDER}/${NAME}_baits_${phiXExt}.bed ] 
then
	a=`grep phiX174 ${INTERVALFOLDER}/${NAME}_baits_${phiXExt}.bed | wc -l`
else
	a=0
fi

if [ -f ${INTERVALFOLDER}/${NAME}_exons_${phiXExt}.bed ]
then
	b=`grep phiX174 ${INTERVALFOLDER}/${NAME}_exons_${phiXExt}.bed | wc -l`
else
	b=0
fi

if [ $a == 0 ] 
then
	if [ "$colcount" == "4" ]
	then
		cat ${INTERVALFOLDER}/${NAME}_baits_${EXTENSION}.bed >> ${baits}.bed.tmp
		echo -e 'NC_001422.1\t151\t5236\tphiX174' >> ${baits}.bed.tmp
	else
                echo -e 'NC_001422.1\t151\t5236\t+\tphiX174' >> ${baits}.bed.tmp
	fi
fi

if [ $b == 0 ]
then
	if [ "$colcount" == "4" ]
	then
		cat ${INTERVALFOLDER}/${NAME}_baits_${EXTENSION}.bed >> ${exons}.bed.tmp
		echo -e 'NC_001422.1\t151\t5236\tphiX174' >> ${exons}.bed.tmp
	else
		echo -e 'NC_001422.1\t151\t5236\t+\tphiX174' >> ${exons}.bed.tmp
	fi
fi

echo "colcount: $colcount"
if [ "${colcount}" == "4" ] 
then
	awk '{print $1"\t"$2"\t"$3"\t+\t"$4}' ${baits}.bed.tmp > ${baits}.bed
	awk '{print $1"\t"$2"\t"$3"\t+\t"$4}' ${exons}.bed.tmp > ${exons}.bed
	echo "added strand information"
else
	cp ${baits}.bed.tmp ${baits}.bed
	cp ${exons}.bed.tmp ${exons}.bed
fi

if [ -f ${baits}.withoutChrX.bed ]
then
	rm ${baits}.withoutChrX.bed 
fi

awk '{
	if ($1 != "X"){
        	print $0 >> "'${baits}'.withoutChrX.bed"
        }
}' ${baits}.bed



if [ -f ${baits}.batch-1.bed ] 
then
        echo "splitting in batches skipped"
elif [[ ${NAME} == *"baits"*  ]] || [[ ${NAME} == *"v37"* ]] || [[ ${NAME} == *"exons"* ]] || [[ ${NAME} == *".bed"* ]]
then
	echo "No need to put extension (baits, exons, v37 or .bed) in the name, only the name of the bed"
	exit 0
fi

TMP="/gcc/groups/gcc/tmp03/tmp"
if [ $COVPERBASE == "true" ] 
then
	if [ ! -f ${baits}.uniq.per_base.bed ]
	then 
		echo "starting to create_per_base_intervals, this may take a while"
		perl /gcc/tools/scripts/create_per_base_intervals.pl -input ${baits}.bed -output ${NAME}_baits_${phiXExt} -outputfolder $TMP

		sort -V -k1 -k2 -k3 ${TMP}/${NAME}_baits_${phiXExt}.per_base.bed | uniq -u > ${baits}.uniq.per_base.bed
		rm ${TMP}/${NAME}_baits_${phiXExt}.per_base.bed

		echo "intervals per base done: ${baits}.uniq.per_base.bed"
	else
		echo "${baits}.uniq.per_base.bed already exists, skipped!"
	fi 
	
	#make interval_list coverage per base
	cat ${phiXRef} > ${baits}.uniq.per_base.interval_list
	cat ${baits}.uniq.per_base.bed >> ${baits}.uniq.per_base.interval_list 

fi

if [ -f ${baits}.batch-1.bed ]
then
        echo "Is this bed file already splitted before? If so, please remove the old ones or do not run this script ;)"
else
	module load picard-tools/1.129

	batchIntervallistDir=${MAP}

	chrXNONPARInterval=${baits}.chrX.nonpar.interval_list
	chrXPARBed=${baits}.chrX.par.bed
	AllWithoutchrXInterval=${baits}.withoutChrX.interval_list

	cat ${phiXRef} > ${AllWithoutchrXInterval}
	cat ${phiXRef} > ${chrXNONPARInterval}

	cat ${baits}.withoutChrX.bed >> ${AllWithoutchrXInterval}
	
	awk '{
		if ($1 == "X"){
	       		if (($2 >= 60001  && $3 <= 2699520 ) || ($2 >= 154931044 && $3 <= 155260560 )){
        	                print $0 >> "'${chrXPARBed}'"
        		} else {
        		        print $0 >> "'${chrXNONPARInterval}'"
        		}
		}
	}' ${baits}.bed

	if [ -f ${chrXPARBed} ]
	then
	        cat ${chrXPARBed} >> ${AllWithoutchrXInterval}
	fi

	#autosomal
	java -jar  -Xmx4g -XX:ParallelGCThreads=4 ${PICARD_HOME}/picard.jar IntervalListTools \
	INPUT=${AllWithoutchrXInterval} \
	OUTPUT=${batchIntervallistDir} \
	PADDING=150 \
	UNIQUE=true \
	SCATTER_COUNT=${BATCHCOUNT}

	echo "AUTOSOMAL DONE"
	#non PAR region
	java -jar  -Xmx4g -XX:ParallelGCThreads=4 ${PICARD_HOME}/picard.jar IntervalListTools \
     	INPUT=${chrXNONPARInterval} \
     	OUTPUT=${batchIntervallistDir} \
     	PADDING=150 \
     	UNIQUE=true \
     	SCATTER_COUNT=${batchCount_X} \

	echo "PAR DONE"
	BATCH_ALL=$((BATCHCOUNT + batchCount_X))
	#move the X chromosome folders
	lengthR=`less ${phiXRef} | wc -l`
	echo "lengthR: $lengthR"
	lengthRef=$(( ${lengthR} + 2 ))
	for i in $(seq 1 ${batchCount_X})
	do
		bi=$(( BATCHCOUNT + i  ))
		ba=${baits}.batch-${bi}X
		echo "ba=$ba bi=$bi"
		if [[ ${i} -lt 10 ]]
               	then
			echo "$i is minder dan 10"
			mv  ${batchIntervallistDir}/temp_000${i}_of_${batchCount_X}/scattered.intervals  ${ba}.interval_list 
			tail -n+${lengthRef} ${ba}.interval_list > ${ba}.bed
		else
			echo "$i is meer dan 10"
			mv  ${batchIntervallistDir}/temp_00${i}_of_${batchCount_X}/scattered.intervals  ${ba}.interval_list
                        tail -n+${lengthRef} ${ba}.interval_list > ${ba}.bed
		fi	
	done

		
	BATCH_Y=$((BATCH_ALL + 1))

	cat ${phiXRef} > ${baits}.batch-${BATCH_Y}Y.interval_list
	#MOVING ALL INTERVAL_LIST FILES 
	for i in $(seq 1 ${BATCHCOUNT})
	do
		if [[ ${i} -lt 10 ]]
		then
	                mv ${batchIntervallistDir}/temp_000${i}_of_${BATCHCOUNT}/scattered.intervals ${baits}.batch-${i}.interval_list
	                tail -n+${lengthRef} ${baits}.batch-${i}.interval_list > ${baits}.batch-${i}.bed

		elif [[ ${i} -gt 99 ]]
		then
			mv ${batchIntervallistDir}/temp_0${i}_of_${BATCHCOUNT}/scattered.intervals ${baits}.batch-${i}.interval_list
                        tail -n+${lengthRef} ${baits}.batch-${i}.interval_list > ${baits}.batch-${i}.bed		
        	else
                        mv ${batchIntervallistDir}/temp_00${i}_of_${BATCHCOUNT}/scattered.intervals ${baits}.batch-${i}.interval_list
                        tail -n+${lengthRef} ${baits}.batch-${i}.interval_list > ${baits}.batch-${i}.bed

        	fi
	done

	rm ${baits}.batch*.interval_list

	for i in $(seq $((${BATCHCOUNT}-10)) ${BATCHCOUNT})
	do
	awk '{
	        if($1 == "Y"){
       	        	print $0
        	}

	}' ${baits}.batch-${i}.bed >> ${baits}.batch-${BATCH_Y}Y.bed
	done

	for i in $(seq $((${BATCHCOUNT}-10)) ${BATCHCOUNT})
	do
		sed '/^Y/ d' ${baits}.batch-${i}.bed > ${baits}.batch-${i}.bed.tmp ; mv ${baits}.batch-${i}.bed.tmp ${baits}.batch-${i}.bed
	done

	echo ""
	for i in $(ls ${MAP}/*batch-*.bed); do cat $i | awk -v var="$i" '{if( $2==$3){print var}}';done > ${MAP}/chompLines.txt	

	while read line
	do
		awk '{if ($2!=$3){ print $0}}' $line > ${line}.tmp
		mv ${line}.tmp $line 
	done<${MAP}/chompLines.txt
	

	echo "batching complete"
	rm -rf ${batchIntervallistDir}/temp_0*
fi

for f in ${MAP}/*_baits_*; do cp $f ${f/_baits_/_exons_}; done

if [ -f ${baits}.interval_list ]
then
	echo "interval_list already exists, skipping"
else
	perl create_interval_listV4.pl -Ref ${phiXRef}  -Exons ${exons}  -Baits ${baits}
	echo "intervals created"
fi
