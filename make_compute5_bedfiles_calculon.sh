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
	-n|--name              	BED name (${bold}SHOULD ALWAYS BE the name 'captured')${normal}
	Optional:

	-c|--coverageperbase	true or false (default: false)
	-d|--data               What kind of data. exome (10 batches) or wgs (20 batches) or chr (per chromosome,splitted X into 2 for par and nonpar region). (default: targeted = batchsize of 6 (3 + 2X + 1Y)
                                in case of wgs samples there should be no padding (this will automatically skipped when wgs is chosen)
	-e|--extension		extension to the bed file (default: human_g1k_v37)
	-o|--intervalfolder	path to intervalfolder (default: this folder)
	-r|--reference 		Which reference file is used (default: /apps/data/1000G/phase1/human_g1k_v37_phiX.dict)
	-t|--tmp		Give tmpfolder location (default: /groups/umcg-gaf/tmp04/tmp)"
}

module load ngs-utils
PARSED_OPTIONS=$(getopt -n "$0"  -o n:o:e:r:c:d:t: --long "name:intervalfolder:extension:reference:coverageperbase:data:tmp:"  -- "$@")

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
if [[ -z "${NAME-}" ]]; then
	usage
	exit 1
fi
if [[ -z "${INTERVALFOLDER-}" ]]; then
	INTERVALFOLDER="./"	
fi
if [[ -z "${EXTENSION-}" ]]; then
        EXTENSION="human_g1k_v37"
fi
if [[ -z "${REFERENCE-}" ]]; then
        REFERENCE="/apps/data/1000G/phase1/${EXTENSION}"
fi
if [[ -z "${COVPERBASE-}" ]]; then
        COVPERBASE="false"
fi
if [[ -z "${DATA-}" ]]; then
        DATA="targeted"
fi
if [[ -z "${TMP-}" ]]; then
	whichHost=$(hostname)
	if [ "$whichHost" == "zinc-finger.gcc.rug.nl" ]
	then
		 TMP="/groups/umcg-gaf/tmp05/tmp"
	elif [ "$whichHost" == "calculon" ]
	then
		TMP="/groups/umcg-gaf/tmp04/tmp"
	else
		echo "unknown host!"
		exit 1
	fi
fi

BATCHCOUNT=3
phiXRef=${REFERENCE}_phiX.dict
phiXExt=${EXTENSION}_phiX
batchCount_X=2

#check which data
if [ $DATA == "wgs" ]
then
        BATCHCOUNT=17
	echo "BATCHCOUNT: $((BATCHCOUNT + 1 + batchCount_X))"
elif [ $DATA == "exome" ]
then
        BATCHCOUNT=7
	echo "BATCHCOUNT: $((BATCHCOUNT + 1 + batchCount_X))"
elif [ $DATA == "chr" ]
then
	awk '{print $1}' ${NAME}.bed | sort | uniq > countChr.tmp
	BATCHCOUNT=$(cat countChr.tmp | wc -l)
	echo "BATCHCOUNT: $BATCHCOUNT"
else
	echo "BATCHCOUNT: $((BATCHCOUNT + 1 + batchCount_X))"
fi

echo "NAME: $NAME"
echo "INTERVALSFOLDER: $INTERVALFOLDER"
echo "EXTENSION: $EXTENSION"
echo "REFERENCE: $REFERENCE"
echo "COVPERBASE: $COVPERBASE"

MAP="${INTERVALFOLDER}"

if [[ ${NAME} == *"baits"*  ]] || [[ ${NAME} == *"v37"* ]] || [[ ${NAME} == *"exons"* ]] || [[ ${NAME} == *".bed"* ]]
then
        echo "No need to put extension (baits, exons, v37 or .bed) in the name, only the name of the bed"
        exit 0
fi

cp /apps/data/1000G/phase1/Mills_and_1000G_gold_standard/1000G_phase1.indels_Mills_and_1000G_gold_standard.indels.b37.human_g1k_v37.* ${MAP}

baits=${MAP}/${NAME}

#sort -V ${baits}.bed > ${baits}.bed.sorted
#mv ${baits}.bed.sorted ${baits}.bed

## If there are 4 columns, it adds an extra column (this is necessary for the GATK batch tool
colcount=`awk '{print NF}' ${baits}.bed | sort | tail -n 1`

#check for the presence of phiX region
if [ -f ${baits}.bed ]
then
	a=`grep phiX174 ${baits}.bed | wc -l`
else
	a=0
fi
cat ${baits}.bed > ${baits}.bed.tmp

if [ $a == 0 ]
then
	if [ "$colcount" == "4" ]
	then
		echo -e 'NC_001422.1\t151\t5236\tphiX174' >> ${baits}.bed.tmp
	else
	        echo -e 'NC_001422.1\t151\t5236\t+\tphiX174' >> ${baits}.bed.tmp
	fi
fi

echo "colcount: $colcount"
if [ "${colcount}" == "4" ]
then
	awk '{print $1"\t"$2"\t"$3"\t+\t"$4}' ${baits}.bed.tmp > ${baits}.bed
	echo "added strand information"
else
	mv ${baits}.bed.tmp ${baits}.bed
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

module load ngs-utils

if [ $COVPERBASE == "true" ]
then
	if [ ! -f ${baits}.uniq.per_base.bed ]
	then
		echo "starting to create_per_base_intervals, this may take a while"
		create_per_base_intervals.pl -input ${baits}.bed -output ${NAME} -outputfolder $TMP

		sort -V -k1 -k2 -k3 ${TMP}/${NAME}.per_base.bed | uniq -u > ${baits}.uniq.per_base.bed
		rm ${TMP}/${NAME}.per_base.bed

		echo "intervals per base done: ${baits}.uniq.per_base.bed"
	else
		echo "${baits}.uniq.per_base.bed already exists, skipped!"
	fi
	if [ ! -f ${baits}.genesOnly ]
	then
		awk '{print $5}' ${baits}.bed > ${baits}.genesOnly
	fi
	#make interval_list coverage per base
	cat ${phiXRef} > ${baits}.uniq.per_base.interval_list
	cat ${baits}.uniq.per_base.bed >> ${baits}.uniq.per_base.interval_list 
	echo "${baits}.uniq.per_base.interval_list created"

	awk '{ if ($0 !~ /^@/){
		minus=($2 + 1);
		print $1"\t"minus"\t"$3"\t"$4"\t"$5
	}
	else{
		print $0
	}}' ${baits}.uniq.per_base.interval_list > ${baits}.uniq.per_base.interval_list.tmp

	mv ${baits}.uniq.per_base.interval_list.tmp ${baits}.uniq.per_base.interval_list

fi

if [ "${DATA}" == "chr" ]
then

	if [ -f ${baits}.batch-1.bed ]
	then
    		echo "Is this bed file already splitted before? If so, please remove the old ones or do not run this script ;)"
	else
    		batchIntervallistDir=${MAP}

        	chrXNONPARBed=${baits}.batch-Xnp.bed
        	chrXPARBed=${baits}.batch-Xp.bed

		##Lastline is always phiX, we want to know whi
		LASTLINE=$(tail -n1 ${baits}.bed)
		ONEBEFORELASTLINE=$(tail -n2 ${baits}.bed | head -1)
		OLDIFS=$IFS
		set $ONEBEFORELASTLINE
		chromo=$1
		position=$2

        	awk '{
        	      	if ($1 == "X"){
				if (($2 == 1) && ($3 == 155270560)){
					print "X\t60001\t2699520\t+\tWGS" > "'${chrXPARBed}'" 
					print "X\t154931044\t155260560\t+\tWGS" > "'${chrXPARBed}'" 
					print "X\t1\t60000\t+\tWGS" > "'${chrXNONPARBed}'"
					print "X\t2699521\t154931043\t+\tWGS" >> "'${chrXNONPARBed}'"
				}	
               	        	else if (($2 >= 60001  && $3 <= 2699520 ) || ($2 >= 154931044 && $3 <= 155260560 )){
                                	print $0 >> "'${chrXPARBed}'"
                        	} else {
	                        	print $0 >> "'${chrXNONPARBed}'"
                        	}
                	}else if ($1 == "NC_001422.1"){
				#do nothing, will added later
			}
			else{
                	      	print $0 >> "captured.batch-"$1".bed"
                	}
        	}' ${baits}.bed
	### Check where to put the phiXref	
	if [ "${chromo}" == "X" ]
	then
		if [[ ${position} -gt 60001 && ${position} -lt 2699520 ]] || [[ $position -gt 154931044 && $position -lt 155260560 ]]
		then
			echo $LASTLINE >> captured.batch-Xp.bed
		else
			echo $LASTLINE >> captured.batch-Xnp.bed
		fi 
	else
		echo $LASTLINE >> captured.batch-${chromo}.bed
	fi
fi
else
	if [ -f ${baits}.batch-1.bed ]
	then
	        echo "Is this bed file already splitted before? If so, please remove the old ones or do not run this script ;)"
	else
		module load picard/1.130-Java-1.7.0_80

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

		awk '{
		if ($0 !~ /^@/){
        	        minus=($2 + 1);
        	        print $1"\t"minus"\t"$3"\t"$4"\t"$5
        	}
        	else
        	        print $0
        	}' ${AllWithoutchrXInterval} > ${AllWithoutchrXInterval}.tmp
		mv ${AllWithoutchrXInterval}.tmp ${AllWithoutchrXInterval}

	#autosomal
	java -jar -Xmx4g -XX:ParallelGCThreads=4 ${EBROOTPICARD}/picard.jar IntervalListTools \
	INPUT=${AllWithoutchrXInterval} \
	OUTPUT=${batchIntervallistDir} \
	UNIQUE=true \
	SCATTER_COUNT=${BATCHCOUNT}

	echo "AUTOSOMAL DONE"
	#non PAR region
	java -jar -Xmx4g -XX:ParallelGCThreads=4 ${EBROOTPICARD}/picard.jar IntervalListTools \
     	INPUT=${chrXNONPARInterval} \
     	OUTPUT=${batchIntervallistDir} \
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

	for i in $(seq $((${BATCHCOUNT}-2)) ${BATCHCOUNT})
	do
	awk '{
	        if($1 == "Y"){
       	        	print $0
        	}

	}' ${baits}.batch-${i}.bed >> ${baits}.batch-${BATCH_Y}Y.bed
	done

	for i in $(seq $((${BATCHCOUNT}-2)) ${BATCHCOUNT})
	do
		sed '/^Y/ d' ${baits}.batch-${i}.bed > ${baits}.batch-${i}.bed.tmp ; mv ${baits}.batch-${i}.bed.tmp ${baits}.batch-${i}.bed
	done

	for i in $(ls ${MAP}/*batch-*.bed); do cat $i | awk -v var="$i" '{if( $2==$3){print var}}';done > ${MAP}/chompLines.txt	

	while read line
	do
		awk '{if ($2!=$3){ print $0}}' $line > ${line}.tmp
		mv ${line}.tmp $line 
	done<${MAP}/chompLines.txt

	##### Because bed is 0-based and intervallist 1-based, do start minus 1
	for i in $(ls ${baits}.batch*.bed)
	do
		awk '{
        	if ($0 !~ /^@/){
                	minus=($2 - 1);
                	print $1"\t"minus"\t"$3"\t"$4"\t"$5
        	}
        	else
        	        print $0
        	}' $i > ${i}.tmp
		mv ${i}.tmp $i
	done

	echo "batching complete"
	rm -rf ${batchIntervallistDir}/temp_0*
fi
fi #end of if/else loop chr
if [ ! -f ${MAP}/captured.femaleY.bed ]
then
	echo -e 'Y\t1\t2\t+\tFake' > ${MAP}/captured.femaleY.bed
fi

#for f in ${MAP}/*_baits_*; do cp $f ${f/_baits_/_exons_}; done

if [ -f ${baits}.interval_list ]
then
	echo "interval_list already exists, skipping"
else
	perl ${EBROOTNGSMINUTILS}/create_interval_listV4.pl -Ref ${phiXRef}  -Exons ${baits}  -Baits ${baits}
	echo "intervals created"
fi

