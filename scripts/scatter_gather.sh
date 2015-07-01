set -e 
set -u


INTERVALS=/gcc/resources/b37/intervals
NAMEBED=Agilent_SureSelect_Human_All_Exon_V5_S04380110_plus_Custom_SCA_0655301

MAP=${INTERVALS}/${NAMEBED}/


if [ ! -d ${MAP} ]
then
	mkdir -p ${MAP}/
fi

EXTENSION=baits_b37_human_g1k_v37
BEDFILE=${INTERVALS}/${NAMEBED}_${EXTENSION}

module load picard-tools/1.102

BATCHCOUNT=45
batchIntervallistDir=${MAP}

chrXNONPARInterval=${MAP}/${NAMEBED}_${EXTENSION}.chrX.nonpar.interval_list
chrXPARBed=${MAP}/${NAMEBED}_${EXTENSION}.chrX.par.bed
AllWithoutchrXInterval=${MAP}/${NAMEBED}_${EXTENSION}.withoutChrX.interval_list

head -85  ${INTERVALS}/GEO_Panel_flanked20bp_target_v1_sorted_baits_b37_human_g1k_v37.interval_list > ${AllWithoutchrXInterval}
head -85  ${INTERVALS}/GEO_Panel_flanked20bp_target_v1_sorted_baits_b37_human_g1k_v37.interval_list > ${chrXNONPARInterval}

cat ${BEDFILE}.withoutChrX.bed >> ${AllWithoutchrXInterval}

awk '{
if ($1 == "X"){
	if (($2 >= 60001  && $3 <= 2699520 ) || ($2 >= 154931044 && $3 <= 155260560 )){
        	        print $0 >> "'${chrXPARBed}'"
       	} else {
                print $0 >> "'${chrXNONPARInterval}'"
       	}
}
}' ${BEDFILE}.bed

if [ -f ${chrXPARBed} ]
then	
	cat ${chrXPARBed} >> ${AllWithoutchrXInterval}
fi

#autosomal
java -jar  -Xmx4g -XX:ParallelGCThreads=4 $PICARD_HOME/IntervalListTools.jar \
INPUT=${AllWithoutchrXInterval} \
OUTPUT=${batchIntervallistDir} \
PADDING=150 \
UNIQUE=true \
SCATTER_COUNT=$BATCHCOUNT \
COMMENT="2Added padding of 150bp and merge overlapping and adjacent intervals to create a list of unique intervals PADDING=150 UNIQUE=true"

#non PAR region
java -jar  -Xmx4g -XX:ParallelGCThreads=4 $PICARD_HOME/IntervalListTools.jar \
     INPUT=${chrXNONPARInterval} \
     OUTPUT=${batchIntervallistDir} \
     PADDING=150 \
     UNIQUE=true \
     SCATTER_COUNT=4 \
     COMMENT="3Added padding of 150bp and merge overlapping and adjacent intervals to create a list of unique intervals PADDING=150 UNIQUE=true"

mv  ${batchIntervallistDir}/temp_0001_of_4 ${batchIntervallistDir}/temp_00$((${BATCHCOUNT}+1))_of_${BATCHCOUNT}
mv  ${batchIntervallistDir}/temp_0002_of_4 ${batchIntervallistDir}/temp_00$((${BATCHCOUNT}+2))_of_${BATCHCOUNT}
mv  ${batchIntervallistDir}/temp_0003_of_4 ${batchIntervallistDir}/temp_00$((${BATCHCOUNT}+3))_of_${BATCHCOUNT}
mv  ${batchIntervallistDir}/temp_0004_of_4 ${batchIntervallistDir}/temp_00$((${BATCHCOUNT}+4))_of_${BATCHCOUNT}

BATCH_ALL=$((${BATCHCOUNT}+4))
BATCH_X1=$((${BATCHCOUNT}+1))
BATCH_X2=$((${BATCHCOUNT}+2))
BATCH_X3=$((${BATCHCOUNT}+3))
BATCH_X4=$((${BATCHCOUNT}+4))
BATCH_Y=$((${BATCHCOUNT}+5))


head -85  ${INTERVALS}/GEO_Panel_flanked20bp_target_v1_sorted_baits_b37_human_g1k_v37.interval_list > ${batchIntervallistDir}/${NAMEBED}_${EXTENSION}.batch-${BATCH_Y}.interval_list

for i in $(seq 1 ${BATCH_ALL})
do
if [[ ${i} -lt 10 ]]
then
                mv ${batchIntervallistDir}/temp_000${i}_of_${BATCHCOUNT}/scattered.intervals ${batchIntervallistDir}/${NAMEBED}_${EXTENSION}.batch-${i}.interval_list
                tail -n+88 ${batchIntervallistDir}/${NAMEBED}_${EXTENSION}.batch-${i}.interval_list > ${batchIntervallistDir}/${NAMEBED}_${EXTENSION}.batch-${i}.bed

        elif [ $i == ${BATCH_X1} ] || [ $i == ${BATCH_X2} ] || [ $i == ${BATCH_X3} ] || [ $i == ${BATCH_X4} ]
        then
                        mv ${batchIntervallistDir}/temp_00${i}_of_${BATCHCOUNT}/scattered.intervals ${batchIntervallistDir}/${NAMEBED}_${EXTENSION}.batch-${i}.interval_list
                        tail -n+88 ${batchIntervallistDir}/${NAMEBED}_${EXTENSION}.batch-${i}.interval_list > ${batchIntervallistDir}/${NAMEBED}_${EXTENSION}.batch-${i}X.bed
        elif [ $i == ${BATCH_Y} ]
	then
                        mv ${batchIntervallistDir}/temp_00${i}_of_${BATCHCOUNT}/scattered.intervals ${batchIntervallistDir}/${NAMEBED}_${EXTENSION}.batch-${i}.interval_list
                tail -n+88 ${batchIntervallistDir}/${NAMEBED}_${EXTENSION}.batch-${i}.interval_list > ${batchIntervallistDir}/${NAMEBED}_${EXTENSION}.batch-${i}Y.bed
        else
                        mv ${batchIntervallistDir}/temp_00${i}_of_${BATCHCOUNT}/scattered.intervals ${batchIntervallistDir}/${NAMEBED}_${EXTENSION}.batch-${i}.interval_list
                        tail -n+88 ${batchIntervallistDir}/${NAMEBED}_${EXTENSION}.batch-${i}.interval_list > ${batchIntervallistDir}/${NAMEBED}_${EXTENSION}.batch-${i}.bed

        fi
done

rm ${batchIntervallistDir}/${NAMEBED}_${EXTENSION}.batch*.interval_list

for i in $(seq $(($BATCHCOUNT-10)) ${BATCHCOUNT})
do
awk '{ 
	if($1 == "Y"){
	        print $0
	} 

}' ${batchIntervallistDir}/${NAMEBED}_${EXTENSION}.batch-${i}.bed >> ${batchIntervallistDir}/${NAMEBED}_${EXTENSION}.batch-${BATCH_Y}.bed
done

echo "batching complete"
rm -rf ${batchIntervallistDir}/temp_00*

