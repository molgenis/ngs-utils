
set -eu
rm 'totalTimes.tsv'
for i in $(ls -d '/groups/umcg-genomescan/tmp16/logs/archive_2026/10'*)
do
	gsBatch=$(basename $i)

	startdate=$(head -1 "${i}/rsync-from-195.log" | awk 'BEGIN {FS=" "}{print $1}')
	starttime=$(head -1 "${i}/rsync-from-195.log" | awk 'BEGIN {FS=" "}{print $2}')
#	endtransferdate=$(tail -1 "${i}/rsync-from-195.log"  | awk 'BEGIN {FS=" "}{print $1}')
#	endtransfertime=$(tail -1 "${i}/rsync-from-195.log"  | awk 'BEGIN {FS=" "}{print $2}')
	projectName=$(head -5 "${i}/rsync-from-195.log" | tail -1  | awk 'BEGIN {FS=" "}{print $5}' | awk 'BEGIN {FS="/"}{print $3}' | awk 'BEGIN {FS="-"}{print $4"-"$5}')

	startEpoch=$(date -d "${startdate} ${starttime}" +%s 2>/dev/null)
#	endTransferEpoch=$(date -d "${endtransferdate} ${endtransfertime}" +%s 2>/dev/null)


	numberOfSamples=$(grep 'hard-filtered.vcf.gz.tbi' "${i}/rsync-from-195.log" | wc -l)


	finishedFile="/groups/umcg-genomescan/tmp16/logs/archive_2026/${projectName}/run01.projectDataCopiedToPrm.finished"

	endDate=$(date -r "${finishedFile}" "+%Y/%m/%d")
	endTime=$(date -r "${finishedFile}" "+%H:%M:%S")
	endEpoch=$(stat -c %Y "${finishedFile}")
	

	if [ -n "$startEpoch" ] && [ -n "${endEpoch}" ]
	then
		durationSec=$((endEpoch - startEpoch))
		durationSecPerSample=$((durationSec / numberOfSamples))
		
        	hours=$((durationSec / 3600))
        	minutes=$(( (durationSec % 3600) / 60 ))
        	seconds=$((durationSec % 60))
        	hoursPerSample=$((durationSecPerSample / 3600))
        	minutesPerSample=$(( (durationSecPerSample % 3600) / 60 ))
        	secondsPerSample=$((durationSecPerSample % 60))
		duration="${hours}h:${minutes}m:${seconds}s"
		durationPerSample="${hoursPerSample}h:${minutesPerSample}m:${secondsPerSample}s"
	else
        	duration="wrong date format"
	fi
	echo -e "${gsBatch}\t${projectName}\t${numberOfSamples}\t${startdate} - ${starttime}\t${endDate} - ${endTime}\t${duration}\t${durationPerSample}" >> 'totalTimes.tsv'


done
echo -e "gsBatch\tprojectName\tnumberofsamples\tstartdate\tstarttime\tenddate\tendDate\tduration\tdurationpersample" > 'totalTimes.sorted.tsv'
sort -rV -k11 'totalTimes.tsv' >> 'totalTimes.sorted.tsv'
