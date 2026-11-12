set -eu

year='2026'
if [[ $(hostname) == 'copperfist' ]]
then
	tmp_dir='tmp16'
	hostname_datastaging='cf-porch+copperfist'
	prm_dir='prm06'
elif [[ $(hostname) == 'betabarrel' ]]
then
	tmp_dir='tmp15'
	hostname_datastaging='bb-porch+betabarrel'
	prm_dir='prm05'
elif [[ $(hostname) == 'wingedhelix' ]]
then
	tmp_dir='tmp07'
	hostname_datastaging='wh-porch+wingedhelix'
	prm_dir='prm07'
else
	echo "unknown host"
	exit 1
fi

tmp_root_dir="/groups/umcg-genomescan/${tmp_dir}/"
archive="${tmp_root_dir}/logs/archive_${year}/"
echo -e "#\n##\n###\n Script that will calculate runtimes of genomescan runs"
echo -e "Files that will be used are the rsync-from-195.log from ${archive}"
echo -e "###\n##\n#\n"
rm -f "${archive}/rsyncToPrmTiming/totalTimes.tsv"

count=0
for i in $(ls -d "${archive}/10"*)
do
	projectName=$(head -5 "${i}/rsync-from-195.log" | tail -1  | awk 'BEGIN {FS=" "}{print $5}' | awk 'BEGIN {FS="/"}{print $3}' | awk 'BEGIN {FS="-"}{print $4"-"$5}')
	if [[ ! -f  "${archive}/rsyncToPrmTiming/${projectName}.rsyncProjectRun.finished" ]]
	then
		if [[ ${count} == 0 ]]
		then
			echo "run these commands on the chaperone machine to get also information about the transfer time of copying the data back to the chaperone machine:"
			count=1
		fi

		echo "rsync -rv \"/groups/umcg-gd/${prm_dir}/logs/${projectName}/run01.rsyncProjectRun.finished\" \"${hostname_datastaging}:${archive}/rsyncToPrmTiming/${projectName}.rsyncProjectRun.finished\""
	fi
done

if [[ ${count} != 0 ]]
then
	read -p "Everything copied? (y/N): " -n 1 -r
	echo "" 

	if [[ ! $REPLY =~ ^[Yy]$ ]]; then
	    echo "Aborted."
	    exit 1
	fi
fi

for i in $(ls -d "${archive}/10"*)
do
	gsBatch=$(basename $i)

	numberOfSamples=$(grep 'hard-filtered.vcf.gz.tbi' "${i}/rsync-from-195.log" | wc -l)

	#
	# starttime
	startdate=$(head -1 "${i}/rsync-from-195.log" | awk 'BEGIN {FS=" "}{print $1}')
	starttime=$(head -1 "${i}/rsync-from-195.log" | awk 'BEGIN {FS=" "}{print $2}')
	startEpoch=$(date -d "${startdate} ${starttime}" +%s 2>/dev/null)

	#
	## transfer from transfer server to diagnostics cluster
	enddateT=$(tail -1 "${i}/rsync-from-195.log"  | awk 'BEGIN {FS=" "}{print $1}')
	endtimeT=$(tail -1 "${i}/rsync-from-195.log"  | awk 'BEGIN {FS=" "}{print $2}')
	endEpochT=$(date -d "${enddateT} ${endtimeT}" +%s 2>/dev/null)
	durationSecT=$((endEpochT - startEpoch))
	hoursT=$((durationSecT / 3600))
	minutesT=$(( (durationSecT % 3600) / 60 ))
	secondsT=$((durationSecT % 60))
	durationT="${hoursT}h:${minutesT}m:${secondsT}s"
	durationSecPerSampleT=$((durationSecT / numberOfSamples))
	minutesPerSampleT=$(( (durationSecPerSampleT % 3600) / 60 ))
	secondsPerSampleT=$((durationSecPerSampleT % 60))
	durationPerSampleT="${minutesPerSampleT}m:${secondsPerSampleT}s"
	##
	#

	projectName=$(head -5 "${i}/rsync-from-195.log" | tail -1  | awk 'BEGIN {FS=" "}{print $5}' | awk 'BEGIN {FS="/"}{print $3}' | awk 'BEGIN {FS="-"}{print $4"-"$5}')

	#
	## Copy to prm
	copyToPrmFile="${archive}/rsyncToPrmTiming/${projectName}.rsyncProjectRun.finished"
	if [[ -f "${copyToPrmFile}" ]]
	then
		startdateC=$(head -3 "${copyToPrmFile}" | tail -1  | awk 'BEGIN {FS=" "}{print $1}')
		starttimeC=$(head -3 "${copyToPrmFile}" | tail -1  | awk 'BEGIN {FS=" "}{print $2}')
		enddateC=$(tail -1 "${copyToPrmFile}"  | awk 'BEGIN {FS=" "}{print $1}')
		endtimeC=$(tail -1 "${copyToPrmFile}"  | awk 'BEGIN {FS=" "}{print $2}')
		startEpochC=$(date -d "${startdateC} ${starttimeC}" +%s 2>/dev/null)
		endEpochC=$(date -d "${enddateC} ${endtimeC}" +%s 2>/dev/null)
		durationSecC=$((endEpochC - startEpochC))
		hoursC=$((durationSecC / 3600))
		minutesC=$(( (durationSecC % 3600) / 60 ))
		secondsC=$((durationSecC % 60))
		durationSecPerSampleC=$((durationSecC / numberOfSamples))
		durationC="${hoursC}h:${minutesC}m:${secondsC}s"
		minutesPerSampleC=$(( (durationSecPerSampleC % 3600) / 60 ))
		secondsPerSampleC=$((durationSecPerSampleC % 60))
		durationPerSampleC="${minutesPerSampleC}m:${secondsPerSampleC}s"
	else
		echo "rsyncToPrm file not found for ${projectName}, skipping copy to prm calculation for this sample"
		durationPerSampleC="NA"
		durationC="NA"
	fi
	##
	#

	finishedFile="${archive}/${projectName}/run01.projectDataCopiedToPrm.finished"

	#
	## finish
	endDate=$(date -r "${finishedFile}" "+%Y/%m/%d")
	endTime=$(date -r "${finishedFile}" "+%H:%M:%S")
	endEpoch=$(stat -c %Y "${finishedFile}")

	durationSec=$((endEpoch - startEpoch))
	durationSecPerSample=$((durationSec / numberOfSamples))
	#totaltime
	hours=$((durationSec / 3600))
	minutes=$(( (durationSec % 3600) / 60 ))
	seconds=$((durationSec % 60))

	#total time per sample
	hoursPerSample=$((durationSecPerSample / 3600))
	minutesPerSample=$(( (durationSecPerSample % 3600) / 60 ))
	secondsPerSample=$((durationSecPerSample % 60))
	duration="${hours}h:${minutes}m:${seconds}s"
	durationPerSample="${hoursPerSample}h:${minutesPerSample}m:${secondsPerSample}s"
	echo -e "${gsBatch}\t${projectName}\t${numberOfSamples}\t${startdate} - ${starttime}\t${endDate} - ${endTime}\t${duration}\t${durationPerSample}\t${durationT}\t${durationPerSampleT}\t${durationC}\t${durationPerSampleC}" >> "${archive}/rsyncToPrmTiming/totalTimes.tsv"
	##
	#
done
echo -e "gsBatch\tprojectName\tnumberofsamples\tstartdatetime\tenddatetime\tduration\tdurationpersample\tdurationTransfer\tdurationPerSampleTransfer\tdurationCopy\tdurationPerSampleCopy" > "${archive}/rsyncToPrmTiming/totalTimes.sorted.tsv"
sort -rV -k10 "${archive}/rsyncToPrmTiming/totalTimes.tsv" >> "${archive}/rsyncToPrmTiming/totalTimes.sorted.tsv"
echo "total times file can be found here: ${archive}/rsyncToPrmTiming/totalTimes.sorted.tsv"
echo -e "gsBatch\tprojectName\tnumberofsamples\tstartdatetime\tenddatetime\tduration\tdurationpersample\tdurationTransfer\tdurationPerSampleTransfer\tdurationCopy\tdurationPerSampleCopy" > "${archive}/rsyncToPrmTiming/totalTimes.tsv.tmp2"
cat "${archive}/rsyncToPrmTiming/totalTimes.tsv" >> "${archive}/rsyncToPrmTiming/totalTimes.tsv.tmp2"
mv "${archive}/rsyncToPrmTiming/totalTimes.tsv.tmp2" "${archive}/rsyncToPrmTiming/totalTimes.tsv"
echo "run: column -t -s $'\t' ${archive}/rsyncToPrmTiming/totalTimes.tsv to see runtimes"
