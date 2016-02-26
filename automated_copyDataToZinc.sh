set -e
set -u
SAMPLESHEETS=$(ssh umcg-rkanninga@gattaca02.gcc.rug.nl "ls /groups/umcg-gaf/scr01/Samplesheets/*.csv" > /groups/umcg-gaf/tmp05/Samplesheets/allSampleSheets.txt)
LOGDIR=/groups/umcg-gaf/tmp05/logs/
echo "Logfiles will be written to /groups/umcg-gaf/tmp05/logs/"
while read line
do
	echo "LINE: $line"
	fileName=$(basename $line)
	fileNameWithoutExt="${fileName%.*}"
	LOGGER=${LOGDIR}/${fileNameWithoutExt}.copyToZinc.logger
echo "1"
	FINISHED="no"
	OLDIFS=$IFS
	IFS=_
	set $fileName
	sequencer=$2
	run=$3
	IFS=$OLDIFS
echo "2"
	isSampleSheetCopied="no"
	isDataCopiedToZinc="no"
	isDataCopiedToPrm="no"
	pipeline="empty"

	## Check if samplesheet is copied
	if [ -f /groups/umcg-gaf/tmp05/Samplesheets/$fileName ]
	then
		isSampleSheetCopied="yes"
	fi

	## Check if data is already copied to tmp05 on zinc-finger
	if [ -f /groups/umcg-gaf/tmp05/rawdata/ngs/$fileName ]
	then
		isRawDataCopiedToZinc="yes"
	fi
echo "3"
	### Check if data is already copied to prm (on calculon)
	if [ -f /groups/umcg-gaf/prm02/rawdata/ngs/$fileName ]
	then
		isRawDataCopiedToPrm="yes"
	fi
echo "4"
	##Check if demultiplexing is finished
	isFINISHED=$(ssh umcg-rkanninga@gattaca02.gcc.rug.nl "if [ -f /groups/umcg-gaf/scr01/logs/$fileName.Demultiplexing.finished ]; then echo 'yes'; else echo 'no';fi")
echo "4.1"
	## Check if rawdata is already on prm
	onPRM=$(ssh umcg-rkanninga@calculon.hpc.rug.nl "if [ -d /groups/umcg-gaf/prm02/rawdata/ngs/${fileName}/ ]; then echo 'yes'; else echo 'no'; fi")

	###Step 1Check if samplesheet already exists
	if [ ! -f /groups/umcg-gaf/tmp05/Samplesheets/$fileName ]
	then
		echo "Copying Samplesheet.." >> ${LOGGER}
		###COPY samplesheet to zinc-finger
		scp umcg-rkanninga@gattaca02.gcc.rug.nl:/groups/umcg-gaf/scr01/Samplesheets/$fileName /groups/umcg-gaf/tmp05/Samplesheets/
		isSampleSheetCopied="yes"
		echo "Samplesheet copied" >> ${LOGGER}
	fi
echo "5"
	if [ "$isFINISHED" == "yes" ] && [ "$isSamplesheetCopied" == "yes" ]
	then
		##Step 2 Check if data is already copied to tmp05 on zinc-finger
		if [ isRawDataCopiedToZinc == "no" ]
		then
			echo "Copying rawdata from gattaca server to zinc-finger" >> ${LOGGER}

			scp -r umcg-rkanninga@gattaca02.gcc.rug.nl:/groups/umcg-gaf/scr01/runs/run_${run}_${sequencers}/results/* /groups/umcg-gaf/tmp05/rawdata/ngs/
			isRawDataCopiedToZinc="yes"
			echo "data copied to zinc-finger" >> ${LOGGER}
		fi
	fi
	### Step 3: Check if data is already copied to prm (on calculon)
echo "6"
	if [ "$isRawDataCopiedToZinc" == "yes" ]
	then
		scp -r /groups/umcg-gaf/tmp05/rawdata/ngs/${fileName}/ umcg-rkanninga@calculon.hpc.rug.nl:/groups/umcg-gaf/prm02/rawdata/ngs/

		### Step 4: Does the pipeline need to run?
		if [ $pipeline == "RNA-Lexogen-reverse" ]
		then
			echo "RNA-Lexogen-reverse" >> ${LOGGER}
		elif [ $pipeline == "RNA-Lexogen" ]
		then
			echo "RNA-Lexogen" >> ${LOGGER}
		elif [ $pipeline == "RNA" ]
		then
			echo "RNA" >> ${LOGGER}
		elif [ $pipeline == "DNA" ]
		then
			echo "DNA" >> ${LOGGER}
		else
			echo "Pipeline is skipped" >> ${LOGGER}
		fi
	fi
echo "7"
done<allSamplesheets.txt
