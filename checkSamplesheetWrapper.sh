set -eu
for i in $(ls /groups/umcg-gd/scr01/Samplesheets/new/*.csv)
do
	python checkSampleSheet_v2.py --input $i --logfile ${i}.logger

        filename=$(basename ${i})
        check=$(cat $i.logger)
        if [ "${check}" == "OK" ]
        then
		echo "samplesheet is OK, moving to the real samplesheetsdir"
                mv ${i} /groups/umcg-gd/scr01/Samplesheets/
                rm ${i}.logger
        else

		echo "samplesheet ${filename} is not correct, see logger"
        fi
done
