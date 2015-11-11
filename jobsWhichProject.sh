jobname=$1

whichTmp=$(ls -1 /groups/umcg-gaf | grep tmp)

fgrep $jobname /groups/umcg-gaf/${whichTmp}/tmp/submits/*.txt
