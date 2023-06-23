#!/bin/bash

set -e
set -u


function showHelp() {
        #
        # Display commandline help on STDOUT.
        #
        cat <<EOH
===============================================================================================================
Script to remove data from prm,scr and diagnostic cluster for ngs
NGS:
prm --> rm -rf rawdata/ngs/{filePrefix}
prm --> rm -rf projects/{project}
prm --> rm -rf logs/{fileprefix}
prm --> rm -rf logs/{project}
prm --> remove (archive) {fileprefix} samplesheet from Samplesheets/
prm --> remove (archive) {project} samplesheet from Samplesheets/
tmp --> rm -rf projects/{project}
tmp --> rm -rf generatedscripts/{project}
tmp --> rm -rf tmp/{project}
tmp --> rm -rf logs/{project}
tmp --> rm {project}.csv samplesheet from Samplesheets/
scr01 --> move {fileprefix} samplesheet from Samplesheets/archive to Samplesheets/
scr01 --> rm -rf generatedscripts/{filePrefix}
scr01 --> rm -rf tmp/{filePrefix}
scr01 --> rm -rf runs/{filePrefix}
scr01 --> rm -rf logs/{filePrefix}
scr01 --> rm -rf rawdata/ngs/{filePrefix}
Usage:
        $(basename $0) OPTIONS
Options:
        -h   Show this help.

   required:
        -p   projectname ngs
        -f   rawdata name ngs (called filePrefix) (e.g. 201101_NB501093_0001_AGFCCVASXB)
        -g   which group
    optional:
        -p   which prm (e.g. prm06) default is based on place of execution of this script (leu-chap-gat1-prm06), (zinc-coe-gat-prm05)
        -d   which gattaca (e.g. gattaca01) default is based on place of execution of this script (leu-chap-gat1-prm06), (zinc-coe-gat2-prm05)
===============================================================================================================
EOH
        trap - EXIT
        exit 0
}

while getopts "p:f:g:p:d:h" opt;
do
	case $opt in h)showHelp;; p)project="${OPTARG}";; f)filePrefix="${OPTARG}";;  g)group="${OPTARG}";; p)prmOther="${OPTARG}";; d)gattacaOther="${OPTARG}";;
esac
done

## setting defaults
if [[ $(hostname -s) == "betabarrel" ]]
then
	prm="prm05"
	prmCluster="bb-chaperone"
	tmpFolder="tmp05"
elif [[ $(hostname -s) == "copperfist" ]]
then
	prm="prm06"
	prmCluster="cf-chaperone"
	tmpFolder="tmp06"
elif [[ $(hostname -s) == "wingedhelix" ]]
then
	prm="prm07"
	prmCluster="wh-chaperone"
	tmpFolder="tmp07"
else
	echo "at this moment we can only run this on betabarrel, copperfist or wingedhelix"
fi

if [[ -z "${project:-}" ]]; then showHelp ; echo "project not defined" ; fi
if [[ -z "${filePrefix:-}" ]]; then showHelp ; echo "filePrefix not specified" ; fi
if [[ -z "${group:-}" ]]; then showHelp ; echo "group not specified" ; fi
if [[ -z "${prmOther:-}" ]]; then prm="${prm}" ; echo "prm i"; else prm="${prmOther}" ; fi

echo "project=${project}"
echo "filePrefix=${filePrefix}"
echo "group=${group}"
echo "prmCluster=${prmCluster}"
echo "prm=${prm}"
echo "tmpFolder=${tmpFolder}"

rm -rvf "/groups/${group}/${tmpFolder}/logs/${project}"
rm -rvf "/groups/${group}/${tmpFolder}/generatedscripts/NGS_DNA/${project}"
rm -rvf "/groups/${group}/${tmpFolder}/projects/NGS_DNA/${project}"
rm -rvf "/groups/${group}/${tmpFolder}/tmp/NGS_DNA/${project}"
if [[ -e "/groups/${group}/${tmpFolder}/Samplesheets/NGS_DNA/${project}.csv" ]]
then
	mv -vf "/groups/${group}/${tmpFolder}/Samplesheets/NGS_DNA/${project}.csv" "/groups/${group}/${tmpFolder}/Samplesheets/NGS_DNA/archive/${project}.csv"
fi

echo "cleaned up tmp on NGS_DNA stuff diagnostics cluster"
echo "now moving to clean up the NGS_Demultiplexing data for ${filePrefix}"
rm -rvf "/groups/${group}/${tmpFolder}/logs/${filePrefix}"
rm -rvf "/groups/${group}/${tmpFolder}/rawdata/ngs/${filePrefix}"
rm -rvf "/groups/${group}/${tmpFolder}/runs/NGS_Demultiplexing/${filePrefix}"
rm -rvf "/groups/${group}/${tmpFolder}/tmp/NGS_Demultiplexing/${filePrefix}"
if [[ -e "/groups/${group}/${tmpFolder}/Samplesheets/NGS_Demultiplexing/archive/${filePrefix}.csv" ]]
then
	mv -vf "/groups/${group}/${tmpFolder}/Samplesheets/NGS_Demultiplexing/archive/${filePrefix}.csv" "/groups/${group}/${tmpFolder}/Samplesheets/NGS_Demultiplexing/"
fi
echo "cleaned up scr on NGS_Demultiplexing stuff"

echo "!it is not possible to delete data via this script on prm. commands will be printed only"
echo -e "
go to ${prmCluster} and execute the following commands:
rm -rvf /groups/${group}/${prm}/logs/${project}\n\
rm -rvf /groups/${group}/${prm}/logs/${filePrefix}\n\
rm -rvf /groups/${group}/${prm}/rawdata/ngs/${filePrefix}\n\
rm -rf /groups/${group}/${prm}/projects/${project}\n\
mv -vf /groups/${group}/${prm}/Samplesheets/${project}.csv /groups/${group}/${prm}/Samplesheets/archive/\n\
mv -vf /groups/${group}/${prm}/Samplesheets/${filePrefix}.csv /groups/${group}/${prm}/Samplesheets/archive/"
