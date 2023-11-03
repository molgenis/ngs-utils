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
        -g   which group
        -p   projectname ngs (required if -f is not set)
        -f   glaasje numbers, when multiple seperator is ',' NO SPACE inbetween (e.g. 1023123,4532101)

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
	case $opt in h)showHelp;; p)project="${OPTARG}";; f)glaasjeNumbers="${OPTARG}";;  g)group="${OPTARG}";; p)prmOther="${OPTARG}";; d)gattacaOther="${OPTARG}";;i)removeIDAT="-i";;
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

if [[ -z "${project:-}" &&  -z "${glaasjeNumbers:-}" ]]
then
	showHelp
	echo "project and filePrefix not defined" 
fi

if [[ -z "${project:-}" ]]; then project="NOT_SET" ; fi
if [[ -z "${glaasjeNumbers:-}" ]]; then glaasjeNumbers="NOT_SET" ; fi
if [[ -z "${group:-}" ]]; then showHelp ; echo "group not specified" ; fi
if [[ -z "${prmOther:-}" ]]; then prm="${prm}" ; echo "prm i"; else prm="${prmOther}" ; fi

echo "project=${project}"
echo "glaasjeNumbers=${glaasjeNumbers}"
echo "group=${group}"
echo "prmCluster=${prmCluster}"
echo "prm=${prm}"
echo "tmpFolder=${tmpFolder}"


if [[ "${glaasjeNumbers}" == "NOT_SET" ]]
then
	echo "glaasjeNumbers not defined" 
else
	echo "now moving to clean up the IDAT/GTC data for ${glaasjeNumbers}"
	IFS=',' read -ra GLAASJES <<< "${glaasjeNumbers}"
	if [[ "${group}" == 'umcg-gap' ]]
	then
		echo "!it is not possible to delete data via this script on prm. commands will be printed only"
		if [[ -n "${removeIDAT:-}" ]]
		then
			for i in "${GLAASJES[@]}" ; do echo -e "rm -rf /groups/${group}/${prm}/rawdata/array/"{IDAT,GTC}"/"${i}"\n" ; done
		else
			for i in "${GLAASJES[@]}" ; do echo -e "rm -rf /groups/${group}/${prm}/rawdata/array/GTC/"${i}"\n" ; done
		fi
	else
		if [[ -n "${removeIDAT:-}" ]]
		then
			for i in "${GLAASJES[@]}" ; do rm -rvf "/groups/${group}/${tmpFolder}/rawdata/array//"{IDAT,GTC}"/${i}" ; done
		else
			for i in "${GLAASJES[@]}" ; do rm -rvf "/groups/${group}/${tmpFolder}/rawdata/array//GTC/${i}" ; done
		fi
	fi
fi

if [[ "${project}" == "NOT_SET" ]]
then
	echo "project not defined" 
else
	rm -rvf "/groups/${group}/${tmpFolder}/logs/${project}"
	rm -rvf "/groups/${group}/${tmpFolder}/generatedscripts/"{AGCT,GAP}"/${project}"
	rm -rvf "/groups/${group}/${tmpFolder}/projects/GAP/${project}"
	rm -rvf "/groups/${group}/${tmpFolder}/runs/AGCT/${project}"
	rm -rvf "/groups/${group}/${tmpFolder}/tmp/"{AGCT,GAP}"/${project}"
	mv -vf "/groups/${group}/${tmpFolder}/Samplesheets/GAP/${project}.csv" "/groups/${group}/${tmpFolder}/Samplesheets/archive/${project}.csv"
	
	echo "!it is not possible to delete data via this script on prm. commands will be printed only"
	echo -e "
	go to ${prmCluster} and execute the following commands:
	rm -rvf /groups/${group}/${prm}/logs/${project}\n\
	rm -rvf /groups/${group}/${prm}/projects/${project}\n\
	mv -vf /groups/${group}/${prm}/Samplesheets/${filePrefix}.csv /groups/${group}/${prm}/Samplesheets/archive/"
fi

