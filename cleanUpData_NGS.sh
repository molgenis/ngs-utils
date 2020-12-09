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
	exit 0
}

while getopts "p:f:g:p:d:h" opt; 
do
	case $opt in h)showHelp;; p)project="${OPTARG}";; f)filePrefix="${OPTARG}";;  g)group="${OPTARG}";; p)prmOther="${OPTARG}";; d)gattacaOther="${OPTARG}";;
esac 
done

## setting defaults
if [[ $(hostname -s) == "leucine-zipper" ]]
then
	prm="prm06"
	prmCluster="chaperone"
	gattaca="gattaca01.gcc.rug.nl"
	tmpFolder="tmp06"
	scr="scr01"
elif [[ $(hostname -s) == "zinc-finger" ]]
then
	prm="prm05"
	prmCluster="coenzyme"
	gattaca="gattac02.gcc.rug.nl"
	tmpFolder="tmp05"
	scr="scr01"
else
	echo "at this moment we can only run this on leucine-zipper (leu-chap-gat1-prm06) or zinc-finger train (zinc-coe-gat2-prm05)"
fi

if [[ -z "${project:-}" ]]; then showHelp ; echo "project not defined" ; fi
if [[ -z "${filePrefix:-}" ]]; then showHelp ; echo "filePrefix not specified" ; fi
if [[ -z "${group:-}" ]]; then showHelp ; echo "group not specified" ; fi
if [[ -z "${prmOther:-}" ]]; then prm="${prm}" ; echo "prm i"; else prm="${prmOther}" ; fi
if [[ -z "${gattacaOther:-}" ]]; then gattaca="${gattaca}" ; else gattaca="${gattacaOther}" ; fi


echo "project=${project}"
echo "filePrefix=${filePrefix}"
echo "group=${group}"
echo "gattaca=${gattaca}"
echo "prmCluster=${prmCluster}"
echo "prm=${prm}"
echo "tmpFolder=${tmpFolder}"
echo "scr=${scr}"

echo "switch user to ${group}-dm to remove data from ${prmCluster} for ${filePrefix} and ${project}"
sudo -u ${group}-dm bash -l << EOF
ssh ${prmCluster} "
rm -rf /groups/${group}/${prm}/logs/${project}
echo \"rm -rf /groups/${group}/${prm}/logs/${project}\"
rm -rf /groups/${group}/${prm}/logs/${filePrefix}
echo \"rm -rf /groups/${group}/${prm}/logs/${filePrefix}\"
rm -rf /groups/${group}/${prm}/rawdata/ngs/${filePrefix}
echo \"rm -rf /groups/${group}/${prm}/rawdata/ngs/${filePrefix}\"
rm -rf /groups/${group}/${prm}/projects/${project}
echo \"rm -rf /groups/${group}/${prm}/projects/${project}\"
mv -f /groups/${group}/${prm}/Samplesheets/${project}.csv /groups/${group}/${prm}/Samplesheets/archive/
echo \"mv -f /groups/${group}/${prm}/Samplesheets/${project}.csv /groups/${group}/${prm}/Samplesheets/archive/\"
mv -f /groups/${group}/${prm}/Samplesheets/${filePrefix}.csv /groups/${group}/${prm}/Samplesheets/archive/
echo \"mv -f /groups/${group}/${prm}/Samplesheets/${filePrefix}.csv /groups/${group}/${prm}/Samplesheets/archive/\"
echo \"cleaned up prm\"
"
exit
EOF

echo "switch user to ${group}-ateambot to remove data from $(hostname -s) for ${project}"
sudo -u ${group}-ateambot bash -l << EOF
rm -rf /groups/${group}/${tmpFolder}/logs/${project}
echo "rm -rf /groups/${group}/${tmpFolder}/logs/${project}"
rm -rf /groups/${group}/${tmpFolder}/generatedscripts/${project}
echo "rm -rf /groups/${group}/${tmpFolder}/generatedscripts/${project}"
rm -rf /groups/${group}/${tmpFolder}/projects/${project}
echo "rm -rf /groups/${group}/${tmpFolder}/projects/${project}"
rm -rf /groups/${group}/${tmpFolder}/tmp/${project}
echo "rm -rf /groups/${group}/${tmpFolder}/tmp/${project}"
rm -f /groups/${group}/${tmpFolder}/Samplesheets/${project}.csv
echo "rm -f /groups/${group}/${tmpFolder}/Samplesheets/${project}.csv"
echo "cleaned up tmp on diagnostics cluster"
echo "now moving to ${gattaca} to clean up ${filePrefix}"
ssh ${gattaca} "
rm -rf /groups/${group}/${scr}/logs/${filePrefix}
echo \"rm -rf /groups/${group}/${scr}/logs/${filePrefix}\"
rm -rf /groups/${group}/${scr}/rawdata/ngs/${filePrefix}
echo \"rm -rf /groups/${group}/${scr}/rawdata/ngs/${filePrefix}\"
rm -rf /groups/${group}/${scr}/runs/${filePrefix}
echo \"rm -rf /groups/${group}/${scr}/runs/${filePrefix}\"
rm -rf /groups/${group}/${scr}/tmp/${filePrefix}
echo \"rm -rf /groups/${group}/${scr}/tmp/${filePrefix}\"
mv /groups/${group}/${scr}/Samplesheets/archive/${filePrefix}.csv /groups/${group}/${scr}/Samplesheets/
echo \"cleaned up scr on gattaca\"

"
exit
EOF

