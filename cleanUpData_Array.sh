#!/bin/bash

set -e
set -u


function showHelp() {
	#
	# Display commandline help on STDOUT.
	#
	cat <<EOH
===============================================================================================================
Script to remove data from prm,scr and diagnostic cluster for ARRAY:
prm --> rm -rf rawdata/array/IDAT/{glaasje}
prm --> rm -rf rawdata/array/GTC/{glaasje}
prm --> rm -rf projects/{project}
prm --> rm -rf logs/{project}
prm --> remove (archive) {project} samplesheet from Samplesheets/
tmp --> rm -rf projects/{project}
tmp --> rm -rf generatedscripts/{project}
tmp --> rm -rf tmp/{project}
tmp --> rm -rf logs/{project}
tmp --> rm {project}.csv samplesheet from Samplesheets/
scr01 --> move {project} samplesheet from Samplesheets/archive to Samplesheets/
scr01 --> rm -rf rawdata/array/GTC/{glaasje}
scr01 --> rm -rf generatedscripts/{project}
scr01 --> rm -rf tmp/{project}
scr01 --> rm -rf projects/{project}
scr01 --> rm -rf logs/{project}

OPTIONAL:
scr01 --> rm -rf rawdata/array/IDAT/{glaasje}

Usage:
	$(basename $0) OPTIONS
Options:
	-h   Show this help.

   required:
	-p   projectname array
	-f   glaasje numbers, when multiple seperator is ',' NO SPACE inbetween (e.g. 1023123,4532101)
	-g   which group
    optional:
	-i   remove IDAT from sourceServer also
	-p   which prm (e.g. prm06) default is based on place of execution of this script (leu-chap-gat1-prm06), (zinc-coe-gat-prm05)
	-d   which gattaca (e.g. gattaca01) default is based on place of execution of this script (leu-chap-gat1-prm06), (zinc-coe-gat2-prm05)
===============================================================================================================
EOH
	exit 0
}

while getopts "p:f:g:p:d:h" opt;
do
	case $opt in h)showHelp;; p)project="${OPTARG}";; f)glaasjeNumbers="${OPTARG}";;  g)group="${OPTARG}";; p)prmOther="${OPTARG}";; d)gattacaOther="${OPTARG}";; i)removeIDAT="-i";;
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
	gattaca="gattaca02.gcc.rug.nl"
	tmpFolder="tmp05"
	scr="scr01"
else
	echo "at this moment we can only run this on leucine-zipper (leu-chap-gat1-prm06) or zinc-finger train (zinc-coe-gat2-prm05)"
fi

if [[ -z "${project:-}" ]]; then showHelp ; echo "project not defined" ; fi
if [[ -z "${glaasjeNumbers:-}" ]]; then showHelp ; echo "glaasjeNumbers not specified" ; fi
if [[ -z "${group:-}" ]]; then showHelp ; echo "group not specified" ; fi
if [[ -z "${prmOther:-}" ]]; then prm="${prm}"; else prm="${prmOther}" ; fi
if [[ -z "${gattacaOther:-}" ]]; then gattaca="${gattaca}" ; else gattaca="${gattacaOther}" ; fi


echo "project=${project}"
echo "glaasjeNumbers=${glaasjeNumbers}"
echo "group=${group}"
echo "gattaca=${gattaca}"
echo "prmCluster=${prmCluster}"
echo "prm=${prm}"
echo "tmpFolder=${tmpFolder}"
echo "scr=${scr}"

IFS=',' read -ra GLAASJES <<< "${glaasjeNumbers}"
if [[ "${group}" == 'umcg-gap' ]]
then
	echo "it is not possible to delete data via this script on prm. commands will be printed only (in case a variable is not expanded and everthing will be deleted by accident)"
	for i in "${GLAASJES[@]}" ; do echo -e "rm -rf /groups/${group}/${prm}/rawdata/array/"{IDAT,GTC}"/"${i}"\n" ; done
	echo -e "
	rm -rvf \"/groups/${group}/${prm}/logs/${project}\"\n\
	rm -rvf \"/groups/${group}/${prm}/projects/${project}\"\n\
	mv -vf \"/groups/${group}/${prm}/Samplesheets/${project}.csv\" \"/groups/${group}/${prm}/Samplesheets/archive/\"
	"
else
	echo "switch user to ${group}-dm to remove data from ${prmCluster} and ${project}"
	echo "[${prmCluster}]"
	sudo -u "${group}-dm" bash << EOT
	ssh "${prmCluster}" "
	for i in \"${GLAASJES[@]}\" ; do rm -rvf \"/groups/${group}/${prm}/rawdata/array/\"{IDAT,GTC}\"/\${i}\" ; done
	echo \"removed ${GLAASJES[@]}\"
	rm -rvf \"/groups/${group}/${prm}/logs/${project}\"
	rm -rvf \"/groups/${group}/${prm}/projects/${project}\"
	mv -vf \"/groups/${group}/${prm}/Samplesheets/${project}.csv\" \"/groups/${group}/${prm}/Samplesheets/archive/\"
	echo \"cleaned up prm" \"
	exit
EOT
fi
echo "switch user to ${group}-ateambot to remove data from $(hostname -s) for ${project}"
sudo -u "${group}-ateambot" bash -l << EOF
rm -rvf "/groups/${group}/${tmpFolder}/logs/${project}"
rm -rvf "/groups/${group}/${tmpFolder}/generatedscripts/${project}"
rm -rvf "/groups/${group}/${tmpFolder}/projects/${project}"
rm -rvf "/groups/${group}/${tmpFolder}/tmp/${project}"
rm -vf "/groups/${group}/${tmpFolder}/Samplesheets/${project}.csv"
echo "cleaned up tmp on diagnostics cluster"
echo "now moving to ${gattaca} to clean up ${project}"

ssh "${gattaca}" "
rm -rvf \"/groups/${group}/${scr}/logs/${project}\"
if [[ -n \"${removeIDAT:-}\" ]]
then
	for i in \"${GLAASJES[@]}\" ; do rm -rvf \"/groups/${group}/${scr}/rawdata/array//\"{IDAT,GTC}\"/\${i}\" ; done
else
	for i in \"${GLAASJES[@]}\" ; do rm -rvf \"/groups/${group}/${scr}/rawdata/array//GTC/\${i}\" ; done
fi
rm -rvf \"/groups/${group}/${scr}/projects/${project}\"
rm -rvf \"/groups/${group}/${scr}/tmp/${project}\"
mv -f \"/groups/${group}/${scr}/Samplesheets/archive/${project}.csv\" \"/groups/${group}/${scr}/Samplesheets/\"
echo \"cleaned up scr on gattaca\"
"
exit
EOF