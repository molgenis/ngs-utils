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
scr01 --> rm -rf rawdata/array/IDAT/{glaasje}
scr01 --> rm -rf generatedscripts/{project}
scr01 --> rm -rf tmp/{project}
scr01 --> rm -rf projects/{project}
scr01 --> rm -rf logs/{project}
Usage:
	$(basename $0) OPTIONS
Options:
	-h   Show this help.

   required:
	-p   projectname ngs
	-f   glaasje numbers, when multiple seperator is ',' NO SPACE inbetween (e.g. 1023123;4532101)
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
if [[ -z "${prmOther:-}" ]]; then prm="${prm}"; else prm="${prmOther}" ; fi
if [[ -z "${gattacaOther:-}" ]]; then gattaca="${gattaca}" ; else gattaca="${gattacaOther}" ; fi


echo "project=${project}"
echo "filePrefix=${filePrefix}"
echo "group=${group}"
echo "gattaca=${gattaca}"
echo "prmCluster=${prmCluster}"
echo "prm=${prm}"
echo "tmpFolder=${tmpFolder}"
echo "scr=${scr}"

IFS=',' read -ra GLAASJES <<< \"${filePrefix}\"

echo "switch user to ${group}-dm to remove data from ${prmCluster} and ${project}"
sudo -u ${group}-dm bash << EOF
ssh ${prmCluster} '
for i in ${GLAASJES[@]} ; do echo \$i ; done
rm -rf /groups/${group}/${prm}/logs/${project}
echo "rm -rf /groups/${group}/${prm}/logs/${project}"
rm -rf /groups/${group}/${prm}/projects/${project}
echo "rm -rf /groups/${group}/${prm}/projects/${project}"
mv -f /groups/${group}/${prm}/Samplesheets/${project}.csv /groups/${group}/${prm}/Samplesheets/archive/ 2&>1 /dev/null
echo "mv -f /groups/${group}/${prm}/Samplesheets/${project}.csv /groups/${group}/${prm}/Samplesheets/archive/"
echo "cleaned up prm"
'
exit
EOF

echo "switch user to ${group}-ateambot to remove data from $(hostname -s) for ${project}"
sudo -u ${group}-ateambot bash -l << EOF
id
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
echo cleaned up tmp on diagnostics cluster
echo now moving to ${gattaca} to clean up ${project}

ssh ${gattaca} "
rm -rf /groups/${group}/${scr}/logs/${project}
echo rm -rf /groups/${group}/${scr}/logs/${project}
rm -rf /groups/${group}/${scr}/rawdata/ngs/${filePrefix}
echo rm -rf /groups/${group}/${scr}/projects/${project}
rm -rf /groups/${group}/${scr}/projects/${project}
rm -rf /groups/${group}/${scr}/tmp/${project}
echo rm -rf /groups/${group}/${scr}/tmp/${project}
mv -f /groups/${group}/${scr}/Samplesheets/archive/${project}.csv /groups/${group}/${scr}/Samplesheets/ 2&>1 /dev/null
echo \"cleaned up scr on gattaca\"
"
exit
EOF

