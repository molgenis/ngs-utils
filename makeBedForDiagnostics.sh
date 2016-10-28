set -e
set -u

thisDir=$(pwd)

if [ $thisDir != "/apps/data/Agilent" ]
then
	echo "you should be in /apps/data/Agilent to run this"
	exit 1

fi

fileName=$1
if [ -z  ${1+x} ]
then
	echo "expecting 2 arguments (filename and newFile name)"
	echo "e.g. sh makeBedForDiagnostics.sh PCS_3004471_+en-20_target_v2.BED PCS_v4"
	exit 1
fi


if [[ $fileName != *target* ]]
then
	echo "expecting target file, not bp.. If it is the target file, than please put in the name"
	exit 1
fi	

newName=$2

if [ -d /apps/data/Agilent/$newName ]
then
	echo "/apps/data/Agilent/$newName already exists"
	exit 1
elif [ -d /apps/data/UMCG/Diagnostics/$newName ]
then
	echo "/apps/data/UMCG/Diagnostics/$newName already exists"
	exit 1
fi

mkdir -p ${newName}/human_g1k_v37/
echo "created ${newName}/human_g1k_v37/"
cp ${fileName} ${newName}/
echo "copied ${fileName} ${newName}/"
## navigate to folder 
cd ${newName}

withoutExtension=${fileName%.*}

awk '{ 
	if ($0 !~ /^@/){
		plus=($3 +1); print $1"\t"$2"\t"plus"\t"$4 
	}
}' ${fileName} > ${withoutExtension}_StopPlusOne.bed
echo "created ${withoutExtension}_StopPlusOne.bed"

cp ${withoutExtension}_StopPlusOne.bed human_g1k_v37/captured.bed
echo "copied ${withoutExtension}_StopPlusOne.bed to human_g1k_v37/captured.bed"
module load ngs-utils

cd human_g1k_v37/

## Run the prepare step
sh ${EBROOTNGSMINUTILS}/prepare_NGS_Bedfiles.sh -n captured -c true -d targeted

## 
cd $thisDir
echo "copied ${newName} to /apps/data/UMCG/Diagnostics/"
cp -r ${newName} /apps/data/UMCG/Diagnostics/

cd /apps/data/UMCG/Diagnostics/${newName}/human_g1k_v37/
echo "renaming captured into ${newName}"
rename captured ${newName} captured.*

echo "FINISHED"

