set -eu
manList=$1

ml ngs-utils
EBROOTNGSMINUTILS=$HOME/github/ngs-utils
if ls -d /apps/data/UMCG/ManVarList*  1> /dev/null 2>&1
then
	version=$(ls -d /apps/data/UMCG/ManVarList* | tail -1 | awk 'BEGIN {FS="_v"}{print $2+1}')
else
	version=1
fi
## UMCG-MVL_VKGLconsensusMVL-2017-11-27.txt
awk '{print $1,$2}' "${manList}" \
| awk -F ':' '{print $1"\t"$2"\t"$3"\t"$4}' | awk '{if ($2 ~ /-/){print $1"\t"$2"\t"$3}else {print $1"\t"$2"-"$2"\t"$3}}' | awk -F '-' '{print $1"\t"$2}' | awk '{print $1"\t"($2-1)"\t"$3"\t"$4}' |  sed '1d' > /apps/data/UMCG/ManVarList_v${version}_target.bed
cd /apps/data/UMCG/

sh ${EBROOTNGSMINUTILS}/makeBedForDiagnostics.sh -b ManVarList_v${version}_target.bed -n ManVarList_v${version} -w /apps/data/UMCG/

cd -

