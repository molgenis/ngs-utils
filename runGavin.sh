#!/bin/bash

set -e
set -u

host=$(hostname)
echo "HOST is $host"

if [[ "${host}" == "zinc-finger.gcc.rug.nl" || "${host}" == *"zf-compute"* ]]
then
    	echo "tmpie05"
        tmpName="tmp05"
elif [[ "${host}" == "leucine-zipper.gcc.rug.nl" || "${host}" == *"lc-compute"* ]]
then
    	echo "tmpie06"

        tmpName="tmp06"
fi
GAVINHOME=/groups/umcg-gd/${tmpName}/GavinStandAlone/

if [ $host == "calculon" ]
then
	tmpName=tmp04
	GAVINHOME=/groups/umcg-gdio/${tmpName}/GavinStandAlone/
fi

ml GavinStandAlone/1.0.0

cd ${GAVINHOME}/input


if ls ../input/*.vcf.gz 1> /dev/null 2>&1
then
	for i in $(ls ../input/*.vcf.gz)
	do
		withoutGZ=${i%.*}
		gzip -d ${i}
	done
fi
if ls ../input/*.vcf 1> /dev/null 2>&1
then

	perl -pi -e 's|ID=GC,Number=1,Type=Integer|ID=GC,Number=1,Type=Float|g' ../input/*.vcf
	perl -pi -e 's|MQ0Fraction,Number=1,Type=Integer|MQ0Fraction,Number=1,Type=Float|g' ../input/*.vcf
	perl -pi -e 's|INFO=\<ID=MQ0,Number=1,Type=Integer|INFO=\<ID=MQ0,Number=1,Type=Float|g' ../input/*.vcf
	for i in $(ls ../input/*.vcf) 
	do		
		echo "select only chromosomes 1 till M (no GL regions)"
		egrep "^#|^[123456789XYM]" $i > $i.recoded
		mv $i.recoded $i
		filepath=$(readlink -f $i)
		a=$(basename $i)
		b=${a%%.*}
		cp ${EBROOTGAVINSTANDALONE}/GavinStandAlone.sh ../tmpGavin/GavinStandAlone_${b}.sh
		cd ../tmpGavin/

		perl -pi -e "s|GavinStandAlone_VULIN|GavinStandAlone_${b}|" GavinStandAlone_${b}.sh
		perl -pi -e "s|vcfFile=VULHIERJEVCFFILEIN|vcfFile=${filepath}|" GavinStandAlone_${b}.sh
		echo "running $i"
		sbatch GavinStandAlone_${b}.sh
		cd -
	done
fi

cd $GAVINHOME/output
find . -maxdepth 1 -name '*.vcf' -mtime +6 -exec mv '{}' ouderdan1week/ \;

cd $GAVINHOME/tmp/
find . -maxdepth 1 -name '*.vcf' -mtime +14 -exec rm -rf '{}' \;

cd $GAVINHOME/tmpGavin/
find . -maxdepth 1 -name '*.vcf' -mtime +14 -exec rm -rf '{}' \;
