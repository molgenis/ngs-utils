#!/bin/bash
set -eu


### simple script that indexes all the available vcf files on prm in the umcg-gd group 
## first 2 find commands are for searching through the GAVIN folder and the regular vcf files
### third 'find' command is to search for all the available manta diploidSV files
#
# Script can be used as input with the pseudo_bam and pseudo_vcf.sh scripts in the ngs-utils repo

mkdir "${HOME}/indexing"

if [ -e "/groups/umcg-gd/prm06" ]
then 
	echo "yes there is prm"
else 
	echo "no prm available here, please run this on a cluster where prm is available" 
	exit 1
fi
	
echo "starting with gavin vcfs"
find /groups/umcg-gd/prm0*/projects/*/run01/results/variants/GAVIN/ -maxdepth 1 -name '*.vcf.gz' > "${HOME}/indexing/AllVCFs.txt"
echo "starting with no Gavin files"
find /groups/umcg-gd/prm0*/projects/*/run01/results/variants/ -maxdepth 1 -name '*final.vcf.gz' -o -name '*.final.vcf'  > "${HOME}/indexing/AllVCFs_noGAVIN.txt"
echo "starting with manta files"
find /groups/umcg-gd/prm0*/projects/*/run01/results/variants/cnv/ -maxdepth 1 -name '*_diploidSV.vcf.gz' > "${HOME}/indexing/AllVCFs_Manta.txt"
