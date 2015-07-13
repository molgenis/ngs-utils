#!/bin/bash


if [ $# -lt 2 ]; then
	echo "Usage: sortVCF in.vcf.gz out.vcf.gz"
	echo "Note1: only takes vcf.gz files (bgzip compressed)"
	echo "Note2: You must have the bgzip util in your PATH. Available from SAMTools tabix package. On Millipede in the /data/gcc/tools/tabix-x/ folder"
	echo "Adapted from VCFTools shell one-liners"
else

	echo `(zcat $1 | head -100 | grep ^#; zcat $1 | grep -v ^# | sort -k1,1d -k2,2n;) | bgzip -c > $2`
	
	

fi

