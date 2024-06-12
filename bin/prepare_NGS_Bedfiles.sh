#!/bin/bash
set -u
set -e
underline=`tput smul`
normal=`tput sgr0`
bold=`tput bold`

function usage () {
echo "
${bold}This script is a small pipeline of 3 seperate scripts creating:
- an interval_list
- per chromosome bed file
- optionally a coverage_per_base bed file and interval_list ${normal}

These scripts are also available to run seperately, please run the scripts below from the /gcc/tools/scripts folder:
To make interval lists          --> create_interval_listV4.pl
To split bed into batches	--> scatter_gather.sh
Make coverage per base file	--> coverage_per_base.sh

${bold}Arguments${normal}

	Required:
	-n|--name		BED name (${bold}SHOULD ALWAYS BE the name 'captured')${normal}
	Optional:

	-c|--coverageperbase	true or false (default: false)
	-d|--data 	What kind of data? Small (targeted) or bigger (chromosome; default)
					Small batchsize of 5 (3 autosomal + 1X + 1Y)
	-e|--extension		extension to the bed file (default: human_g1k_v37)
	-o|--intervalfolder	path to intervalfolder (default: this folder)
	-r|--reference		Which reference file is used (default: /apps/data/1000G/phase1/human_g1k_v37_phiX.dict)
	-t|--tmp		Give tmpfolder location (default: /groups/umcg-gaf/tmp04/tmp)"
}

module load ngs-utils
PARSED_OPTIONS=$(getopt -n "$0"  -o n:o:e:r:c:d:t: --long "name:intervalfolder:extension:reference:coverageperbase:data:tmp:"  -- "$@")

#
# Bad arguments, something has gone wrong with the getopt command.
#
if [ $? -ne 0 ]; then
        usage
        echo "FATAL: Wrong arguments."
        exit 1
fi

eval set -- "$PARSED_OPTIONS"

#
# Now goes through all the options with a case and using shift to analyse 1 argument at a time.
# $1 identifies the first argument, and when we use shift we discard the first argument, so $2 becomes $1 and goes again through the case.
#
while true; do
  case "$1" in
        -n|--name)
                case "$2" in
                "") shift 2 ;;
                *) NAME=$2 ; shift 2 ;;
            esac ;;
	-c|--coverageperbase)
                case "$2" in
                *) COVPERBASE=$2 ; shift 2 ;;
            esac ;;
	-d|--data)
                case "$2" in
                *) DATA=$2 ; shift 2 ;;
            esac ;;
	-o|--INTERVALFOLDER)
                case "$2" in
                *) INTERVALFOLDER=$2 ; shift 2 ;;
            esac ;;
	-e|--extension)
                case "$2" in
                *) EXTENSION=$2 ; shift 2 ;;
            esac ;;
	-r|--reference)
                case "$2" in
                *) REFERENCE=$2 ; shift 2 ;;
            esac ;;
	-t|--tmp)
                case "$2" in
                *) TMP=$2 ; shift 2 ;;
            esac ;;
        --) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
  esac
done

empty=""
#
# Check required options were provided.
if [[ -z "${NAME-}" ]]; then
	usage
	exit 1
fi
if [[ -z "${INTERVALFOLDER-}" ]]; then
	INTERVALFOLDER="."
fi
if [[ -z "${EXTENSION-}" ]]; then
        EXTENSION="human_g1k_v37"
fi
if [[ -z "${REFERENCE-}" ]]; then
        REFERENCE="/apps/data/1000G/phase1/${EXTENSION}"
fi
if [[ -z "${COVPERBASE-}" ]]; then
        COVPERBASE="false"
fi
if [[ -z "${DATA-}" ]]; then
        DATA="chr"
fi
if [[ -z "${TMP-}" ]]; then
	THISDIR=$(pwd)
	TMP=${THISDIR}/TMP/
	mkdir -p ${TMP}
fi

BATCHCOUNT=3
phiXRef="${REFERENCE}_phiX.dict"
phiXExt="${EXTENSION}_phiX"
batchCount_X=2

echo "NAME: ${NAME}"
echo "INTERVALSFOLDER: ${INTERVALFOLDER}"
echo "EXTENSION: ${EXTENSION}"
echo "REFERENCE: ${REFERENCE}"
echo "COVPERBASE: ${COVPERBASE}"
echo "TMPDIR: ${TMP}"

MAP="${INTERVALFOLDER}"

if [[ "${NAME}" == *"baits"*  || "${NAME}" == *"v37"* || "${NAME}" == *"exons"* || "${NAME}" == *".bed"* ]]
then
	echo "No need to put extension (baits, exons, v37 or .bed) in the name, only the name of the bed"
	exit 0
fi

baits="${MAP}/${NAME}"

if [ -f "${baits}.batch-1.bed" ]
then
	echo "splitting in batches skipped"
elif [[ "${NAME}" == *"baits"*  || "${NAME}" == *"v37"* || "${NAME}" == *"exons"* || "${NAME}" == *".bed"* ]]
then
	echo "No need to put extension (baits, exons, v37 or .bed) in the name, only the name of the bed"
	exit 0
fi

module load ngs-utils
module load BEDTools

## sorting bed
sort -V "${baits}.bed" > "${baits}.bed.sorted"
mv "${baits}.bed.sorted" "${baits}.bed"

bedtools merge -i "${baits}.bed" -c 4 -o distinct > "${baits}.merged.bed"
perl -pi -e "s/\r//g" "${baits}.merged.bed"
wc -l  "${baits}.bed"
wc -l  "${baits}.merged.bed"

#
##
### MAKE INTERVALLIST OUT OF BED FILE (0-BASED) and 4th column is strand
##
#

cat "${phiXRef}" > "${baits}.interval_list"
cat "${phiXRef}" > "${baits}.autosomal.interval_list"

awk '{print $1"\t"$2+1"\t"$3"\t+\t"$4}' "${baits}.merged.bed" >> "${baits}.interval_list"
awk '{if ($1 ~ /^[0-9]+$/){print $1"\t"$2+1"\t"$3"\t+\t"$4}}' "${baits}.merged.bed" >> "${baits}.autosomal.interval_list"

#
##
### Make genesOnly file
##
#
if [[ ! -f "${baits}.genesOnly" ]]
then
	awk '{print $4}' "${baits}.merged.bed" > "${baits}.genesOnly"
fi

if [[ "${COVPERBASE}" == "true" ]]
then
	if [ ! -f "${baits}.uniq.per_base.bed" ]
	then
		echo "starting to create_per_base_bed, this may take a while"
		create_per_base_bed.pl -input "${baits}.merged.bed" -output "${NAME}" -outputfolder "${TMP}"
		awk '{print $1"\t"$2"\t"($3+1)"\t"$5}' "${TMP}/${NAME}.per_base.bed" > "${baits}.uniq.per_base.bed"
		wc -l "${TMP}/${NAME}.per_base.bed"
		#

		#sort -V -k1 -k2 -k3 "${TMP}/${NAME}.per_base.intervals" | uniq > "${baits}.uniq.per_base.intervals.tmp"
		#head -n 86 "${baits}.interval_list" > "${baits}.uniq.per_base.interval_list"
		#sort -V "${baits}.uniq.per_base.intervals.tmp" >> "${baits}.uniq.per_base.interval_list"
		#tail -n+87 "${baits}.uniq.per_base.interval_list" |  awk '{print $1"\t"$2"\t"($3+1)"\t"$5}' > "${baits}.uniq.per_base.bed"
	fi
fi


if [ "${DATA}" == "chr" ]
then
	## PER CHROMOSOME
	if [ -f "${baits}.batch-1.bed" ]
	then
		echo "Is this bed file already splitted before? If so, please remove the old ones or do not run this script ;)"
	else
		batchIntervallistDir="${TMP}"

		awk '{ print $0 >> "captured.batch-"$1".bed"}' "${baits}.merged.bed"
	fi
else
	## BATCHING
	if [ -f "${baits}.batch-1.bed" ]
	then
		echo "Is this bed file already splitted before? If so, please remove the old ones or do not run this script ;)"
	else
		module load picard

		batchIntervallistDir="${TMP}"

		#autosomal: split up in 3 batches	
		java -jar -Xmx4g -XX:ParallelGCThreads=4 ${EBROOTPICARD}/picard.jar IntervalListTools \
		INPUT="${baits}.autosomal.interval_list" \
		OUTPUT="${batchIntervallistDir}" \
		UNIQUE=true \
		SCATTER_COUNT='3'

		for i in '1' '2' '3'
		do
			mv "${batchIntervallistDir}/temp_000${i}_of_3/scattered.interval_list" "${baits}.batch-${i}.interval_list"
			grep -v '^@' "${baits}.batch-${i}.interval_list" | awk '{print $1"\t"$2-1"\t"$3"\t"$5}' > "${baits}.batch-${i}.bed"
		done	
		rm "${baits}.batch"*".interval_list"
		echo "AUTOSOMAL DONE"
		
		##chromosomes X and Y
		awk '{ if($1 == "X"){print $0}}' "${baits}.merged.bed" > "captured.batch-X.bed"
		awk '{ if($1 == "Y"){print $0}}' "${baits}.merged.bed" > "captured.batch-Y.bed"
		echo -e 'Y\t1\t2\tFake' > "${MAP}/captured.femaleY.bed"
		echo "batching complete"
		rm -rf "${batchIntervallistDir}/temp_0"*
	fi
	

	if grep phiX174 "${baits}.merged.bed"
	then
		if [ ! -f "${baits}.batch-NC_001422.1.merged.bed" ]
		then
			echo -e 'NC_001422.1\t0\t5386\tphiX174' > "${baits}.batch-NC_001422.1.bed"
		fi
		echo "phiX already inside bed file" 
	else
		echo -e 'NC_001422.1\t0\t5386\tphiX174' >> "${baits}.bed"
		echo -e 'NC_001422.1\t0\t5386\tphiX174' >> "${baits}.merged.bed"
		echo -e 'NC_001422.1\t0\t5386\tphiX174' > "${baits}.batch-NC_001422.1.bed"
	fi
fi
bedtools intersect -a "${baits}.merged.bed" -b '/apps/data/GSAarray/GSA_sorted.bed' | uniq > 'GSA_SNPS.bed'

