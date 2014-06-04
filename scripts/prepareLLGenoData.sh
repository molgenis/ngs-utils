#!/bin/sh

set -u
set -e

module load jdk

if [[ "$(hostname)" != "scheduler02" ]];
then 
	echo "You should be on scheduler02 to run this script"
fi
if [  ${1-:} ];
then
	echo "
	Required:
        -r|--runid		Put here the runnumber
                        	e.g. OV10-0051

        -d|--data		Give the type of Data (Unimputed,1000G,GONL,HapMap2)

	Optional:
	-p|--pseudodir		Put here the absolute path of the directory of the pseudofile
				pseudo example name: OV10-0051.txt
				default: /gcc/groups/lifelines/home/rsnieders/GenoPseudo
									
	" 		 	
	 
else
 

PARSED_OPTIONS=$(getopt -n "$0"  -o r:d:p: --long "runid:,data:,pseudodir:"  -- "$@")

#Bad arguments, something has gone wrong with the getopt command.
if [ $? -ne 0 ];
then
  exit 1
fi

eval set -- "$PARSED_OPTIONS"

RunId=""
Data=""
PseudoDir="/gcc/groups/lifelines/home/rsnieders/GenoPseudo"
GeneratedScript="/gcc/groups/lifelines/tmp03/lifelines/generatedscripts/${RunId}"
Compute5="/gcc/groups/lifelines/tmp03/lifelines/tools/molgenis_compute5/molgenis-compute-core-0.0.1-SNAPSHOT-20130920/"
SampleSheet="/gcc/groups/lifelines/prm02/samplesheets/"

# Now goes through all the options with a case and using shift to analyse 1 argument at a time.
#$1 identifies the first argument, and when we use shift we discard the first argument, so $2 becomes $1 and goes again through the case.
while true; do
  case "$1" in
	-p|--pseudodir)
        	case "$2" in
		"") shift 2 ;;
                *) PseudoDir=$2 ; shift 2 ;;
            esac ;;
   	-r|--runid)
		case "$2" in
                *) RunId=$2 ; shift 2 ;;
            esac ;;
   	-d|--data)
		case "$2" in 
		*) Data=$2 ; shift 2 ;;
            esac ;; 
    	--) shift ; break ;;
	*) echo "Internal error!" ; exit 1 ;;
  esac
done


echo ${PseudoDir}
echo ${RunId}
echo ${Data}

if [[ "${Data}" == "Unimputed" || "${Data}" == "HapMap2" || "${Data}" == "GONL" || "${Data}" == "1000G" ]];
then 
	echo "Data is okay."
else
       echo "please fill in 'Unimputed', 'HapMap2', 'GONL' or '1000G'" 
fi 


#go to directory to get the Geno Data

cd ${PseudoDir}


#copy the mapping file tmp03
if [ -f ${RunId}.txt ];
	then
	cp ${RunId}.txt /gcc/groups/lifelines/tmp03/lifelines/pseudoFiles 
else
	echo "file does not exist"
fi

#go to samplesheets on permanent storage
cd ${SampleSheet}

#copy an already existing old run mapping file 
cp example_${Data}.csv ${RunId}_${Data}.csv

#rename the old run number in the file
perl -pi -e "s|example|${RunId}|g" ${RunId}_${Data}.csv




#make directory in generatedscripts and load jdk
if [ ! -e ${GeneratedScript} ];
then
	mkdir ${GeneratedScript}
	echo "created folder: ${GeneratedScript}"
	cd ${GeneratedScript}
	mkdir ${Data}
elif [[ -d ${GeneratedScript} && -r ${GeneratedScript} && -w ${GeneratedScript} && -x ${GeneratedScript}  ]]; then
	echo "directory already exists, and it is r/w/x" 
	cd ${GeneratedScript}
	if [ ! -d ${Data} ];
	then
		mkdir ${Data}
	fi
	
else
	echo "error, the ${GeneratedScript} already exists, but you don't have read, write or executable rights"
fi

WORKFLOW=""

if [ ${Data} == "Unimputed" ];
then 
	WORKFLOW="workflowUnimputedData.csv"
else
	WORKFLOW="workflow.csv"
fi

#generate code
sh ${Compute5}/molgenis_compute.sh -b pbs -p ${SampleSheet}/${RunId}_${Data}.csv \
-p ${Compute5}/pipelines/LifeLines_update_genotype_Data/parameters.csv \
-w ${Compute5}/pipelines/LifeLines_update_genotype_Data/${WORKFLOW} \
-rundir ${GeneratedScript}/${Data}/

#go to the generatedscripts folder 
cd ${GeneratedScript}/${Data}/

#now change the queue from ‘PBS’ to ‘devel’ since you are on scheduler02
perl -pi -e 's/#PBS -q gaf/#PBS -q devel/g' *.sh

#be sure that you are the lifelines group
newgrp lifelines

fi
