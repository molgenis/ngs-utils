#!/bin/bash

set -e
set -u


function showHelp() {
        #
        # Display commandline help on STDOUT.
        #
        cat <<EOH
===============================================================================================================
Script to pseudo anonimize vcf files. 

There is an option to search the database for samples, or you can give a folder containing all the data to be pseudo anonimized 
Usage:
        $(basename $0) OPTIONS
Options:
        -h   Show this help.

   required:
        -s   search database (file containing DNA numbers and mapping in second column, tab seperated) (not in combination with -m)
        -i   input folder containing vcf files (not in combination with -s) 
        -m   mapping file (in combination with -i)
        -g   which group (default = umcg-gd)

		
    optional:
        -e   manta data (default disabled) 
        -p   which prm (e.g. prm06) default is both prm05 and prm06




===============================================================================================================
EOH
        exit 0
}

while getopts "i:g:mp:h" opt;
do
        case $opt in h)showHelp;; i)input="${OPTARG}";; g)group="${OPTARG}";; s)search="${OPTARG}";; m)mapping="${OPTARG}";;  e)manta;; p)prm="${OPTARG}";;
		esac
done


if [[ -z "${input:-}" && -z "${search:-}" ]] 
then  
	showHelp 
	echo "neither input as search option is selected, one of them should be chosen" 
elif [[ -n "${input:-}" && -n "${search:-}" ]]
then
	showHelp 
	echo "input and search option is selected, one of them should be chosen" 
elif [ -z "${input:-}" ]
then
	##search
	echo "search=${search}"
	
elif [ -z "${search:-}" ]
	##input
	echo "hier komt een pad naar de files ${input}"
fi

exit 1
if if [[ -z "${group:-}" ]]; then showHelp ; echo "group not specified" ; fi
if [[ -z "${mapping:-}" ]]; then showHelp ; echo "input is not defined" ; fi
if [[ -z "${prm:-}" ]]; then prm="${prm}"; else prm="${prmOther}" ; fi


echo "project=${project}"
echo "search=${filePrefix}"
echo "group=${group}"
echo "prm=${prm}"


