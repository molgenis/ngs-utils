#!/bin/bash

##
# MD5SumCheck.sh
##
# WBKoetsier 20120815
##
# this script performs md5sum checks using the input md5 files and their originals
# does not follow symlinks!
# input: dir(s) and/or file(s), md5s or other
# output: stdout
##

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ $# == 0 ]; then
	echo "Please provide input files or directories."
	exit 0
fi

# enable process substitution
set +o posix

# functions must precede calls

## input file checks: is it an md5? Do both the md5 and its original exist?
# after checking, passes the file to CHECK_MD5_SUM
# input: one existing file
function CHECK_AND_PASS_FILE () {
	infile=$1

	cd $(dirname $infile)

	# if infile is md5, check if original exists
	if [[ ${infile/*./} == "md5" ]]; then
		# check if original exists
		# ${infile%.*} -> strip extension
		if [[ -e "${infile%.*}" ]] || [[ -e "${infile%.*}.gz" ]]; then
			CHECK_MD5_SUM $infile
		else
			echo "No matching file found for $md5file"
		fi

	# infile is not md5, check if matching md5 file exists
	else
		if [[ -e "$infile.md5" ]]; then
			CHECK_MD5_SUM "$infile.md5"

		elif [[ -e "${infile%.*}.md5" ]]; then
                        CHECK_MD5_SUM "${infile%.*}.md5"

		else
                        echo "No matching md5 file found for $infile"

                fi

	fi
}

## basically performs md5sum -c file.md5
# takes into account fq.gz (zcat)
# input: existing md5 file of which original also exists
function CHECK_MD5_SUM () {

	md5file=$1

	# fq.gz or not?
	if [[ $md5file == *.fq.* ]]; then
		integrity=$(zcat ${md5file%.*}.gz | md5sum -c <(sed 's/[[:space:]][^[:space:]].*$/ -/' $md5file))
		echo "Checking integrity of unzipped ${md5file%.*}.gz:$integrity"
		
#		if [[ $integrity == *FAILED* ]]; then
#			touch ${md5file%.*}.gz.md5sumcheck.FAILED
#		fi
	else
		integrity=$(md5sum -c $md5file)
		echo "Checking integrity of ${md5file%.*}:$integrity"
		
#		if [[ $integrity == *FAILED* ]]; then
#			touch ${md5file%.*}.md5sumcheck.FAILED
#		fi
	fi

}

## the function calls

for i in $@; do

	# get canonical path to be safe but do not dereference symlinks
	item=$(cd $(dirname "$i"); pwd -P)/$(basename "$i")

	if [ ! -e $item ]; then
		# wrong input.
		echo "$item: no such file or directory, skipping."
		
	elif [ -L $item ]; then
		# it's a symlink.
		echo "$item is a symlink to $(readlink -f $item), skipping."
		
	elif [ -f $item ]; then
		# it's a file.
		CHECK_AND_PASS_FILE $item
	
	elif [ -d $item ]; then
		# it's a directory.
		for f in $(find $item | grep .md5); do
			if [ ! -L $f ]; then
				CHECK_AND_PASS_FILE $f
			else
				echo "$f is a symlink to $(readlink -f $f), skipping."
			fi
		done
	  
	else
		# que?
		echo "$item is not a regular file or directory. Please check. Skipping."

	fi
done

# return to starting point.
cd $DIR
