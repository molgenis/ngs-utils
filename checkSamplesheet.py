#!/usr/bin/env python

import argparse
import os
import csv
import sys
import re
from collections import defaultdict
from os.path import basename
parser = argparse.ArgumentParser(description='Process commandline opts.')
parser.add_argument("--input")
parser.add_argument("--log")
args = parser.parse_args()

columns = defaultdict(list)
f = open(args.input, 'r') # opens the samplesheet file.
print("INFO: input = " + args.input)
reader = csv.DictReader(f)  # creates the reader object.
inputFileName=(basename(args.input))
inputFileNameBase = re.sub('\..*$', '', inputFileName)

#
# Parse meta-data from the filename. 
#
inputFileNameComponents = inputFileNameBase.split('_')
sequencingStartDate = inputFileNameComponents[0]
sequencer= inputFileNameComponents[1]
run = inputFileNameComponents[2]
flowcell = inputFileNameComponents[3]

if len(inputFileNameComponents) > 4:
	for i in range(4,len(inputFileNameComponents)):
		flowcell+="_"+ str(inputFileNameComponents[i])

w = open(args.log, 'w')
print("INFO: log   = " + args.log)
sanityCheckOk=True
alreadyErrored=False
hasRows = False
listOfErrors=[]

#
# Iterate over the rows of the file.
#
for number, row in enumerate(reader,1):
	hasRows = True
	#
	# Check if the required columns are present.
	#
	for columnName in ('externalSampleID','project','sequencer','sequencingStartDate','flowcell','run','flowcell','lane','seqType','prepKit','capturingKit','barcode','barcodeType'):
		if columnName not in row.keys():
			sanityCheckOk=False
			if not alreadyErrored:
				listOfErrors.extend('ERROR: Required column is missing (or has a trailing space): ' + columnName)
				alreadyErrored=True
		else:
			if row[columnName] == "":
				sanityCheckOk=False
				if columnName in ('capturingKit','barcode','barcodeType'):
					listOfErrors.append('ERROR on line ' + str(number) + ': Variable ' + columnName + ' is empty! Please fill in "None" (to make sure it is not missing).')
				else:
					listOfErrors.append('ERROR on line ' + str(number) + ': Variable ' + columnName + ' is empty!')
	#
	# Check if the data inside the file matches the expected filename.
	#
	if row['sequencer'] != sequencer and 'sequencer' in row.keys():
		sanityCheckOk=False
		listOfErrors.append('ERROR on line ' + str(number) + ': sequencer value in samplesheet (' + row['sequencer'] + ') does not match sequencer in filename (' + sequencer + ').')
	if row['sequencingStartDate'] != sequencingStartDate and 'sequencingStartDate' in row.keys():
		sanityCheckOk=False
		listOfErrors.append('ERROR on line ' + str(number) + ': sequencingStartDate value in samplesheet (' + row['sequencingStartDate'] + ') does not match sequencingStartDate in filename (' + sequencingStartDate + ').')
	if row['run'] != run  and 'run' in row.keys():
		sanityCheckOk=False
		listOfErrors.append('ERROR on line ' + str(number) + ': run value in samplesheet (' + row['run'] + ') does not match run in filename (' + run + ').')
	if row['flowcell'] != flowcell and 'flowcell' in row.keys():
		sanityCheckOk=False
		listOfErrors.append('ERROR on line ' + str(number) + ': flowcell value in samplesheet ' + row['flowcell'] + ' does not match flowcell in filename (' + flowcell + ').')

f.close()

if not hasRows:
	sanityCheckOk=False
	print("File is empty?!")
	listOfErrors.append("File is empty?!")

if sanityCheckOk:
	w.write("OK")
	w.close()
	sys.exit(0)
else:
	print('\n'.join(listOfErrors))
	w.write('\n'.join(listOfErrors))
	w.close()
	sys.exit(1)

