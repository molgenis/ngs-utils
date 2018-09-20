#!/usr/bin/env python

import argparse
import os
import csv
import sys
from collections import defaultdict
from os.path import basename
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument("--input")
parser.add_argument("--logfile")
args = parser.parse_args()

columns = defaultdict(list)
f = open(args.input, 'r') # opens the csv file
print("inputfile:" + args.input)
reader = csv.DictReader(f)  # creates the reader object
sampleName=(basename(os.path.splitext(args.input)[0]))

#
# Parse meta-data from the filename. 
#
splittedFileName = sampleName.split('_')
sequencestartdate = splittedFileName[0]
sequencer= splittedFileName[1]
runid = splittedFileName[2]
flowcell = splittedFileName[3]

if len(splittedFileName) > 4:
	for i in range(4,len(splittedFileName)):
		flowcell+="_"+ str(splittedFileName[i])

w = open(args.logfile, 'w')
print("logfile:" + args.logfile)
stopRun="false"
alreadyErrored="false"
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
			if alreadyErrored == "false":
				listOfErrors.extend("One required column is missing (or has a trailing space): " + columnName)
				print("One required column is missing (or has a trailing space): " + columnName)
				alreadyErrored="true"
		else:
			if row[columnName] == "":
				if columnName in ('capturingKit','barcode','barcodeType'):
					if alreadyErrored == "false":
						listOfErrors.append("The variable " + sleutel + " on line " + str(number) +  " is empty! Please fill in None (this to be sure that it is not missing)")
						stopRun="true"
						alreadyErrored="true"
				else:
					if alreadyErrored == "false":
						listOfErrors.append("The variable " + columnName + " on line " + str(number) +  " is empty!")
						stopRun="true"
						alreadyErrored="true"
	#
	# Check if the data inside the file matches the expected filename.
	#
	if row['sequencer'] != sequencer and 'sequencer' in row.keys():
		stopRun="true"
		listOfErrors.append("the sequencer in the samplesheet is not matching the sequencer in the filename on line: " + str(number + 1))
		print("the sequencer in the samplesheet is not matching the sequencer in the filename on line: " + str(number + 1))
	if row['sequencingStartDate'] != sequencestartdate and 'sequencingStartDate' in row.keys():
		stopRun="true"
		listOfErrors.append("the sequencingStartDate in the samplesheet is not matching the sequencingStartDate in the filename on line: " + str(number + 1))
		print("the sequencingStartDate in the samplesheet is not matching the sequencingStartDate in the filename on line: " + str(number + 1))
	if row['run'] != runid  and 'run' in row.keys():
		stopRun="true"
		listOfErrors.append("the run in the samplesheet is not matching the run in the filename on line: " + str(number + 1))
		print("the run in the samplesheet is not matching the run in the filename on line: " + str(number + 1))
	if row['flowcell'] != flowcell and 'flowcell' in row.keys():
		stopRun="true"
		listOfErrors.append("the flowcell in the samplesheet is not matching the flowcell in the filename on line: " + str(number + 1))
		print("the flowcell in the samplesheet is not matching the flowcell in the filename on line: " + str(number + 1))


if not hasRows:
	print("The complete file is empty?! in ")
	listOfErrors.append("The complete file is empty?!")

if stopRun == "true":
	w.write('\n'.join(listOfErrors))
	sys.exit(1)
else:
	w.write("OK")

w.close()
f.close()
