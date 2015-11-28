#this script converts .dose files formatted in beagle 'single dosage' format, value 0..2.00
#in this format the header has 2 columns per individual for fam+ind id
#each data row has only one column per individual
#using plink you can read this via plink --fam d.fam --dosage a.txt list format=1
#the a.txt file contains a list of all .dose files

#######################
#create argument parser
#######################
import optparse
parser = optparse.OptionParser(usage='usage: %prog -s subsfile -d dosagefile -o outputfile\nChoose option -h for extensive help')
parser.add_option('-s', '--subsetFile', metavar='FILE', help='tab delimited file with two columns defining selectedIds and pseudoIds')
parser.add_option('-d', '--doseFile', metavar='FILE', help='.dose file in beagle single dosage format, i.e. header: SNP A1 A2 F1 I1 F2 I2.. and data rs123 A T 1.99 0.11 (so value per 2 headers)')
parser.add_option('-o', '--outFile', metavar='FILE', help='where the output will be written to')
parser.add_option('-t', '--test', help='when this option is set to True then only 1000 lines will be converted',default=False)
(options, args) = parser.parse_args()

print options
if options.subsetFile==None or options.doseFile==None or options.outFile==None:
	parser.print_help()
	exit()

###################################################
#read subsetFile into two lists, of old and new ids
###################################################
import csv
csvReader = csv.reader(open(options.subsetFile, 'rb'), delimiter='\t')
selectedIds = []
pseudoIds = {}
for row in csvReader:
	selectedIds.append(row[1])
	pseudoIds[row[1]]=row[3]

if options.test:
	for el in selectedIds:
		print el+'='+pseudoIds[el]

##########################################################################################
#iterate through dose file, and rename columns in first row remember what indexes to keep.
#expected header= SNP A1 A2 F1 I1 F2 F2 (two columns per individual)
#expected data=rs1234 A T 1.99 0.69 (one column per individual)
##########################################################################################

csvReader = csv.reader(open(options.doseFile, 'rb'), delimiter='\t', skipinitialspace=True,quoting=csv.QUOTE_NONE)
csvWriter = csv.writer(open(options.outFile, 'wb'), delimiter='\t')
count = 0

selectedColHeaders = [0,1,2]
selectedCols = [0,1,2]
noElements = 0
#check if a familySEPid is split by " " or "\t"
#tab = True
for row in csvReader:

	if count == 0:
	        #if first row then find out what columns to keep and rename those with new ids
		print 'filtering data headers'
	
		noElements = len(row)
		idx = 0

		#select cols
		myheader = [row[0],row[1],row[2]]
		for col in row:

			#print 'testing col: '+col+' on idx '+str(idx)

			#[0,1,2] are always included
			if idx > 2:

				print 'business '+str(idx)
                        	#are the family and identifiers seperated by tab (skip=true) or space?
                        	#family = row[idx-1]
                        	#id = row[idx]
                        	#if " " in col:
				#	tab = False
                                id = col.partition(" ")[2]
                                family = col.partition(" ")[0]
                        	
				print 'testing id ['+id+']'

				if id in selectedIds:
					#because of single dose format the data column idx is half of header column idx
					#if tab:
					#	selectedCols.append( 3 + (idx-3)/2)
					#else:
					selectedCols.append(idx)
					
					print 'mapping '+id
					
					#if tab:
					#	selectedColHeaders.append(idx-1)
					selectedColHeaders.append(idx)

					myheader.append(family + " " +pseudoIds[id])
					#myheader.append(pseudoIds[id])

				else:
					print 'skipping '+id
			idx = idx+1
			
			#progress monitoring
			if idx % 1000 == 0:
				print 'filtered header index: index'+str(idx)

		#write the header
                csvWriter.writerow(myheader)

		#debug info
		#if tab:
		#		print 'selected samples: '+str(len(selectedIds))+', available geno samples: '+str((idx-3)/2)+', filtered geno samples: '+str( (len(selectedCols)-3) /2 )
		#else:
		print 'selected samples: '+str(len(selectedIds))+', available geno samples: '+str(idx - 3)+', filtered geno samples: '+str(len(selectedCols) - 3)

	else:		
		#filter the row into 'myvalues' and write to csv
		myvalues = []
 
		for col in selectedCols:
			myvalues.append(row[col])
		csvWriter.writerow(myvalues)

		#debug info
		if count % 1000 == 0:
			print 'converting row '+str(count)

	#for testing purpose you can break early
	if options.test and count >= 100:
		break
	count=count+1

##############################################################################################
#QC: check if the headers in outfile match our selectedIds (and do not contain any previousIds)
##############################################################################################
print 'check of output whether all ids are properly converted'

import os
csvReader = csv.reader(open(options.outFile, 'rb'), delimiter='\t', skipinitialspace=True,quoting=csv.QUOTE_NONE)

foundPseudoIds = []
expectedPseudoIds = pseudoIds.values()
count = 0
for row in csvReader:

	if count == 0:
		#check headers
		idx=0
		for col in row:

			#for each col >2 should be pseudoIds
			if idx >2:
				#if not found we delete outfile and give error
				id = col.partition(" ")[2]
				if id not in expectedPseudoIds:
					#os.remove(options.outFile)
					print 'conversion FAILED: id \''+col+'\' not a pseudoId' 
					exit()

				#remember id so we can count
				foundPseudoIds.append(id)
			else:
				#verify no original ids are in there
				if col in selectedIds:
					os.remove(options.outFile)
                                        print 'conversion FAILED: id \''+col+'\' not a pseudoId'
                                        exit()

			idx=idx+1
					
	if count == 1:
		#check values todo
		print 'have to check values still!'

	if count > 1:
		break

	count=count+1

#check for missing 
if len(foundPseudoIds) != len(expectedPseudoIds):
	for key in pseudoIds.keys():
		if pseudoIds[key] not in foundPseudoIds:
			print 'WARNING: mapping  '+key+'='+pseudoIds[key]+' not in dosage file'

print 'conversion completed'

