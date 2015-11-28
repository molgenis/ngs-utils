#this script converts the PCA.txt 

#######################
#create argument parser
#######################
import optparse
parser = optparse.OptionParser(usage='usage: %prog -s subsfile -d pcafile -o outputfile\nChoose option -h for extensive help')
parser.add_option('-s', '--subsetFile', metavar='FILE', help='tab delimited file with two columns defining selectedIds and pseudoIds')
parser.add_option('-d', '--pcaFile', metavar='FILE', help='the pca file)')
parser.add_option('-o', '--outFile', metavar='FILE', help='where the output will be written to')
(options, args) = parser.parse_args()

print options
if options.subsetFile==None or options.pcaFile==None or options.outFile==None:
	parser.print_help()
	exit()

###################################################
#read subsetFile into two lists, of old and new ids
###################################################
import csv
csvReader = csv.reader(open(options.subsetFile, 'rU'), delimiter='\t')
selectedIds = []
pseudoIds = {}
for row in csvReader:
	selectedIds.append(row[1])
	pseudoIds[row[1]]=row[3]

##########################################################################################
#iterate through pca file
#keep line if the id is in the selectedIds
#then rename to pseudoId
##########################################################################################

csvReader = csv.reader(open(options.pcaFile, 'rU'), delimiter='\t')
f = open(options.outFile, 'wb')
csvWriter = csv.writer(f, delimiter='\t')
count = 0

for row in csvReader:
	if count == 0:
		csvWriter.writerow(row)
	else:
		if row[0] in selectedIds:
			#filter the row into 'myvalues' and write to csv
			myvalues = []
			for col in row:
				myvalues.append(col)
			myvalues[0] = pseudoIds[row[0]]
			csvWriter.writerow(myvalues)

			#debug info
			if count % 1000 == 0:
				print 'converting row '+str(count)
	count=count+1

f.close()

##############################################################################################
#QC: check if the rows[0] in outfile match our selectedIds (and do not contain any previousIds)
##############################################################################################
print 'check of output whether all ids are properly converted'

import os
csvReader = csv.reader(open(options.outFile, 'rU'), delimiter='\t')

foundPseudoIds = []
expectedPseudoIds = pseudoIds.values()
count = 0
for row in csvReader:
	if count > 0:
		#check for illegal ids
		if row[0] not in expectedPseudoIds:
			os.remove(options.outFile)
			print 'conversion FAILED: id \''+row[0]+'\' not a pseudoId'
			exit()
		if row[0] in selectedIds:
			os.remove(options.outFile)
			print 'conversion FAILED: id \''+row[0]+'\' not a pseudoId'
			exit()

		#remember id so we can count
		foundPseudoIds.append(row[0])

		#debug info
		if count % 1000 == 0:
			print 'checked row '+str(count)

	count=count+1

print str(len(foundPseudoIds)) +' versus '+str(len(expectedPseudoIds))

#check for missing 
f = open(options.outFile+'.missing', 'w')
if len(foundPseudoIds) != len(expectedPseudoIds):
	for key in pseudoIds:
		if pseudoIds[key] not in foundPseudoIds:
			f.write('WARNING: mapping  '+key+'='+pseudoIds[key]+' not in pca file\n')
			print 'WARNING: mapping  '+key+'='+pseudoIds[key]+' not in pca file'

print 'conversion completed'

