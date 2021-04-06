#!/usr/bin/env python
import sys
import math
import csv
from collections import defaultdict

columns = defaultdict(list)

dictAvgCoverage=dict()
dictCount=dict()
dictSD=dict()
dictM10=dict()
dictM20=dict()
dictM30=dict()
dictM50=dict()
dictM100=dict()
#print sys.argv[1]
reader = csv.DictReader(open(sys.argv[1], "rb"), delimiter="\t")
for row in reader:
	gene=row['Gene']
	cov=row['AvgCoverage']
	sd=row['SD']
	moreThan10x=row['moreThan10x']
	moreThan20x=row['moreThan20x']
	moreThan30x=row['moreThan30x']
	moreThan50x=row['moreThan50x']
	moreThan100x=row['moreThan100x']
	if row['Gene'] not in dictAvgCoverage:
		dictAvgCoverage[gene] = 0 
		dictSD[gene] = 0
		dictM10[gene] = 0
		dictM20[gene] = 0
		dictM30[gene] = 0
		dictM50[gene] = 0
		dictM100[gene] = 0
		dictCount[gene] = 0

	dictAvgCoverage[gene] += float(cov)
	dictSD[gene] += float(sd)
	if moreThan10x:
		dictM10[gene] += float(moreThan10x)
	if moreThan20x:
		dictM20[gene] += float(moreThan20x)
	if moreThan30x:
		dictM30[gene] += float(moreThan30x)
	if moreThan50x:
		dictM50[gene] += float(moreThan50x)
	if moreThan100x:
		dictM100[gene] += float(moreThan100x)
	dictCount[gene] += 1


for gene in dictAvgCoverage:
	c=dictAvgCoverage[gene]/dictCount[gene]
	de=dictSD[gene]/dictCount[gene]
	e=dictM10[gene]/dictCount[gene]
	f=dictM20[gene]/dictCount[gene]
	g=dictM30[gene]/dictCount[gene]
	h=dictM50[gene]/dictCount[gene]
	i=dictM100[gene]/dictCount[gene]
	print gene + "\t"+  str(c) + "\t" + str(dictCount[gene]) + "\t"+ str(de) +"\t"+str(e) + "\t"+str(f) + "\t"+ str(g)+ "\t"+str(h) + "\t"+str(i)
