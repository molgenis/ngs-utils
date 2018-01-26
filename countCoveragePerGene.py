#!/usr/bin/env python
import sys
import math
import csv
from collections import defaultdict

columns = defaultdict(list)

dictAvgCoverage=dict()
dictCount=dict()
dictMedian=dict()
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
	median=row['Median']
	sd=row['SD']
	moreThan10=row['moreThan10']
	moreThan20=row['moreThan20']
	moreThan30=row['moreThan30']
	moreThan50=row['moreThan50']
	moreThan100=row['moreThan100']
	if row['Gene'] not in dictAvgCoverage:
		dictAvgCoverage[gene] = 0 
		dictMedian[gene] = 0
		dictSD[gene] = 0
		dictM10[gene] = 0
		dictM20[gene] = 0
		dictM30[gene] = 0
		dictM50[gene] = 0
		dictM100[gene] = 0
		dictCount[gene] = 0

	dictAvgCoverage[gene] += float(cov)
	dictSD[gene] += float(sd)
	dictMedian[gene] += float(median)
	if moreThan10:
		dictM10[gene] += float(moreThan10)
	if moreThan20:
		dictM20[gene] += float(moreThan20)
	if moreThan30:
		dictM30[gene] += float(moreThan30)
	if moreThan50:
		dictM50[gene] += float(moreThan50)
	if moreThan100:
		dictM100[gene] += float(moreThan100)
	dictCount[gene] += 1


for gene in dictAvgCoverage:
	c=dictAvgCoverage[gene]/dictCount[gene]
	de=dictSD[gene]/dictCount[gene]
	d=dictMedian[gene]/dictCount[gene]
	e=dictM10[gene]/dictCount[gene]
	f=dictM20[gene]/dictCount[gene]
	g=dictM30[gene]/dictCount[gene]
	h=dictM50[gene]/dictCount[gene]
	i=dictM100[gene]/dictCount[gene]
	print gene + "\t"+  str(c) + "\t" + str(dictCount[gene]) + "\t"+ str(d) +"\t"+ str(de) +"\t"+str(e) + "\t"+str(f) + "\t"+ str(g)+ "\t"+str(h) + "\t"+str(i)
