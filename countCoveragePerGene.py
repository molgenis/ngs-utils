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
dictU10=dict()
dictU20=dict()
dictU50=dict()
dictU100=dict()
#print sys.argv[1]
reader = csv.DictReader(open(sys.argv[1], "rb"), delimiter="\t")
for row in reader:
	gene=row['Gene']
	cov=row['AvgCoverage']
	median=row['Median']
	sd=row['SD']
	u10=row['u10']
	u20=row['u20']
	u50=row['u50']
	u100=row['u100']
	if row['Gene'] not in dictAvgCoverage:
		dictAvgCoverage[gene] = 0 
		dictMedian[gene] = 0
		dictSD[gene] = 0
		dictU10[gene] = 0
		dictU20[gene] = 0
		dictU50[gene] = 0
		dictU100[gene] = 0
		dictCount[gene] = 0

	dictAvgCoverage[gene] += float(cov)
	dictSD[gene] += float(sd)
	dictMedian[gene] += float(median)
	if u10:
		dictU10[gene] += float(u10)
	if u20:
		dictU20[gene] += float(u20)		
	if u50:
		dictU50[gene] += float(u50)
	if u100:
		dictU100[gene] += float(u100)
	dictCount[gene] += 1


for gene in dictAvgCoverage:
	c=dictAvgCoverage[gene]/dictCount[gene]
	de=dictSD[gene]/dictCount[gene]
	d=dictMedian[gene]/dictCount[gene]
	e=dictU10[gene]/dictCount[gene]
	f=dictU20[gene]/dictCount[gene]
	g=dictU50[gene]/dictCount[gene]
	h=dictU100[gene]/dictCount[gene]
	print gene + "\t"+  str(c) + "\t" + str(dictCount[gene]) + "\t"+ str(d) +"\t"+ str(de) +"\t"+str(e) + "\t"+str(f) + "\t"+ str(g)+ "\t"+str(h) 
