#!/usr/bin/env python
import sys
import csv
from collections import defaultdict

columns = defaultdict(list)

dictAvgCoverage=dict()
dictCount=dict()
dictMedian=dict()
dictU10=dict()
dictU20=dict()
dictU50=dict()
dictU100=dict()
reader = csv.DictReader(open("/groups/umcg-gaf/tmp04/umcg-rkanninga/CoverageOverview.txt", "rb"), delimiter="\t")
for row in reader:
	gene=row['Gene']
	cov=row['AvgCoverage']
	median=row['Median']
	u10=row['u10']
	u20=row['u20']
	u50=row['u50']
	u100=row['u100']

	if row['Gene'] not in dictAvgCoverage.keys():
		dictAvgCoverage[row['Gene']] = float(cov)
		dictMedian[row['Gene']] = float(median)
                dictU10[row['Gene']] = float(u10)
                dictU20[row['Gene']] = float(u20)
                dictU50[row['Gene']] = float(u50)
                dictU100[row['Gene']] = float(u100)
		dictCount[row['Gene']] = 1
	else:
		dictAvgCoverage[row['Gene']] += float(cov)
		dictMedian[row['Gene']] += float(median)
		dictU10[row['Gene']] += float(u10)
		dictU20[row['Gene']] += float(u20)		
		dictU50[row['Gene']] += float(u50)
		dictU100[row['Gene']] += float(u100)

		dictCount[row['Gene']] += 1


for gene in dictAvgCoverage:
	c=dictAvgCoverage[gene]/dictCount[gene]
	d=dictMedian[gene]/dictCount[gene]
	e=dictU10[gene]/dictCount[gene]
	f=dictU20[gene]/dictCount[gene]
	g=dictU50[gene]/dictCount[gene]
	h=dictU100[gene]/dictCount[gene]
	print gene + "\t"+  str(c) + "\t" + str(dictCount[gene]) + "\t"+ str(d) + "\t"+str(e) + "\t"+str(f) + "\t"+ str(g)+ "\t"+str(h) 
