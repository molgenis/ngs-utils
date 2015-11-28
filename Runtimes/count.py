import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument("--input")
parser.add_argument("--output")

args = parser.parse_args()

print "bla" + args.input

file = open(args.input, 'r')
fileW = open(args.output, 'w')
tekst =file.read()
count=0
my_hash={}
for line in tekst.split("\n")[:-1]:
	splitted=line.split('\t')
	if splitted[0] in my_hash.keys():
		my_hash[splitted[0]] += [int(splitted[1])]  
		
	else: 
		my_hash.update({splitted[0]:[int(splitted[1])]})
totalMax=0
totalMean=0
for i in my_hash.keys():
	fileW.write(i +"\t" + str( max(my_hash[i])) + "\t"+  str(sum(my_hash[i])/len(my_hash[i]))+"\n")
	if i == "s03_FastQC" or i == "s14b_CollectHSMetrics" or i == "s14c_CollectGCBiasMetrics" or i == "s14d_CollectBamIndexMetrics" or i == "s11_MakeDedupBamMd5" or i == "s12_SequenomConcordanceCheck" or i == "s13_CoveragePerBase":
		print "skipped " + i
	else:
		if i == "s10a_Delly" or i == "s09a_Delly":
			totalDelly=(totalMax+max(my_hash[i]))
		else:
			totalMax+=max(my_hash[i])
			totalMean+=sum(my_hash[i])/len(my_hash[i])

fileW.write("total time (max): " + str(totalMax) +"\n")
fileW.write("total time (max), skipping steps: " + str(totalMax) + " minutes. That is " + str(totalMax/60) + " hours. That is "+ str((float(totalMax)/float(1440)))  + " days \n")
fileW.write("total time till Delly, including Delly: " + str(totalDelly) + " minutes. That is " + str(totalDelly/60) + " hours. That is "+ str((float(totalDelly)/float(1440)))  + " days \n")
fileW.write("total time (mean): " + str(totalMean))
file.close()
fileW.close()
