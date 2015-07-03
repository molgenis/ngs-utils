
"""

Sorts a VCF (http://vcftools.sourceforge.net/specs.html) file according to position of the variants.
All the sorting happens in memory (tested for 150MB of VCF)
Usage: 
#> python SortVCFFilename.py <inputVCFFilename.vcf> <outputVCFFilename>


alexandros.kanterakis@gmail.com
Genomics Coordination Center / UMC Groningen

"""

import os
import sys

def DeleteFilename(filename=None):
	os.remove(filename)

def ChrToInt(chr):
        if chr == "X": return 23
        if chr == "Y": return 24
        if chr == "M" or chr == "MT": return 25

	try:
        	return int(chr)
	except:
		print "Invalid chromosome value: ", chr
		return 0


def sortFunction(x,y):
        splittedX = x.split("\t")
        splittedY = y.split("\t")

        if splittedX[0] != splittedY[0]:
                return ChrToInt(splittedX[0]) - ChrToInt(splittedY[0])

	if len(splittedX) < 2 and len(splittedY) < 2:
		return 1 if splittedX[0] < splittedY[0] else -1
	if len(splittedX) < 2: return 1
	if len(splittedY) < 2: return -10

	ret = 0

	try:
		ret = int(splittedX[1]) - int(splittedY[1])
        	return ret
	except:
		pass
	
	try:
		a = int(splittedX[1])
	except:
		print "WARNING: Do not know how to sort these fields:", splittedX[1]
		return 1
	
	try:
		a = int(splittedY[1])
	except:
		print "WARNING: Do not know how to sort these fields:", splittedY[1]
		return -1


	

def SortVCFFilename(inputFilename, outputFilename):

        print "Removing the header.."
        command = 'cat ' + inputFilename + ' | grep -v "#"  > ' + outputFilename + '.toSort'
        print command
        os.system(command)

        command = 'cat ' + inputFilename + ' | grep "#"  > ' + outputFilename + '.header'
        print command
        os.system(command)
        print "..Done"

        print "Loading.."
        filelines = open(outputFilename + '.toSort').readlines()
        print "..Done"

        print "Sorting.."
        fileLinesSorted = sorted(filelines, sortFunction)
        print "..Done"

        print "Outputing..."
        fileOutput = open(outputFilename + '.sorted', "w")
        for line in fileLinesSorted: fileOutput.write(line)
        fileOutput.close()
        print "..Done"

        command = 'cat ' + outputFilename + '.header' + ' ' + outputFilename + '.sorted > ' + outputFilename
        print command
        os.system(command)

        DeleteFilename(outputFilename + '.header')
        DeleteFilename(outputFilename + '.toSort')
        DeleteFilename(outputFilename + '.sorted')

if __name__ == "__main__":
	SortVCFFilename(sys.argv[1], sys.argv[2])
