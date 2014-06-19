#! /usr/bin/env python
import sys,os,commands
from optparse import OptionParser
from optparse import OptionGroup

##################################################################################################
# AUTHOR:       M.G. Elferink
# DATE:         11-06-2014
# USE:		./Filter_BAMs_NIPT.py  (in the root folder op a NIPT run)
# Purpose:	Filter raw BAM files for tags: X1==0, CM==0, XT==U + removal of duplicated reads. 
#		Filtering and deduplication is performed with the Sambamba tool.
#		Script is modified to structure of NIPT (UMCU)
##################################################################################################


def filter(list,wkdir): 
	if os.path.isfile(str(wkdir)+"dedup.txt") == True:
		os.system("rm dedup.txt")
	for item in list:
			print "removing non-unique and reads with mismatch, suboptimal hits, and removing duplicates"+"\t" + str(item)
			os.system(str(opt.sambamba) +" view -f bam --filter=\"[X1]==0 and [CM]==0 and [XT]==\'U\'\" " +str(item) + " -o "+ str(item)[0:-4]+"temp.bam -t "+str(threads))
			os.system("echo "+str(item[0:-4])+ " >> "+str(wkdir)+"dedup.txt")
			os.system(str(opt.sambamba)+ " markdup -r "+str(item)[0:-4]+"temp.bam"+ " -t "+str(threads)+ " " +str(item[0:-4]) + "_u0mm_subopt.bam 2>> "+str(wkdir)+"dedup.txt")
			os.system("rm "+str(item)[0:-4]+"temp.bam")
			os.system(str(opt.sambamba)+ " index -t "+str(threads)+ " " +str(item[0:-4]) + "_u0mm_subopt.bam")
	return list

def listbams(wkdir):
	list=[]
	list_dir = os.walk(str(wkdir)).next()[1]
	if len(list_dir) != 0:
		for item in list_dir:
			if "_L0" not in item:
				try:
					for file in os.listdir(str(wkdir)+str(item)):
						if "bam" in file and "bai" not in file and "temp" not in file and "subopt" not in file:
							list+= [str(wkdir)+str(item)+"/"+str(file)]
				except:
					pass
	if len(list) == 0:
		print "No BAMs detected"
		sys.exit()
        else:
		filter(list,wkdir)
	
####################################################################################################

if __name__ == "__main__":
	parser = OptionParser();
        group = OptionGroup(parser, "Main options")
        group.add_option("-x", default="./", dest="wkdir", metavar="[STRING]", help="Working directory [default = \"./\"] ")
        group.add_option("-s", default="./sambamba_v0.4.5", dest="sambamba", metavar="[STRING]", help="full path to SamBamBa binary [default = ./sambamba_v0.4.5]")
	group.add_option("-t", default=4, dest="threads", metavar="[INT]", help="number of threads for sambamba [default = 4]")
	parser.add_option_group(group)
        (opt, args) = parser.parse_args()
	
	if int(opt.threads) < 0:
		threads=1
	else:
		threads=int(opt.threads)	
	
	wkdir=str(opt.wkdir)
	listbams(wkdir)

sys.exit("Finished")

