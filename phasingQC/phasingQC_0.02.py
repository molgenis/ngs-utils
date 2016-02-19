
import vcf
import pysam
import numpy as np
import re
import time
import cProfile, pstats,StringIO

### FUNCTIONS
def swap(gt):
    # Returns swapped genotype
    #Input Example: "ACC|A" output: "A|ACC"
    if "/" in gt:
        print("Warning: it is meaningless to swap unphased genotypes")
        return
    if len(gt) != 3:
        return "|".join(gt.split("|")[::-1])
    else:
        return gt[::-1]

def homozygote(gt):
    # Is homozygote?
    #Input Example: "A|A" O: True
    if gt[0] != gt[2]:
        return False
    return True

def double_mismatch(gt1,gt2):
    # Is it double mismatcht: "A|T" "T|A"
    #Input Example: "A|T" "T|A" O: True
    if "/" in gt1 or "/" in gt2:
        print("Only phased data can be matched")
        return
    if ( gt1.split("|")[0] == gt2.split("|")[0] or
         gt1.split("|")[1] == gt2.split("|")[1]
        ):
        return False
    return True

def mismatch_type(ref_gt, eva_gt):
    # Function says the type of mismatch: currently we swich only #2
    # Input Example: "A|A" "A|T" False
    # Output Example:  5
    # 0 = A|T A|T, 1 = A|A A|A, 2 = A|T T|A,
    # 3 = A|T A|X, 4 = A|T X|A, 5 = A|T T|T && A|A A|X
    # 6 = (A|T or A|A) X|X
    # current limitations: Trialellic SNPs

    # Diagnostic:
    if len(ref_gt.split("|")) != 2 and  len(eva_gt.split("|")) != 2:
        print(ref_gt , eva_gt)
        print("mismatch_type: Error: Unexpected genotype format... skipping")
        return
    swap_eva_gt = swap(eva_gt)
    if ref_gt == swap_eva_gt:
        return 2 #Faster
    if double_mismatch(ref_gt,eva_gt):
        if double_mismatch(ref_gt,swap_eva_gt):
            return 6
    #    if ref_gt == swap_eva_gt:
    #        return 2
        return 4
    #if ref_gt == eva_gt:
    #    if homozygote(eva_gt):
    #        return 1
    #    return 0
    if homozygote(ref_gt) or homozygote(eva_gt):
        return 5
    return 3

###############################################################################
###MAIN
# profiling Python code: https://docs.python.org/2/library/profile.html
pr = cProfile.Profile()
pr.enable()
ref_reader = vcf.Reader(open('testRef.vcf.gz', 'r' ))
chk_reader = vcf.Reader(open('testCheck.vcf.gz', 'r'))
## where do the closes go?

checkSamples = np.asarray(chk_reader.samples)

# seven cases of match/mismatch
metric = [[0 for i in range(7)] for j in  range(len(checkSamples))]
rev = [False for i in range(len(checkSamples))]

start_time = time.time()##TIMER

for record in chk_reader:
    refsnp = ref_reader.fetch(record.CHROM, record.POS) # if fails to get it then wt? try!
    sampleCounter = 0
    for sample in checkSamples:
        refgt = refsnp.genotype(sample).gt_bases
        chkgt = record.genotype(sample).gt_bases
        if "/" in (refgt + chkgt):
            print(record.ID,sample)
            print("Only phased data can be matched... skipped")
            print("")
            continue
        if rev[sampleCounter]: #should I swap?
            chkgt = swap(chkgt)
        if refgt != chkgt: # if mismatch
            mismatch_number = mismatch_type(refgt,chkgt) # get mismatch type
            if mismatch_number == None:
                print(record.ID,sample)
                print("")
                continue
            if mismatch_number == 2:
                rev[sampleCounter] = not rev[sampleCounter] # only reverse when A|T T|A
            metric[sampleCounter][mismatch_number] += 1
        elif homozygote(refgt):
            metric[sampleCounter][1] += 1
        else:
            metric[sampleCounter][0] += 1
        sampleCounter += 1
else:
    print("Nothing was skipped :)")
print("Program%f seconds" % (time.time() - start_time))
i = 0
output = open( "PhasingQC_Output.txt", "w" )
for each in checkSamples:
    output.write(checkSamples[i])
    output.write("\t")
    for unc in metric[i]:
        output.write(str(unc))
        output.write("\t")
    output.write("\n")
    i += 1
output.close()
print("written%f seconds" % (time.time() - start_time))
# put off profiling
pr.disable()
s = StringIO.StringIO()
sortby = 'cumulative'
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()
print s.getvalue()
##TIMER
###END
###############################################################################
