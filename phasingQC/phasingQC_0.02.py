
import vcf
import pysam
import numpy as np
import re
import time

### FUNCTIONS
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

def mismatch_type(ref_gt, eva_gt, reverse=False):
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
    if reverse:
        eva_gt = eva_gt[::-1]

    if double_mismatch(ref_gt,eva_gt):
        if double_mismatch(ref_gt,eva_gt[::-1]):
            return 6
        if ref_gt == eva_gt[::-1]:
            return 2
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


ref_reader = vcf.Reader(open('/groups/umcg-wijmenga/tmp04/umcg-raguirre/QCPhasing/testRef.vcf.gz', 'r' ))
chk_reader = vcf.Reader(open('/groups/umcg-wijmenga/tmp04/umcg-raguirre/QCPhasing/testCheck.vcf.gz', 'r'))
output = open( "PhasingQC_Output.txt", "w" )


checkSamples = np.asarray(chk_reader.samples)
refSamples = np.asarray(ref_reader.samples)

metric = [[0 for i in range(7)] for j in  range(len(checkSamples))]
rev = [False for i in range(len(checkSamples))]

start_time = time.time()##TIMER

for record in chk_reader:
    refsnp = ref_reader.fetch(record.CHROM, record.POS) # if fails to get it then wt?
    sampleCounter = 0
    for sample in checkSamples:
        refgt = refsnp.genotype(sample).gt_bases ##add verification  
        chkgt = record.genotype(sample).gt_bases
        if "/" in refgt or "/" in chkgt:
            print(record.ID,sample)
            print("Only phased data can be matched... skipped", end ="\n\n")
            next
        if refgt != chkgt: # if mismatch
            mismatch_number = mismatch_type(refgt,chkgt, rev[sampleCounter]) # get mismatch type
            if mismatch_number == None:
                print(record.ID,sample, end ="\n\n")
                next
            if mismatch_number == 2:
                rev[sampleCounter] = not rev[sampleCounter] # only reverse when A|T T|A
            metric[sampleCounter][mismatch_number] += 1
        if homozygote(refgt):
            metric[sampleCounter][1] += 1
        else:
            metric[sampleCounter][0] += 1
        sampleCounter += 1
else:
    print("UNSUCCESFUL")
i = 0
for each in checkSamples:
    print(checkSamples[i]+ " ", metric, file=output)

print("%f seconds" % (time.time() - start_time))##TIMER
###END
###############################################################################

