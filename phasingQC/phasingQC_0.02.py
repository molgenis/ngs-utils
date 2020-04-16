
import vcf
import pysam
import numpy as np
import time
import sys
import argparse

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
    # 0 = A|T A|T, 1 = A|A A|A, 2 = A|T T|N,
    # 3 = A|T A|N, 4 = A|T N|A, 5 = A|T T|T && A|A A|N
    # 6 = (A|T or A|A) N|N
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
##I/O

__author__ = 'Raul,Carlos:UMCG'
same_names = True
fast_mode = False # fast mode compares 0|1 and 0|1 instead of A|G A|G...CT|CT T|T is not detected!

parser = argparse.ArgumentParser(description='This is phasing QC.')
parser.add_argument('-i','--reference', help='Reference file VCF',required=True)
parser.add_argument('-I','--toevaluate', help='Check file VCF',required=True)
parser.add_argument('-o','--output',help='Output file table', required=True)
parser.add_argument('-c','--coupling',help='Coupling file', required=False)
parser.add_argument('-f','--fast',help='Fast Mode', action='store_true')
args = parser.parse_args()

if args.coupling:
    same_names = False
if args.fast is not None:
    fast_mode = args.fast

if args.coupling is not None:
    same_names = False
if args.fast is not None:
    fast_mode = True
fast_mode = False # fast mode compares 0|1 and 0|1 instead of A|G A|G...CT|CT T|T is not detected!

# If the sample names are not the same provide a two column tab delimited coupling file CheckID RefID
ref_reader = vcf.Reader(filename=args.reference)
chk_reader = vcf.Reader(filename=args.toevaluate)

checkSamples = np.asarray(chk_reader.samples)
links = []
if same_names:
    for sample in checkSamples:
        links.append([sample,sample])
else:
    for line in open(args.coupling):
        links.append(line[:-1].split("\t")[::-1])
get_ref_id = dict(links)
## where do the closes go?

metric = [[0 for i in range(7)] for j in  range(len(checkSamples))]
rev = [False for i in range(len(checkSamples))]

output = open('metric_per_snp', "w" )
start_time = time.time()##TIMER
recordCounter = 0
previousPosition = None
for record in chk_reader:
    #refsnp = ref_reader.fetch(record.CHROM, record.POS) # if fails to get it then wt? try!
    try:
        refsnp = ref_reader.fetch(record.CHROM, record.POS) # if fails to get it then wt? try!
    except:
        print('Position not found in REF: Skipping')
        print(record.POS)
        continue
    if refsnp.POS != record.POS:
        print('Incorrect fetching ocurred: Skipping')
        print(refsnp.POS,record.POS)
        continue
    if previousPosition == record.POS:
        print('Duplicate position encoutnered: Skipping')
        print(refsnp.POS,record.POS)
        continue
    sampleCounter = 0
    metric_per_snp = [0 for i in range(7)]
    for sample in checkSamples:
        if fast_mode:
            refgt = refsnp.genotype(get_ref_id[sample])['GT']
            chkgt = record.genotype(sample)['GT']
        else:
            refgt = refsnp.genotype(get_ref_id[sample]).gt_bases
            chkgt = record.genotype(sample).gt_bases
        if "/" in (refgt + chkgt):
            print(record.POS,sample)
            print("Only phased data can be matched... skipped")
            print("")
            continue
        if rev[sampleCounter]: #should I swap?
            chkgt = swap(chkgt)
        if refgt != chkgt: # if mismatch
            mismatch_number = mismatch_type(refgt,chkgt) # get mismatch type
            if mismatch_number == None:
                print("")
                continue
            if mismatch_number == 2:
                rev[sampleCounter] = not rev[sampleCounter] # only reverse when A|T T|A
            metric[sampleCounter][mismatch_number] += 1
            metric_per_snp[mismatch_number] += 1
        elif homozygote(refgt):
            metric[sampleCounter][1] += 1
            metric_per_snp[1] += 1
        else:
            metric[sampleCounter][0] += 1
            metric_per_snp[0] += 1
        sampleCounter += 1
    ### Prints a metric per SNPS to detect mismatches that ocurr globally... techical probbly
    previousPosition = record.POS
    output.write(str(record.CHROM)+"\t"+str(record.POS)+"\t"+str(record.ID)+"\t")
    for j in metric_per_snp[:-1]:
        output.write(str(j)+"\t")
    output.write(str(metric_per_snp[-1])+"\n")
    ###
    recordCounter += 1
    
print("Program%f seconds" % (time.time() - start_time))
output.close()

i = 0
with open(args.output, "w" ) as output:
    for each in checkSamples:
        output.write(checkSamples[i]+"\t")
        for j in metric[i][:-1]:
            output.write(str(j)+"\t")
        output.write(str(metric[i][-1])+"\n")
        i += 1
print("written%f seconds" % (time.time() - start_time))
##TIMER
###END
###############################################################################
