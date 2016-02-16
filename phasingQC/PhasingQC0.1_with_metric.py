import vcf
import pysam

### FUNCTIONS
def homozygote(gt):
    # Is homozygote?
    #Input Example: "A|A" O: True
    if gt[0] != gt[1]:
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
    if len(ref_gt) != 3 or len(eva_gt) != 3:
        print("mismatch_type: Error: Unexpectedly long reference genotype")
        return
    if reverse:
        eva_gt = eva_gt[::-1]
        
    if double_mismatch(ref_gt,eva_gt):
        if double_mismatch(ref_gt,eva_gt[::-1]):
            return 6
        if ref_gt == eva_gt[::-1]:
            return 2
        return 4
    if ref_gt == eva_gt:
        if homozygote(eva_gt):
            return 1
        return 0
    if homozygote(ref_gt) or homozygote(eva_gt):
        return 5
    return 3

tests = ["A|A","A|T","T|T","A|X","X|X","X|T"]
for reference in tests[:2]:
    for evaluation in tests:
        print("R: " + reference)
        print("E: " + evaluation + " "
              + str(mismatch_type(reference, evaluation, False))
              )
        print("E: " + evaluation[::-1] + " "
              + str(mismatch_type(reference, evaluation[::-1], False))
              )
        print("")

###############################################################################

ref_reader = vcf.Reader(open('Example1.vcf.gz', 'r'))
chk_reader = vcf.Reader(open('Example2.vcf.gz', 'r'))

metric = [0] * 7
rev = False
for record in chk_reader:
    refsnp = vcf_reader.fetch(record.CHROM, record.POS) # if fails to get it then wt?
    refgt = refsnp.genotype('SampleName2').gt_bases
    chkgt = record.genotype('SampleName1').gt_bases
    if "/" in refgt or "/" in chkgt:
        print("Only phased data can be matched")
        return
    if refgt != chkgt: # if mismatch
        mismatch_number = mismatch_type(refgt,chkgt, rev) # get mismatch type
        if mismatch_number == 2:
            rev = not rev # only reverse when A|T T|A
        n[ mismatch_number ] += 1
    if homozygote(refgt):
        n[1] += 1
print (metric)
