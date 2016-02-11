
### FUNCTIONS
def homozygote(gt):
    # Is homozygote?
    #Input Example: "A|A" O: True
    if gt.split("|")[0] != gt.split("|")[1]:
        return False
    return True

def double_mismatch(gt1,gt2):
    # Is it double mismatcht: "A|T" "T|A"
    #Input Example: "A|T" "T|A" O: True
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

def switch_authority(ref_gt, eva_gt, reverse=False):
    # Function says if the mismatch deserves to switch the phasing
    # Input Example: "A|A" "A|T" False
    # Output Example: False
    # current limitations: Trialellic SNPs, cannot detect all swaps
    # Diagnostic:
    if len(ref_gt) != 3 or len(eva_gt) != 3:
        print("switch authority: Error: Unexpectedly long reference genotype")
        return
    if reverse:
        eva_gt = eva_gt[::-1]
    if ref_gt == eva_gt:
        return False 
    if ref_gt == eva_gt[::-1]:
        return True



###############################################################################
