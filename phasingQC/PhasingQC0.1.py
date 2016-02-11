
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
    ####
    if ref_gt == eva_gt:
        return False #Speedup
    ####
    if not homozygote(ref_gt) and not homozygote(ref_gt):
        if ref_gt == eva_gt[::-1]:
            return True #Speedup
        elif not double_mismatch(ref_gt,eva_gt):
            return False
        elif not double_mismatch(ref_gt,eva_gt[::-1]):
            return True
        return False
    return False

###############################################################################
