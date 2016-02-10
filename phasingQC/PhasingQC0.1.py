import sys
import optparse
import pyvcf
import pysam
import subprocess
import time

### FUNCTIONS
def Homocygote(gt):
    # Is Homocygote?
    #Input Example: "A|A"
    if gt.split("|")[0] != gt.split("|")[1]:
        return False
    else:
        return True

def Double_mismatch(gt1,gt2):
    # Is it double mismatcht: "A|T" "T|A"
    #Input Example: "A|T" "T|A"
    if ( gt1.split("|")[0] == gt2.split("|")[0] or
         gt1.split("|")[1] == gt2.split("|")[1]
        ):
        return False
    else:
        return True

def Switch_authority(ref_gt, eva_gt, reverse=False):
    # Function says if the mismatch deserves to switch the phasing
    # Input Example: "A|A" "A|T" True
    # Output Example: False
    # current limitations: Trialellic SNPs, cannot detect all swaps

    # Diagnostic:
    if len(ref_gt) != 3:
        print("switch authority: Error: Unexpectedly long reference gt")
        return
    if len(eva_gt) != 3:
        print("switch authority: Error: Unexpectedly long evaluated gt")
        return

    if reverse == True:
        eva_gt = eva_gt[::-1]
    ####
    if ref_gt == eva_gt:
        return False #Speedup
    ####
    if not Homocygote(ref_gt) and not Homocygote(ref_gt):
        if ref_gt == eva_gt[::-1]:
            return True #Speedup
        elif not Double_mismatch(ref_gt,eva_gt):
            return False
        elif not Double_mismatch(ref_gt,eva_gt[::-1]):
            return True
        else:
            return False
    else:
        return False

###############################################################################

parser = optparse.OptionParser()
parser.add_option('-f1', '--file1', dest='file1', help='file #1 for comparisson')
parser.add_option('-f2', '--file2', dest='file2', help='file #2 for comparisson')
parser.add_option('-O', '--outputDir', dest='outputDir', help='path to desired output')

(options, args) = parser.parse_args()

vcf_RNA = vcf.Reader(open(f1, 'r'))
