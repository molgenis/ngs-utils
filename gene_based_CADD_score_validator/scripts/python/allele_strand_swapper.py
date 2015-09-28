#! /usr/bin/local/python3
"""
Author  : Bram Miedema
Version : 0.1
date    : 22-September-2015

This script 

This script needs two arguments -o [path_to_output_file]> and a one or more files
any argument after the options below will be considered file names <two or more vcf batch files to be merged>]

Below is a list of options available:
    -h [--help] shows this screen
    -v [--verbose] verbose mode test purpose
    -o [--output] specifies the output file path and file name e.g /homes/output.vcf
    -t [--test] enables test run no output file will be created

report bugs:
    a.t.miedema@st.hanze.nl

"""

from urllib.request import urlopen
import getopt
import sys
import csv
import re
import os

settings = dict()

def parseOpt():
    '''
    parses commandline arguments and set these as global settings
    '''

    opts, args = getopt.getopt(sys.argv[1:], 'vo:ht', ['verbose','output=', 'help', 'test' ])
    for o, a in opts:
        if o in ('-v', '--verbose'):
            settings['verbose'] = 1
        if o in ('-o', '--output'):
            settings['outputFile'] = a
        if o in ('-h', '--help'):
            usage()            
        if o in ('-t', '--test'):
            settings['test'] = True
    return args

def verbose(message):
    '''
    prints verbose messages when verbose is enabled
    '''
    if settings['verbose'] > 0:
        print(message)

def usage():
    ''' 
    prints usage
    '''
    print(__doc__, file = sys.stderr)
    exit(1)

def setup():
    '''
    init setup 
    '''
    settings['verbose'] = 0
    settings['outputFile'] = ""
    settings['test'] = False

def error(err):
    '''
    prints error stream
    '''
    print('Het volgende gaat fout:\n{}'.format(err.args), file = sys.stderr)
    exit(1)
                        
def getHeader():
    header = '''##fileformat=VCFv4.1
##contig=<ID=Un|NT_113961.1,length=166566>
##contig=<ID=Un|NT_113923.1,length=186858>
##contig=<ID=Un|NT_167208.1,length=164239>
##contig=<ID=Un|NT_167209.1,length=137718>
##contig=<ID=Un|NT_167210.1,length=172545>
##contig=<ID=Un|NT_167211.1,length=172294>
##contig=<ID=Un|NT_167212.1,length=172149>
##contig=<ID=Un|NT_113889.1,length=161147>
##contig=<ID=Un|NT_167213.1,length=179198>
##contig=<ID=Un|NT_167214.1,length=161802>
##contig=<ID=Un|NT_167215.1,length=155397>
##contig=<ID=Un|NT_167216.1,length=186861>
##contig=<ID=Un|NT_167217.1,length=180455>
##contig=<ID=Un|NT_167218.1,length=179693>
##contig=<ID=Un|NT_167219.1,length=211173>
##contig=<ID=Un|NT_167221.1,length=128374>
##contig=<ID=Un|NT_167222.1,length=129120>
##contig=<ID=Un|NT_167223.1,length=19913>
##contig=<ID=Un|NT_167224.1,length=43691>
##contig=<ID=Un|NT_167225.1,length=27386>
##contig=<ID=Un|NT_167228.1,length=40531>
##contig=<ID=Un|NT_167230.1,length=41934>
##contig=<ID=Un|NT_167232.1,length=39939>
##contig=<ID=Un|NT_167235.1,length=42152>
##contig=<ID=Un|NT_167236.1,length=43523>
##contig=<ID=Un|NT_167237.1,length=43341>
##contig=<ID=Un|NT_167238.1,length=39929>
##contig=<ID=Un|NT_167240.1,length=38154>
##contig=<ID=Un|NT_167220.1,length=15008>
##contig=<ID=Un|NT_167226.1,length=40652>
##contig=<ID=Un|NT_167227.1,length=45941>
##contig=<ID=Un|NT_167229.1,length=34474>
##contig=<ID=Un|NT_167231.1,length=45867>
##contig=<ID=Un|NT_167233.1,length=33824>
##contig=<ID=Un|NT_167234.1,length=41933>
##contig=<ID=Un|NT_167239.1,length=36651>
##contig=<ID=Un|NT_167241.1,length=36422>
##contig=<ID=Un|NT_167242.1,length=39786>
##contig=<ID=Un|NT_167243.1,length=38502>
##reference=C:\\Program_Files_(x86)\\SoftGenetics\\NextGENe\\References\\GRCh37
##contig=<ID=1,length=249240621>
##contig=<ID=1|NT_113878.1,length=106433>
##contig=<ID=1|NT_167207.1,length=547496>
##contig=<ID=2,length=243189373>
##contig=<ID=3,length=197962430>
##contig=<ID=4,length=191044276>
##contig=<ID=4|NT_113885.1,length=189789>
##contig=<ID=4|NT_113888.1,length=191469>
##contig=<ID=5,length=180905260>
##contig=<ID=6,length=171055067>
##contig=<ID=7,length=159128663>
##contig=<ID=7|NT_113901.1,length=182896>
##contig=<ID=8,length=146304022>
##contig=<ID=8|NT_113909.1,length=38914>
##contig=<ID=8|NT_113907.1,length=37175>
##contig=<ID=9,length=141153431>
##contig=<ID=9|NT_113915.1,length=187035>
##contig=<ID=9|NT_113911.1,length=36148>
##contig=<ID=9|NT_113914.1,length=90085>
##contig=<ID=9|NT_113916.2,length=169874>
##contig=<ID=10,length=135524747>
##contig=<ID=11,length=134946516>
##contig=<ID=11|NT_113921.2,length=40103>
##contig=<ID=12,length=133841895>
##contig=<ID=13,length=115109878>
##contig=<ID=14,length=107289540>
##contig=<ID=15,length=102521392>
##contig=<ID=16,length=90294753>
##contig=<ID=17,length=81195210>
##contig=<ID=17|NT_113943.1,length=81310>
##contig=<ID=17|NT_113930.1,length=174588>
##contig=<ID=17|NT_113941.1,length=37498>
##contig=<ID=17|NT_113945.1,length=41001>
##contig=<ID=18,length=78017248>
##contig=<ID=18|NT_113947.1,length=4262>
##contig=<ID=19,length=59118983>
##contig=<ID=19|NT_113949.1,length=159169>
##contig=<ID=19|NT_113948.1,length=92689>
##contig=<ID=20,length=62965520>
##contig=<ID=21,length=48119895>
##contig=<ID=21|NT_113950.2,length=27682>
##contig=<ID=22,length=51244566>
##contig=<ID=X,length=155260560>
##contig=<ID=Y,length=59363566>
##contig=<ID=UN|NT_113961.1,length=166566>
##contig=<ID=UN|NT_113923.1,length=186858>
##contig=<ID=UN|NT_167208.1,length=164239>
##contig=<ID=UN|NT_167209.1,length=137718>
##contig=<ID=UN|NT_167210.1,length=172545>
##contig=<ID=UN|NT_167211.1,length=172294>
##contig=<ID=UN|NT_167212.1,length=172149>
##contig=<ID=UN|NT_113889.1,length=161147>
##contig=<ID=UN|NT_167213.1,length=179198>
##contig=<ID=UN|NT_167214.1,length=161802>
##contig=<ID=UN|NT_167215.1,length=155397>
##contig=<ID=UN|NT_167216.1,length=186861>
##contig=<ID=UN|NT_167217.1,length=180455>
##contig=<ID=UN|NT_167218.1,length=179693>
##contig=<ID=UN|NT_167219.1,length=211173>
##contig=<ID=UN|NT_167221.1,length=128374>
##contig=<ID=UN|NT_167222.1,length=129120>
##contig=<ID=UN|NT_167223.1,length=19913>
##contig=<ID=UN|NT_167224.1,length=43691>
##contig=<ID=UN|NT_167225.1,length=27386>
##contig=<ID=UN|NT_167228.1,length=40531>
##contig=<ID=UN|NT_167230.1,length=41934>
##contig=<ID=UN|NT_167232.1,length=39939>
##contig=<ID=UN|NT_167235.1,length=42152>
##contig=<ID=UN|NT_167236.1,length=43523>
##contig=<ID=UN|NT_167237.1,length=43341>
##contig=<ID=UN|NT_167238.1,length=39929>
##contig=<ID=UN|NT_167240.1,length=38154>
##contig=<ID=UN|NT_167220.1,length=15008>
##contig=<ID=UN|NT_167226.1,length=40652>
##contig=<ID=UN|NT_167227.1,length=45941>
##contig=<ID=UN|NT_167229.1,length=34474>
##contig=<ID=UN|NT_167231.1,length=45867>
##contig=<ID=UN|NT_167233.1,length=33824>
##contig=<ID=UN|NT_167234.1,length=41933>
##contig=<ID=UN|NT_167239.1,length=36651>
##contig=<ID=UN|NT_167241.1,length=36422>
##contig=<ID=UN|NT_167242.1,length=39786>
##contig=<ID=UN|NT_167243.1,length=38502>
##contig=<ID=M,length=16569>
##INFO=<ID=CLINVAR_CLNSIG,Number=.,Type=String,Description="Value representing clinical significant allele 0 means ref 1 means first alt allele etc.">
##INFO=<ID=CLINVAR_CLNALLE,Number=.,Type=String,Description="Value representing the clinical significanct according to ClinVar">
##INFO=<ID=GENE_SYMBOL,Number=1,Type=String,Description="HGNC gene symbol">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\t'''
    return header 
            
def inflateAggregates(args):
    header = getHeader()
    print(header)
    for arg in args:
        
        try:
            fileName = arg
            fh = open(fileName, 'r')
            
            csvReader = csv.reader(open(arg), delimiter='\t', quotechar='"')
            for row in csvReader:
                
                if "Position" in row[0]:
                    
                    continue
                   
                   # header.append( ",".join(str(x) for x in row))
                    #variant line
                else:
                    # we are only interrested in SNPs and non UTR regions
                    
                        
                    ref = getBaseFromOtherStrand(row[3])
                    alt = getBaseFromOtherStrand(row[4])
                    chr = row[0]
                    position = row[1]
                    info = row[7]
                    id = row[2]
                    quality = row[5]
                    filter = row[6]
                    
                    line = chr + "\t" + position + "\t" + id + "\t" + ref + "\t" + alt + "\t" + quality + "\t" + filter + "\t" + info
                    
                    print(line)

               #18:21153422-21153422
            #   18      21153422        .       C       T       .       .       CLINVAR_CLNSIG="Pathogenic";CLINVAR_CLNALLE=1;GENE_SYMBOL=NPC1;
        except Exception as err:
            error(err)
    #header line to print:  print("\n".join(str(x) for x in header))

def getBaseFromOtherStrand(base):
    if (base == "A"):
        base = "T"
    elif (base == "C"):
        base = "G"
    elif (base == "G"):
        base = "C"
    elif (base =="T"):
        base = "A"
    
    return base


def main():
    setup()
    args = parseOpt()
   # if(settings['outputFile'] == ""):
    #    usage()
    
    if(len(args) > 0):
        inflateAggregates(args)

main()