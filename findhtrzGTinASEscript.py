#import vcf


# read in GTF, with the positions of each gene
# gtf = '/apps/data/ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf'
# TODO:
# 1.1. open the file (google -> with open(x)
# 1.2. loop through the file
# 1.3. select the gene regions
# 1.4. save the gene region positions

import glob
import sys
#import numpy
import fileinput
import csv
import io
import vcf
import argparse

parser = argparse.ArgumentParser(description='Process chromosome of intrest, and use the correct gtf, vcf and patient files')
parser.add_argument('--chr', help='chromosome to use', default='22')
parser.add_argument('--gtf', help='gtf file to use', default='/Users/kdelange/Documents/VoorNiekASE/Homo_sapiens.GRCh37.75.gtf')
parser.add_argument('--vcf', help='vcf file to use', default='/Users/kdelange/Documents/VoorNiekASE/genotypes/GEUVADIS.chr22.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.noAnnot.genotypes.2000lines.vcf.gz')
parser.add_argument('--sd', '--sampledirectory', help='use directory with ASE.txt files per sample', default= '/Users/kdelange/Documents/VoorNiekASE/geneAE/')

# nargs= '+' gathers all command line arguments into a list. There must also be one or more arguments or an error message will be generated. 1 argument is necessary, if '*' is used 0 or more files are possible

args = parser.parse_args()
print(args.chr)
print(args.gtf)
print(args.vcf)
print(args.sd)

# thinks to parese: geneAElist (for the sampele names line 63, fileHG ASE files line 74,  

#gtf = 'args.gtf'
dctchr22geneloc = {}


with open(args.gtf) as f:
    for line in f:
        genelocation = []
        if line.startswith(args.chr):
            splitline = line.split('\t')
            geneID = splitline[8]
            ENSG = geneID.split('gene_id "')
            ENSGA = ENSG[1].split('"; ')
            ENSGnr = ENSGA[0]
            if splitline[2] == 'exon':
                start = int(splitline[3])
                stop = int(splitline[4]) + 1 # stop position + 1 
                chr = int(splitline[0])
                genelocation.append(chr)
                genelocation.append(start)
                genelocation.append(stop)
                genelocation.append(ENSGnr)
                genelocation.append(range(start, stop)) 
                dctchr22geneloc [ENSGnr] = genelocation
        
                
#print (dctchr22geneloc)

geneAElist = glob.glob(args.sd+'/HG00*.txt')
geneAEpersample ={}
samplenamelist = []
listASEgenes = []

for name in geneAElist:
    nameA = name.split('/')
    nameB = nameA[-1].split('.')
    samplename = nameB[0]
    print('samplename:', samplename)
    print('name:', name)
    print('nameB:', nameB)
    samplenamelist.append(samplename)

    fileHG = args.sd+'/'+samplename+'.geneAE.txt'
    
    with open(fileHG) as g:
         g.readline()
         list = []
         for line in g:
            splitlineB = line.split('\t')
            ENSG = splitlineB[3]
            list.append(ENSG)
            listASEgenes.append(ENSG)
            #Log2_aFC = splitlineB[7]  
            #list.append(Log2_aFC)
    geneAEpersample [samplename] = list #dictionary met als key de sample name en de ENSG nummers als informatie per sample
setASEgenes = set(listASEgenes)
print ('setASEgenes:', setASEgenes)

print ('geneAEpersample:', geneAEpersample) # print alle ENSG nummers van sample HG...... with an allele specific expression            


# vcf reader is helpfull for reading VCF files in Python
#vcf = "/Users/kdelange/Documents/VoorNiekASE/genotypes/GEUVADIS.chr22.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes.vcf.gz"
filename = args.vcf #'/Users/kdelange/Documents/VoorNiekASE/genotypes/GEUVADIS.chr22.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.noAnnot.genotypes.2000lines.vcf.gz' # file with the fisrt 2000 lines, make sure to use the whole file when the script is finisched!!
print ('only the first 2000 lines of the vcf file are used')
vcf_reader = vcf.Reader(filename = filename)

records = vcf_reader.fetch(args.chr) # get all SNPs on chromosome 22

dctHetGT = {}


# loop over all the SNPs
for record in records:# loop over alle records = in de vcf file alle SNPs op chr 22
    position = record.POS # link de positie van de SNP aan de SNP
    for call in record.samples: # kijk naar 1 SNP loop over de samples, begint met 1 SNP op 1 sample
        sample = call.sample # sample = de snp in de record aan een samples gelinkt
        genotype = record.genotype(sample) ['GT']  # neem van de SNPs gelinkt aan een sample, in de vcf file het genotype
        if sample in samplenamelist: #if sample is in de samplelist ->
            if sample not in dctHetGT: # als het sample niet in AllGT (dictionary) zit:
                dctHetGT[sample] = [] #   voeg sample toe als key met lege lijst als value
            if genotype == '0/1':
                #dctHetGT[sample].append(genotype) # append genotype van het sample aan de lijst
                dctHetGT[sample].append(position) # append position of GT to list

print('dctHetGT:', dctHetGT)    

# the script below creates a dictionary with hetrozygous snps per gene with key patients 
dctHetGTperGene = {}
listhtrzgenes = []

for sample in dctHetGT: # per sample uit dctHetGT
    #print ('patientnummer:', sample)
    for SNP in dctHetGT[sample]: # per SNP per sample
        #print ('SNP:', SNP)
        for genes in dctchr22geneloc:# de genes from the dctchr22geneloc
            
            #print (genes)
            if SNP in dctchr22geneloc[genes][4]: # if the SNP from 1 patient is in the range of specific gene:
                #print ('dict htrzGenes:', dctchr22geneloc[genes][3]) # print the gene
                listhtrzgenes.append(dctchr22geneloc[genes][3])
                #print (dctchr22geneloc[genes][4]) # and print the range of the gene
                #print (SNP)
                if sample not in dctHetGTperGene:
                    dctHetGTperGene[sample] = []
                    #print ('lege lijst:', dctHetGTperGene[sample])
                    dctHetGTperGene[sample].append(dctchr22geneloc[genes][3])
                    #print (dctchr22geneloc[genes][3])
                    #dctHetGTperGene[sample].append(SNP)
                else:
                    dctHetGTperGene[sample].append(dctchr22geneloc[genes][3])
                    #dctHetGTperGene[sample].append(SNP)  
                    #print ('volle lijst?:', dctHetGTperGene[sample])   

sethtrzgenes = set(listhtrzgenes) # a set of the hetrozygous genes of all patients

#print ('sethtrzgenes:', sethtrzgenes)           
#print (dctHetGTperGene) #dictionary with hetrozygous snps per gene per patients 
#______________________________________________________________________

dctENSGht = {} #dict with key gene ID, values, hets patients and ASE patients and hets vs ASE ratio

for gene in sethtrzgenes:
    dctENSGht[gene] = []
    sampleHT = []
    sampleASE = []
    #HetASEratio = []
    dctENSGht[gene].append(gene)
    #print (dctENSGht[genes])
    for sample in dctHetGTperGene:
        if gene in dctHetGTperGene[sample]:
            sampleHT.append(sample)
        if gene in geneAEpersample[sample]:
            sampleASE.append(sample)
    dctENSGht[gene].append(sampleHT)        
    dctENSGht[gene].append(sampleASE)
    dctENSGht[gene].append(float(len(sampleASE)/len(sampleHT)))  
    #print ('HetASEration:', len(sampleASE)/len(sampleHT))
    #print ('hets:', len(sampleHT))
    #print ('ASE:', len(sampleASE))         
            
print (dctENSGht) #dict with key gene ID, values, hets patients and ASE patients and hets vs ASE ratio


with open('/Users/kdelange/Documents/VoorNiekASE/GTvsASE.csv', 'w') as outfile:
    outfile.write('ENSGid\tHetrozygous\tASE\tratioHTRZvsASE')
    for gene in dctENSGht:
        outfile.write('\n'+gene+'\t'+','.join(dctENSGht[gene][1])+'\t'+','.join(dctENSGht[gene][2])+'\t'+(str(dctENSGht[gene][3])))
print('/Users/kdelange/Documents/VoorNiekASE/GTvsASE.csv')


#print ('ASEfile HG00107 use original when script is finisched')
for sample in dctHetGTperGene:
    #print ('dctHetGTperGene:', dctHetGTperGene[sample])
    #print ('geneAEpersample:', geneAEpersample[sample])
    #set(dctHetGTperGene[sample])
    #set(geneAEpersample[sample])
    #print (sample)
    #print ('intersection:', (set(geneAEpersample[sample])).intersection(set(dctHetGTperGene[sample])))
    genes = ((set(geneAEpersample[sample])).union(set(dctHetGTperGene[sample])))



#dctHetGTperGene.items()
print ('finished')

# output file, kolom 1: gene hetr, koom 2: samples die 

# genes hets              ASE       ratioASE
# ENSG1 sampleA                     0
# ENSG2 sampleA,sampleB   sampleA   0.5
# ENSG3 sampleB           sampleB   1

# argparse = add all argumets needed to run the script, files needed, chromosoms enz enz..

# 1:hier de hetrozygote SNPs per sample de positie noteren, met behulp van de positie uit de chr22geneloc dictionary de gen naam halen
# 2: de gen naam vergelijken met de geneAEpersample dictionary. noteer gen namen van de heterozygote SNPs
# 3: count for each sample how many genes have at least 1 hetrozygous SNP, save these genes
# 4: note the genes with an ASE (allel specific expresion) 

#chr22geneloc = dictonary with gene name, start, stop of chr 22 
#geneAEpersample = dictionary with sample name and gene name



        # TODO:
        # 3.1 Only use samples that you got from step 2.2 (google continue)
        # 3.2 If the SNP is heterozygous, Look up in which gene the SNP falls
        # 3.3 Save this positie
# TODO
# 4.1 Count for each sample how many genes had at least 1 heterozygous SNP, and for how many of those there was an ASE effect
#
#