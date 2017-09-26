
# Comparers samples with an heterozygous snp in a particulair gene to the same sample with an ASE effect. 
## Then creases a file with the genes containing a heterozygous snp, samples with heterozygous snp in the gene, 
### sample with an ASE effect in the gene and the ratio of the hetrozygous snps vs ASE effect.
##

import glob, sys, fileinput, csv, io, vcf, argparse

parser = argparse.ArgumentParser(description='Process chromosome of intrest, and use the correct gtf, vcf and patient files, to compare GT to ASE effect')
parser.add_argument('--chr', help='chromosome to use')
parser.add_argument('--gtf', help='gtf file to use')
parser.add_argument('--vcf', help='vcf file to use')
parser.add_argument('--sd', '--sampledirectory', help='use directory with ASE.txt files per sample')
parser.add_argument('--od', '--outputdirectory', help='the directory to store your output file')
parser.add_argument('--of', '--outputfile', help='the name of the outputfile.csv')

args = parser.parse_args()
print('chr:', args.chr)
print('gtf file:', args.gtf)
print('vcf file:', args.vcf)
print('sample directories:', args.sd)
print('outputfile directory:', args.od)
print('output file name:', args.of)

#
# creates a dictionary with ENSG number as key, with chr, and start and stop location of this gene as value.
#
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
                stop = int(splitline[4]) + 1
                chr = int(splitline[0])
                genelocation.append(chr)
                genelocation.append(start)
                genelocation.append(stop)
                genelocation.append(ENSGnr)
                genelocation.append(range(start, stop))
                dctchr22geneloc [ENSGnr] = genelocation
#
# creates a dictionary with sample name as key of samples with an ASE effect, and ENSG numbers as value per sample.
#
geneAElist = glob.glob(args.sd+'/HG00*.txt')
geneAEpersample ={}
samplenamelist = []
listASEgenes = []

for name in geneAElist:
    nameA = name.split('/')
    nameB = nameA[-1].split('.')
    samplename = nameB[0]
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
    geneAEpersample [samplename] = list
setASEgenes = set(listASEgenes)

#
# Creates a dictionary with key sample name of samples in the vcf file. 
# The positions of the heterozygous genotypes are values.
#
filename = args.vcf
vcf_reader = vcf.Reader(filename = filename)

records = vcf_reader.fetch(args.chr)

dctHetGT = {}

for record in records:
    position = record.POS 
    for call in record.samples:
        sample = call.sample
        genotype = record.genotype(sample) ['GT']
        if sample in samplenamelist:
            if sample not in dctHetGT:
                dctHetGT[sample] = []
            if genotype == '0/1':
                dctHetGT[sample].append(position)

#
# creates a set with all the genes from the vcf file with an heterozygous snp.
#
dctHetGTperGene = {}
listhtrzgenes = []

for sample in dctHetGT:
    for SNP in dctHetGT[sample]:
        for genes in dctchr22geneloc:
            if SNP in dctchr22geneloc[genes][4]:
                listhtrzgenes.append(dctchr22geneloc[genes][3])
                if sample not in dctHetGTperGene:
                    dctHetGTperGene[sample] = []
                    dctHetGTperGene[sample].append(dctchr22geneloc[genes][3])
                else:
                    dctHetGTperGene[sample].append(dctchr22geneloc[genes][3])

sethtrzgenes = set(listhtrzgenes) 

#
# Creates dictionary with key: gene ID of the genese with heterozygous snp. 
# Values: list of samples with heterozygous snp in the gene, a list of sample with an ASE effect in the gene, and the ratio of the hetrozygous snps vs ASE.
#
dctENSGht = {}

for gene in sethtrzgenes:
    dctENSGht[gene] = []
    sampleHT = []
    sampleASE = []
    dctENSGht[gene].append(gene)
    for sample in dctHetGTperGene:
        if gene in dctHetGTperGene[sample]:
            sampleHT.append(sample)
        if gene in geneAEpersample[sample]:
            sampleASE.append(sample)
    dctENSGht[gene].append(sampleHT)
    dctENSGht[gene].append(sampleASE)
    dctENSGht[gene].append(float(len(sampleASE)/len(sampleHT)))

#
# creases a csv file with the genes with a heterozygous snp, samples with heterozygous snp in the gene, 
# sample with an ASE effect in the gene and the ratio of the hetrozygous snps vs ASE.
#
with open(args.od+'/'+args.of, 'w') as outfile:
    outfile.write('ENSGid\tHetrozygous\tASE\tratioHTRZvsASE')
    for gene in dctENSGht:
        outfile.write('\n'+gene+'\t'+','.join(dctENSGht[gene][1])+'\t'+','.join(dctENSGht[gene][2])+'\t'+(str(dctENSGht[gene][3])))

print ('finished')

