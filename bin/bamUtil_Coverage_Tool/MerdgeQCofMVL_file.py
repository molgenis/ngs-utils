
#note 10-01-2018
#gem cov, 0x en onder de 30x
#per gen een grafiekje maken met onder en boven lijn voor cov en cov van samples/positie zelf.

import glob, fileinput, math, argparse, statistics

parser = argparse.ArgumentParser(description='Process chromosome of intrest, and use the correct gtf, vcf and patient files, to compare GT to ASE effect')
parser.add_argument('--inn', help='BAM util output file, input for python script, ${TMPDIR}/output_final/')
parser.add_argument('--ff', help='directory to store final output txt files ${TMPDIR}/output_final_txt_files/')

args = parser.parse_args()
print('inputfiles:', args.inn)
print ('directory output files:', args.ff)

outputlist = glob.glob((args.inn)+'/*_gene.txt')
print ('number of files:', len(outputlist))
## Makes a dictionary with key variant position and values the coverage of all patient in inputlist.
QCdict = {}

for inputQCfile in outputlist:
    with open(inputQCfile) as h:
        for line in h:
            splitlineA = line.rstrip('\n')
            splitlineB = splitlineA.split('\t')
            VariantPosition = splitlineB[3]+'_'+splitlineB[0]+':'+splitlineB[1]
            if VariantPosition not in QCdict:
                QCdict[VariantPosition] = []
                QCdict[VariantPosition].append(int(splitlineB[2]))
            else:
                QCdict[VariantPosition].append(int(splitlineB[2]))
#print (QCdict)

## makes a dictionary with key variant position and the sum of the coverage of all patients in the input list divided by 50.
## and print the variants with coverage < 1, <10, <30
gemdict = {}
Covdict30 = {}


for pos in QCdict:
    count = 0
    SDposition = statistics.pstdev(QCdict[pos])
    SEposition = SDposition/math.sqrt(len(outputlist))
    medianpos = statistics.median(QCdict[pos])
    tellerpos = 0
    teller30 = 0
    teller20 = 0
    teller10 = 0
    for cov in QCdict[pos]:
        count+=int(cov)
        if cov > 30:
            teller30+=1
        if cov > 20:
            teller20+=1
        if cov > 10:
            teller10+=1
        tellerpos+=1
    per30x = (float(teller30)/float(tellerpos))*100
    per20x = (float(teller20)/float(tellerpos))*100
    per10x = (float(teller10)/float(tellerpos))*100
    minSDcov = (count/len(outputlist))-(SDposition*3)
    maxSDcov = (count/len(outputlist))+(SDposition*3)
    gemcount = count/len(outputlist)
    gemdict[pos] = ["%.2f" % gemcount, "%.2f" % SDposition, round(minSDcov, 0), round(maxSDcov, 0), "%.2f" % SEposition, "%.2f" % int(per30x), "%.2f" % int(per20x), "%.2f" % int(per10x), medianpos]

#print (gemdict)


# it will write a txt file with the gene, position, average coverage, SD , and percentage of samples with coverage above 30x, 20x and 10x.
with open((args.ff)+'/Avrg_SD_moreThan30x20x10x.txt', 'w') as averagefile:
    averagefile.write('Gene\tPosition\tAverageCov\tSD\tMoreThan30x\tMoreThan20x\tMoreThan10x')
    for pos in gemdict:
        name = pos.split('_')
        gene = name[0]
        position = name[1]
        averagefile.write('\n'+gene+'\t'+position+'\t'+str(gemdict[pos][0])+'\t'+str(gemdict[pos][1])+'\t'+str(gemdict[pos][5])+'\t'+str(gemdict[pos][6])+'\t'+str(gemdict[pos][7]))

#file with all the position with average coverage of below 1x, 10x and 30x. 
with open((args.ff)+'/PosBelow_30x_10x_1x.txt', 'w') as outfile:
    outfile.write('Gene\tPosition\tAverageCov\tCov<')
    for pos in gemdict:
        name = pos.split('_')
        gene = name[0]
        position = name[1]
        if float(gemdict[pos][0]) < 1.00:
            outfile.write('\n'+gene+'\t'+position+'\t'+str(gemdict[pos][0])+'\t'+'<1x')
        elif float(gemdict[pos][0]) < 10.00:
            outfile.write('\n'+gene+'\t'+position+'\t'+str(gemdict[pos][0])+'\t'+'<10x')
        elif float(gemdict[pos][0]) < 30.00:
            outfile.write('\n'+gene+'\t'+position+'\t'+str(gemdict[pos][0])+'\t'+'<30x')

# file for Roan in Molgenis
with open((args.ff)+'/Avrg_SD_morethan20x10x_Molgenis.txt', 'w') as averagefile:
    averagefile.write('Chr\tStart\tStop\tGene\tAvgCoverage\tSD\tMoreThan10x\tMoreThan20x\tMedian')
    for pos in gemdict:
        name = pos.split('_')
        gene = name[0]
        positionchr = name[1].split(':')
        chr = positionchr[0]
        positionstart = int(positionchr[1])+1
        positionstop = int(positionchr[1])+1
        averagefile.write('\n'+chr+'\t'+str(positionstart)+'\t'+str(positionstop)+'\t'+gene+'\t'+str(gemdict[pos][0])+'\t'+str(gemdict[pos][1])+'\t'+str(gemdict[pos][7])+'\t'+str(gemdict[pos][6])+'\t'+str(gemdict[pos][8]))

# makes a txt file with position of patient with a minimal or maximal coverage 3xSD outside the average coverage.
with open((args.ff)+'/MVL_outliers.txt', 'w') as SDfile:
    SDfile.write('Gene\tPosition\tminCov3xSD\tmaxCov3xSD\tCoverage\tavgCov\tSD\tSE')
    for pos in gemdict:
        name = pos.split('_')
        gene = name[0]
        position = name[1]
        for cov in QCdict[pos]:
            if cov not in range(int(gemdict[pos][2]), int(gemdict[pos][3])):
                #print (cov, pos, range(int(gemdict[pos][2]), int(gemdict[pos][3])))
                SDfile.write('\n'+gene+'\t'+position+'\t'+str(gemdict[pos][2])+'\t'+str(gemdict[pos][3])+'\t'+str(cov)+'\t'+str(gemdict[pos][0])+'\t'+str(gemdict[pos][1])+'\t'+str(gemdict[pos][4]))





