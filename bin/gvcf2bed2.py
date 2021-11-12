#
# DATA in /groups/umcg-atd/tmp04/projects/QXTR_276-Exoom_v1/run01/results/variants/gVCF
#

import argparse
from collections import namedtuple
import pandas as pd
import vcf
#import csv
import cyvcf2
from cyvcf2 import VCF

# Get GQ value from a cyvcf2 record
def get_format_value(record, format_field, sample_idx):
    """
    Get GQX value from a cyvcf2 record
    :param record: the record
    :param sample_idx: sample index
    :return: float
    """
    try:
        gq_arr = record.format(format_field)
    except KeyError:
        if record.QUAL is not None:
            return record.QUAL
        return None
    if record.QUAL is not None and gq_arr is not None:
        return min(int(gq_arr[sample_idx][0]), record.QUAL)
    elif gq_arr is not None:
        return gq_arr[sample_idx][0]
    elif record.QUAL is not None:
        return record.QUAL
    else:
        return None

# Converts a GQ into a colors(green, orange , red) bedfile row based on set threshold.
def get_color_code(GQ):
    dictColor = {'low': '255,0,0', 'medium': '255,255,0', 'high': '0,255,0'}

    if GQ < dictGQThresholds['low']:
        return dictColor['low']
    elif GQ > dictGQThresholds['medium']:
        return dictColor['high']
    else:
        return dictColor['medium']

def get_Thresholds_GQ_group(qualityScore):
    #dictThresholds = {'low': 20, 'medium': 50, 'high': 50}

    if qualityScore < dictGQThresholds['low']:
        return str("low")
    elif qualityScore > dictGQThresholds['medium']:
        return str("high")
    else:
        return str('medium')

def get_Thresholds_DP_group(depth):

    if depth < dictDPThresholds['low']:
        return str("low")
    else:
        return str('normal')


# Calculate avarge GQ and DP for the tagert given bases on the given gVCF.
def calculate_stats_per_target(chr,start,stop,geneName,gVCF):
    summedDP=0
    summedGQ=0
    summedRegionSize=0
    
    region = chr + ':' + start + '-' + stop
    targetSize=(int(stop)-int(start))
    targetStart=int(start)
    targetEnd=int(stop)
    regionEnd=0
    countZeros=0
    v=None
    dictGQThresholds = {'low': {
                            'GQ':[],
                            'DP':[],
                            'bases':[]},
                      'medium': {
                            'GQ':[],
                            'DP': [],
                            'bases':[]},
                      'high': {
                            'GQ':[],
                            'DP': [],
                            'bases':[]}
                      }
    dictDPThresholds = {'low': {
                            'DP':[],
                            'bases':[]}, 
                      'normal': {
                            'DP': [],
                            'bases':[]}
                      }

    # find variant and for reference or variant rows.
    for v in gVCF(region):
        if v is None:
            continue
        else:
            if "END" in v.INFO:
                regionEnd=v.INFO['END']
                #print("regionEnd: ",regionEnd)
            else:
                regionStart=int(v.start)
                regionEnd=int(v.end)
                if regionStart < targetStart:
                    if regionEnd > targetEnd:
                    ### targetStart and targetEnd are within 1 region --> targetStart and targetEnd are the boundaries
                        numberBases=(targetEnd-targetStart)
                    else:
                    ### targetStart is higher than regionStart --> targetStart is start, regionEnd is end
                        numberBases=(regionEnd-targetStart)
                elif regionEnd > targetEnd:
                ### regionEnd is higher than targetEnd --> regionStart is start, targetEnd is end
       	       	    numberBases=(targetEnd-regionStart)
                else:
                #### regionStart and region stop are within 1 target --> regionStart and regionEnd are the boundaries     
                    numberBases=(regionEnd-regionStart)
                fieldValue=int(get_format_value(v,'DP',0))
                if fieldValue == 0:
                    countZeros=countZeros+numberBases
                summedDP = summedDP + (int(get_format_value(v,'DP',0))*numberBases)
                summedGQ = summedGQ + (int(get_format_value(v, 'GQ', 0))*numberBases)
                summedRegionSize+=numberBases
                thresholdsGQGroup=get_Thresholds_GQ_group(get_format_value(v, 'GQ', 0))
                thresholdsDPGroup=get_Thresholds_DP_group(get_format_value(v, 'DP', 0))
                gq = get_format_value(v, 'GQ', 0)
                dp = get_format_value(v, 'DP', 0)
                dictGQThresholds[thresholdsGQGroup]['GQ'].append(gq)
                dictGQThresholds[thresholdsGQGroup]['DP'].append(dp)
                dictGQThresholds[thresholdsGQGroup]['bases'].append(numberBases)
                dictDPThresholds[thresholdsDPGroup]['DP'].append(dp)
                dictDPThresholds[thresholdsDPGroup]['bases'].append(numberBases)
    if v is not None:
        avgDP=str(summedDP/targetSize)
        for key in dictGQThresholds.keys():
            percentBasesGQInCategoryLow = ((sum(dictGQThresholds['low']['bases'])/targetSize)*100)
            percentBasesGQInCategoryMedium = ((sum(dictGQThresholds['medium']['bases']) / targetSize) * 100)
            percentBasesGQInCategoryHigh = ((sum(dictGQThresholds['high']['bases']) / targetSize) * 100)

        for key in dictDPThresholds.keys():
            percentBasesDPInCategoryLow = ((sum(dictDPThresholds['low']['bases'])/targetSize)*100)
            percentBasesDPInCategoryNormal = ((sum(dictDPThresholds['normal']['bases']) / targetSize) * 100)
            numberBasesDPInCategoryLow = (sum(dictDPThresholds['low']['bases']) )
            return(percentBasesGQInCategoryLow,percentBasesGQInCategoryMedium,percentBasesGQInCategoryHigh,avgDP,percentBasesDPInCategoryLow,percentBasesDPInCategoryNormal,summedDP,summedRegionSize,numberBasesDPInCategoryLow,countZeros)
    else:
        return(0,0,0,0,0,0,0,targetSize,targetSize,targetSize)


#def vcf_record_to_bed(record, bedgraph=False, val=0):
def vcf_record_to_bedgraph(variant):
    if "END" in variant.INFO:

        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}".format("chr"+variant.CHROM, variant.start, variant.INFO['END'],get_format_value(variant,'DP',0),get_format_value(variant,'GQ',0))

    else:
        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}".format("chr"+variant.CHROM, variant.start, variant.end, str(get_format_value(variant,'DP',0)), str(get_format_value(variant,'GQ',0)),'+',variant.start, variant.end,get_color_code(int(get_format_value(variant,'GQ',0))))
        #chr7    127471196  127472363  Pos1  0  +  127471196  127472363  255,0,0

#reads targets from bedfile and feads is it to calculate_stats_per_target(chr, start, stop, gVCF)
def read_targets(bed,gVCF,outputFile):
#   gVCF = cyvcf2.VCF(gVCF)
    content=[]
    with open(outputFile, "w") as ohandle:
        ohandle.write(str("chr\tstart\tstop\tgene\tGQ_low\tGQ_medium\tGQ_high\tavgDp\tpercentage DP<20\tpercentage DP>20\tsummedDP\ttargetSize\tnumberBasesDPInCategoryLow\tnumberBasesZeroCoverage" + "\n" ))
        with open(bed)as f:
            for line in f:
                content= line.strip().split()
                region=content[0]+':'+ content[1]+'-'+content[2]

                if len(content)>3:
                    geneName =content[3]
                else:
                    geneName="unknown"

                response=calculate_stats_per_target(content[0], content[1],content[2],geneName,gVCF)
                ohandle.write(str("{0}\t{1}\t{2}\t{3}\t{4:.2f}%\t{5:.2f}%\t{6:.2f}%\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}".format(content[0],content[1],content[2], geneName, response[0],response[1],response[2],response[3],response[4],response[5],response[6],response[7], response[8],response[9]) + "\n"))

        ohandle.close()
        f.close()

def targets_to_bedgraph(bed,gVCF,outputFile):
    first=True
    with open(outputFile, "w") as ohandle:
        ohandle.write(str("browser position chr22:16258166-16258172" + "\n"))
        ohandle.write(str("browser hide all" + "\n"))
        ohandle.write(str("track name=\"ItemRGBDemo\" description=\"Item RGB demonstration\" itemRgb=\"On\"" + "\n"))
        with open(bed)as f:
            for line in f:
                content= line.strip().split()
                region=content[0]+':'+ content[1]+'-'+content[2]

                if len(content)>3:
                    geneName =content[3]
                else:
                    geneName="unknown"

                if first is True:
                    ohandle.write("browser position chr"+ region + "\n")
                    ohandle.write(str("browser hide all" + "\n"))
                    ohandle.write(str("track name=\"ItemRGBDemo\" description=\"Item RGB demonstration\" itemRgb=\"On\"" + "\n"))
                    first=False

                for variant in gVCF(region):
                    ohandle.write(str(vcf_record_to_bedgraph(variant)) + "\n")



def main():
    desc = """
    Create a BED file from a gVCF. Takes a list of parameters which to parse for into
    output bed.
    """

    # global Thresholds
    global dictGQThresholds
    global dictDPThresholds
    dictGQThresholds = {'low': 20, 'medium': 50, 'high': 50}
    dictDPThresholds = {'low': 20}

    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument("-I", "--input", type=str,
                        required=True, help="Input gVCF")
    parser.add_argument("-O", "--output", type=str,
                        required=True, help="Output bed file")
    parser.add_argument("-s", "--sample", type=str,
                        required=False,
                        help="Sample name in gVCF file to use. "
                             "Will default to first sample ")
    parser.add_argument("-q", "--qualityThresholds", type=list, default=[20,50],
                        help="Genotype quality range per category low, medium, high. default: [20,50]"
                            "This results in GQ categories: low: 0 <-> 20, medium: 21 <-> 50, high: > 51")
    parser.add_argument("-d", "--depthThresholds", type=list, default=[20],
                        help="DepthThreshold, default=[20] this results in GQ categories: low: 0 <-> 20,high: > 20")
    parser.add_argument("-b", "--bedfile", type=str,required=True,help="target bedfile to calculate percentage per quality threshold."
                        "Bedfile must contain colunms: chr \t start \t stop \t gene")
    parser.add_argument("-g", "--bedgraph", type=bool, default=True,
                        help="Output in bedgraph file for visualisation.")
    args = parser.parse_args()

    #Set Global thresholds.
    dictDPThresholds['low']=args.depthThresholds[0]
    dictGQThresholds['low']=args.qualityThresholds[0]
    dictGQThresholds['medium'] = args.qualityThresholds[1]
#    dictGQThresholds['high'] = args.qualityThresholds[2]

    gVCF = cyvcf2.VCF(args.input)
    bed = args.bedfile

    read_targets(bed,gVCF,args.output)
    if not args.sample:
        sample = 0
    else:
        sample = gVCF.samples.index(args.sample)
    outputBedgraphFiles=('.').join(args.output.split('.')[:-1]) +".bedgraph"
    #targets_to_bedgraph(bed,gVCF,outputBedgraphFiles)
    

if __name__ == "__main__":
    main()
