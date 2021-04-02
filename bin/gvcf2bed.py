import argparse
from collections import namedtuple
import vcf
import csv
import logging
import cyvcf2
from cyvcf2 import VCF

# Get format value from a record
def get_format_value(record, format_field, sample_idx):
    try:
        format_field_arr = record.format(format_field)
    except KeyError:
        return None
    if format_field_arr is not None:
        return format_field_arr[sample_idx][0]
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

# Returns the quality Threshold based of global hash dictGQThresholds
def get_Thresholds_GQ_group(qualityScore):

    if qualityScore < dictGQThresholds['low']:
        return str("low")
    elif qualityScore > dictGQThresholds['medium']:
        return str("high")
    else:
        return str('medium')

# Returns the DP Threshold based of global hash dictDPThresholds
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
    regionSize=(int(stop)-int(start))
    regionEnd=0
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
            else:
                regionEnd=v.end

                numberBases=(regionEnd-v.start)
                summedDP = summedDP + (int(get_format_value(v,'DP',0))*numberBases)
                summedGQ = summedGQ + (int(get_format_value(v, 'GQ', 0))*numberBases)
                summedRegionSize=summedRegionSize+numberBases
                thresholdsGQGroup=get_Thresholds_GQ_group(get_format_value(v, 'GQ', 0))
                thresholdsDPGroup=get_Thresholds_DP_group(get_format_value(v, 'DP', 0))
                gq = get_format_value(v, 'GQ', 0)
                dp = get_format_value(v, 'DP', 0)
                dictGQThresholds[thresholdsGQGroup]['GQ'].append(gq)
                dictGQThresholds[thresholdsGQGroup]['DP'].append(dp)
                dictGQThresholds[thresholdsGQGroup]['bases'].append(numberBases)
                dictDPThresholds[thresholdsDPGroup]['DP'].append(dp)
                dictDPThresholds[thresholdsDPGroup]['bases'].append(numberBases)

                logger.debug ("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format( v.CHROM, v.start, regionEnd,get_format_value(v,'DP',0),get_format_value(v,'GQ',0),int(get_format_value(v, 'GQ', 0))*numberBases, numberBases))
    if v is not None:
        logger.debug("average DP over regio:" +region +"=" +str(summedDP/summedRegionSize))
        avgDP=summedDP/summedRegionSize

        logger.debug("average GQ over regio:" + region + "="+ str(summedGQ/summedRegionSize))
        logger.debug("region size = " +str(summedRegionSize))

        for key in dictGQThresholds.keys():
            percentBasesGQInCategoryLow = ((sum(dictGQThresholds['low']['bases'])/summedRegionSize)*100)
            percentBasesGQInCategoryMedium = ((sum(dictGQThresholds['medium']['bases']) / summedRegionSize) * 100)
            percentBasesGQInCategoryHigh = ((sum(dictGQThresholds['high']['bases']) / summedRegionSize) * 100)

        for key in dictDPThresholds.keys():
            percentBasesDPInCategoryLow = ((sum(dictDPThresholds['low']['bases'])/summedRegionSize)*100)
            percentBasesDPInCategoryNormal = ((sum(dictDPThresholds['normal']['bases']) / summedRegionSize) * 100)

        return(percentBasesGQInCategoryLow,percentBasesGQInCategoryMedium,percentBasesGQInCategoryHigh,avgDP,percentBasesDPInCategoryLow,percentBasesDPInCategoryNormal)
    else:
        return("Nothing to do here.")


#def vcf_record_to_bed(record, bedgraph=False, val=0):
def vcf_record_to_bedgraph(variant):
    if "END" in variant.INFO:
        regionEnd=variant.INFO['END']
    else:
        regionEnd=variant.end

        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}".format("chr"+variant.CHROM, variant.start, regionEnd, int(get_format_value(variant,'DP',0)), int(get_format_value(variant,'GQ',0)),'+',variant.start, regionEnd,get_color_code(int(get_format_value(variant,'GQ',0))))

#reads targets from bedfile and feads is it to calculate_stats_per_target(chr, start, stop, gVCF)
def read_targets(bed,gVCF,outputFile):

    content=[]
    with open(outputFile, "w") as ohandle:
        ohandle.write(str("chr \t start \t stop \t gene \t GQ_low \t GQ_medium \t GQ_high \t DP_avg \t DP_low\t DP_normal" + "\n"))
        logger.info("chr \t start \t stop \t gene \t GQ_low \t GQ_medium \t GQ_high \t DP_avg \t DP_low\t DP_normal")
        with open(bed)as f:
            for line in f:
                content= line.strip().split()
                region=content[0]+':'+ content[1]+'-'+content[2]
                #print(region)

                if len(content)>3:
                    geneName =content[3]
                else:
                    geneName="unknown"

                response=calculate_stats_per_target(content[0], content[1],content[2],geneName,gVCF)

                ohandle.write(str("{0} \t {1} \t {2} \t {3} \t{4:.2f}% \t {5:.2f}% \t {6:.2f}% \t {7:.2f} \t {8:.2f}% \t {9:.2f}%".format(content[0],content[1],content[2], geneName, response[0],response[1],response[2],response[3],response[4],response[5]) + "\n"))
                logger.info("{0} \t {1} \t {2} \t {3} \t{4:.2f}% \t {5:.2f}% \t {6:.2f}% \t {7:.2f} \t {8:.2f}% \t {9:.2f}%".format(content[0],content[1],content[2], geneName, response[0],response[1],response[2],response[3],response[4],response[5]))
        ohandle.close()
        f.close()

def targets_to_bedgraph(bed,gVCF,outputFile):
    first=True
    logger.info("Creating bedfile for gVCF visualistion.")
    with open(outputFile, "w") as ohandle:
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
                    logger.debug(str(vcf_record_to_bedgraph(variant)))
                    ohandle.write(str(vcf_record_to_bedgraph(variant)) + "\n")


def main():
    desc = """
    Create a BED file from a gVCF. Takes a list of parameters which to parse for into
    output bed.
    """

    # global Thresholds
    global dictGQThresholds
    global dictDPThresholds
    dictGQThresholds = {'low': 20, 'medium' :50, 'high':50}
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
                        help="DepthThreshold, default=[20] this results in Depth categories: low: 0 <-> 20, normal: > 20")
    parser.add_argument("-b", "--bedfile", type=str,required=True,help="target bedfile to calculate percentage per quality threshold."
                        "Bedfile must contain colunms: chr \t start \t stop \t gene")
    parser.add_argument("-g", "--bedgraph", type=bool, default=True,
                        help="Output in bedgraph file for visualisation.")
    parser.add_argument("-l", "--logLevel", type=str , default='WARN',
                        help="Loglevels: DEBUG, INFO, WARN. default: WARN")
    args = parser.parse_args()

    # Overwrite Global thresholds.
    if args.depthThresholds:
        #Overwrite Global thresholds.
        dictDPThresholds['low']=args.depthThresholds[0]
    if args.qualityThresholds:
        dictGQThresholds['low']=args.qualityThresholds[0]
        dictGQThresholds['medium'] = args.qualityThresholds[1]
        dictGQThresholds['high'] = args.qualityThresholds[1]

    # Setup logger
    loglevel = args.logLevel
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % loglevel)
    logging.basicConfig(format='%(levelname)s :: %(funcName)s :: lineNr:%(lineno)s :: %(message)s',level=numeric_level)
    global logger
    logger = logging.getLogger()

    gVCF = cyvcf2.VCF(args.input)
    bed = args.bedfile

    read_targets(bed,gVCF,args.output)
    if not args.sample:
        sample = 0
    else:
        sample = gVCF.samples.index(args.sample)

    outputBedgraphFiles=('.').join(args.output.split('.')[:-1]) +".bedgraph"
    targets_to_bedgraph(bed,gVCF,outputBedgraphFiles)


if __name__ == "__main__":
    main()
