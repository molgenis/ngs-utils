import argparse
import os
import os.path
import math
import sys, getopt
import json
import functools
import locale
import subprocess
import threading
import time
import tempfile
import warnings
import Queue
import pprint

# valid column names
# from http://picard.sourceforge.net/picard-metric-definitions.shtml#RnaSeqMetrics
COL_NAMES = {
    'PF_BASES': 'pfBases',
    'PF_ALIGNED_BASES': 'pfAlignedBases',
    'RIBOSOMAL_BASES': 'ribosomalBases',
    'CODING_BASES': 'codingBases',
    'UTR_BASES': 'utrBases',
    'INTRONIC_BASES': 'intronicBases',
    'INTERGENIC_BASES': 'intergenicBases',
    'IGNORED_READS': 'ignoredReads',
    'CORRECT_STRAND_READS': 'correctStrandReads',
    'INCORRECT_STRAND_READS': 'incorrectStrandReads',
    'PCT_RIBOSOMAL_BASES': 'pctRibosomalBases',
    'PCT_CODING_BASES': 'pctCodingBases',
    'PCT_UTR_BASES': 'pctUtrBases',
    'PCT_INTRONIC_BASES': 'pctIntronicBases',
    'PCT_INTERGENIC_BASES': 'pctIntergenicBases',
    'PCT_MRNA_BASES': 'pctMrnaBases',
    'PCT_USABLE_BASES': 'pctUsableBases',
    'PCT_CORRECT_STRAND_READS': 'pctCorrectStrandReads',
    'MEDIAN_CV_COVERAGE': 'medianCvCoverage',
    'MEDIAN_5PRIME_BIAS': 'median5PrimeBias',
    'MEDIAN_3PRIME_BIAS': 'median3PrimeBias',
    'MEDIAN_5PRIME_TO_3PRIME_BIAS': 'median5PrimeTo3PrimeBias',
}

COL_NAMES_MDUP = {
    'LIBRARY': 'library',
    'UNPAIRED_READS_EXAMINED': 'unpairedReadsExamined',
    'READ_PAIRS_EXAMINED': 'readPairsExamined',
    'UNMAPPED_READS': 'unmappedReads',
    'UNPAIRED_READ_DUPLICATES': 'unmappedReadsDuplicates',
    'READ_PAIR_DUPLICATES': 'readPairsDuplicates',
    'READ_PAIR_OPTICAL_DUPLICATES': 'readPairsOpticalDuplicates',
    'PERCENT_DUPLICATION': 'percentDuplication',
    'ESTIMATED_LIBRARY_SIZE': 'estimatedLibrarySize'}

def meanStd(d):
    """
    Calculate the mean and standard deviation of histogram d.

    @arg d: Histogram of real values.
    @type d: dict[int](float)

    @returns: The mean and standard deviation of d.
    @rtype: tuple(float, float)
    """
    sum_l = 0
    sumSquared_l = 0
    n = 0

    for i in d:
        sum_l += i * d[i]
        sumSquared_l += (d[i] * (i * i))
        n += d[i]
    #for
    
    mean = sum_l / float(n)
    return {'Mean.insertsize' : mean, 'stDev' :+math.sqrt((sumSquared_l / float(n)) - (mean * mean))}
    
def parse_metrics_file(metrics_path):
    """Given a path to a Picard CollectRnaSeqMetrics output file, return a
    dictionary consisting of its column, value mappings.
    """
    data_mark = 'PF_BASES'
    tokens = []
    with open(metrics_path, 'r') as source:
        line = source.readline().strip()
        fsize = os.fstat(source.fileno()).st_size
        while True:
            if not line.startswith(data_mark):
                # encountering EOF before metrics is an error
                if source.tell() == fsize:
                    raise ValueError("Metrics not found inside %r" % \
                            metrics_path)
                line = source.readline().strip()
            else:
                break

        assert line.startswith(data_mark)
        # split header line and append to tokens
        tokens.append(line.split('\t'))
        # and the values (one row after)
        tokens.append(source.readline().strip().split('\t'))
    data = {}
    for col, value in zip(tokens[0], tokens[1]):
        if not value:
            data[COL_NAMES[col]] = None
        elif col.startswith('PCT') or col.startswith('MEDIAN'):
            if value != '?':
                data[COL_NAMES[col]] = float(value)
            else:
                warnings.warn("Undefined value for %s in %s: %s" % (col,
                    metrics_path, value))
                data[COL_NAMES[col]] = None
        else:
            assert col in COL_NAMES, 'Unknown column: %s' % col
            data[COL_NAMES[col]] = int(value)

    return data
    
def parse_mdup_metrics_file(metrics_path):
    """Given a path to a Picard DuplicationMetrics output file, return a
    dictionary consisting of its column, value mappings.
    """
    data_mark = 'LIBRARY'
    tokens = []
    with open(metrics_path, 'r') as source:
        line = source.readline().strip()
        fsize = os.fstat(source.fileno()).st_size
        while True:
            if not line.startswith(data_mark):
                # encountering EOF before metrics is an error
                if source.tell() == fsize:
                    raise ValueError("Metrics not found inside %r" % \
                            metrics_path)
                line = source.readline().strip()
            else:
                break

        assert line.startswith(data_mark)
        # split header line and append to tokens
        tokens.append(line.split('\t'))
        # and the values (one row after)
        tokens.append(source.readline().strip().split('\t'))
    data = {}
    for col, value in zip(tokens[0], tokens[1]):
        if not value:
            data[COL_NAMES_MDUP[col]] = None
        elif col.startswith('PCT') or col.startswith('PERCENT'):
            if value != '?':
                data[COL_NAMES_MDUP[col]] = float(value)
            else:
                warnings.warn("Undefined value for %s in %s: %s" % (col,
                    metrics_path, value))
                data[COL_NAMES_MDUP[col]] = None
        elif col.startswith('LIBRARY'):
            if value != '?':
                data[COL_NAMES_MDUP[col]] = value
            else:
                warnings.warn("Undefined value for %s in %s: %s" % (col,
                    metrics_path, value))
                data[COL_NAMES[col]] = None 
        else:
            assert col in COL_NAMES_MDUP, 'Unknown column: %s' % col
            data[COL_NAMES_MDUP[col]] = int(value)

    return data

def getFlagstat(fsFile):
  """
  Get the number of mapped reads from flagstat.
  """

  with open (fsFile, 'r') as f:
    for line in f:
      line = line.rstrip('\n')
      if "total" in line:
        total = int(line.split(' ')[0])
      if "mapped (" in line:
        mapped = int(line.split(' ')[0])
      if "duplicates" in line:
        dup = int(line.split(' ')[0])
    return {'mapped.reads' : mapped,'total.reads' : total ,'duplicates' : dup}

def parse_Star_Log_File(starLog):
  """
  Parse starLog file into hash.
  """

  with open (starLog, 'r') as f:
    for line in f:
      line = line.rstrip('\n')
      if str(line) != '':
      	total = line.split('|')
      	
      	
      	if len(total) > 1:
      	  #print total[0].strip() +total[1]
      	  print("{0:<40s}{1:<11}".format(total[0].strip(), total[1]))
      	else:
      	   print("{0:<40s}\t".format(total[0].strip()))
    #return {'mapped.reads' : mapped,'total.reads' : total ,'duplicates' : dup}


def getHist(file, begin, end):
  """
  Get the GC content per sequence from the fastqc data file given a sampleid.
  """

  collect = False
  data = {}

  with open (file, 'r') as f:
    for line in f:
      line = line.rstrip('\n')
      if collect and end in line:
        collect = False
      if collect and len(line)>0:
        data[int(line.split('\t')[0])]=float(line.split('\t')[1])
      if begin in line:
        collect = True
  return data

def main(argv):
  """
  Main entry point.
  """

  fastQC1 = ''
  fastQC2 = ''
  insertSizeMetrics= ''
  RnaSeqMetrics= ''
  flagstats = ''
   
  try:
    opts, args = getopt.getopt(argv,"h:1:2:i:r:f:d:s:",["fastQC1=","fastQC2=","insertSizeMetrics=","RnaSeqMetrics=","flagstats=","starLog="])
         
  except getopt.GetoptError:
    print 'pull_RNA_Seq_Stats.py\n -1 <FastQC_1>\n -2 <FastQC_2>\n -r <RnaSeqMetrics>\n -i <insertSizeMetrics>\n -f <flagstats>\n -d  <dupMatrics>\n -s <starLogFile>\n'
    sys.exit(2)
  for opt, arg in opts:
    if opt == '-h':
       print 'pull_RNA_Seq_Stats.py -1 <FastQC_1>\n -2 <FastQC_2>\n -r <RnaSeqMetrics>\n -i <insertSizeMetrics>\n -f <flagstats>\n -d  <dupMatrics>\n -s <starLogFile>\n'
       sys.exit(2)
    elif opt in ("-1", "--fastQC1"):
       fastQC1 = arg
    elif opt in ("-2", "--fastQC2"):
       fastQC2 = arg
    elif opt in ("-r", "--RnaSeqMetrics"):
	   RnaSeqMetrics = arg
    elif opt in ("-i", "--insertSizeMetrics"):
       insertSizeMetrics = arg
    elif opt in ("-d", "--dupMatrics"):
       dupMatrics = arg
    elif opt in ("-s", "--starLog"):
       starLog = arg   
    elif opt in ("-f", "--flagstats"):
       flagstats = arg   
         
 # print 'fastQC1 file is "', fastQC1
 # print 'fastQC2 file is "', fastQC2
 # print 'insertSizeMetrics file is "', insertSizeMetrics
 # print 'flagStatsFile file is "', flagstats

  fqcFileR1Raw = fastQC1
  fqcFileR2Raw = fastQC2
  insertSizeFile = insertSizeMetrics
  RnaSeqMetrics = RnaSeqMetrics
  fsFileFull = flagstats
  dupMatrics = dupMatrics
  starLog = starLog
  data = {}
    
  if os.path.isfile(insertSizeFile) and os.access(insertSizeFile, os.R_OK):
    data['insertSizeHist'] = getHist(insertSizeFile, 'insert_size', 'EOF')
  if os.path.isfile(fqcFileR1Raw) and os.access(fqcFileR1Raw, os.R_OK):
	data['R1_raw_GC'] = getHist(fqcFileR1Raw, '#GC', "END_MODULE")
  if os.path.isfile(fqcFileR2Raw) and os.access(fqcFileR2Raw, os.R_OK):
    data['R2_raw_GC'] = getHist(fqcFileR2Raw, '#GC', "END_MODULE")
  data['map2Full'] = getFlagstat(fsFileFull)
  data['RnaSeqMetrics'] = parse_metrics_file(RnaSeqMetrics)
  data['dupMatrics'] = parse_mdup_metrics_file(dupMatrics)
  data['starLog'] = parse_Star_Log_File(starLog)  
  
  #print FastQC_1 output in tablular format
  if os.path.isfile(fqcFileR1Raw) and os.access(fqcFileR1Raw, os.R_OK):	
    print "\n## FASTQC:READ1 ##\n"
    R1_raw_GC = meanStd(data['R1_raw_GC'])
    for key in R1_raw_GC.keys():
      print("{0:<40s}\t{1:<3.3f}".format(key, R1_raw_GC[key]))
  
  #if seqtype is paired end, print FastQC_2 output in tablular format 
  if os.path.isfile(fqcFileR2Raw) and os.access(fqcFileR2Raw, os.R_OK):  
    print "\n## FASTQC:READ2 ##\n"
    R2_raw_GC = meanStd(data['R2_raw_GC'])
    for key in R2_raw_GC.keys():
      print("{0:<40s}\t{1:<3.3f}".format(key, R2_raw_GC[key]))
  
  #print flagstat stats in tablular format 
  print "\n## SAMTOOLS:FLAGSTAT ##\t\n"
  map2Full = data['map2Full']
  for key in map2Full.keys():
    print("{0:<40s}\t{1:<11}".format(key, map2Full[key]))
  
  #print CollectRnaSeqMetrics stats in tablular format
  print "\n## PICARD:COLLECTRNASEQMETRICS ##\t\n"
  RnaSeqMetrics = data['RnaSeqMetrics']
  for key in RnaSeqMetrics.keys():
    print("{0:<40s}\t{1:<30}".format(key, RnaSeqMetrics[key]))
  
  #print dupMatrics stats in tablular format
  print "\n ## PICARD:MARKDEDUPMATRICS ##\t\n"
  dupMatrics = data['dupMatrics']
  for key in dupMatrics.keys():
    print("{0:<40s}\t{1:<30}".format(key, dupMatrics[key]))
       
  if os.path.isfile(insertSizeFile) and os.access(insertSizeFile, os.R_OK):
    print "\n## PICARD:INSERTSIZEMERTICS ##\t\n"
    insertSizeHist = meanStd(data['insertSizeHist'])  
    for key in insertSizeHist.keys():
      print("{0:<40s}\t{1:<11.3f}".format(key, insertSizeHist[key]))
  
  
if __name__ == "__main__":
  main(sys.argv[1:])

