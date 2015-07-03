#!/usr/bin/env Rscript
DOC = "This script does a GC correction as described in Fan, H. C., & Quake, S. R. (2010). Sensitivity of noninvasive prenatal detection of fetal aneuploidy from maternal plasma using shotgun sequencing is limited only by counting statistics"

# Retrieve command line parameters
suppressPackageStartupMessages(library("argparser"))

# Create parser object
parser <- arg.parser(DOC, name="Fan & Quake GC correct")

parser <- add.argument(parser, "-gc",  help = ".tsv file with GC percentage per bin.")
parser <- add.argument(parser, "-s",  help = ".RDS file containing a list with indices of GC percentage bins in steps of 0.1%, pre-calculated")
parser <- add.argument(parser, "-c",  help = "Directory where tsv files are stored")
parser <- add.argument(parser, "-o",  help = "Directory where ouput files will be stored")
parser <- add.argument(parser, "-rds", help = "name of RDS file")

#Parse command line arguments
args <- parse.args(parser, argv = commandArgs(trailingOnly = TRUE))


#
## Functions
#
GetFiles <- function(controlDir, strand)
{
  setwd(controlDir)
  files <- list.files( pattern = paste("*.", strand,".*", ".tsv", sep = ""))
  return (sort(files))
}

# Return list with bins of best control files
GetControlFiles = function(control.file.base.name)
{
  control.file.bin = list()
  for (i in 1:length(control.file.base.name) ){
    binfile = as.matrix(read.delim(control.file.base.name[i], header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE))
    
    # Only select relevant chromosomes
    control.file.bin[[i]] = binfile
  }
  
  return (control.file.bin)
}
#returns vector with average number of reads per bin in a GC interval
GetMeanNumberOfReadsGcStep <- function(indices.gc.percentages, sample)
{
  avg.reads.gc.set <- NULL  
  for (i in 1:length(indices.gc.percentages))
  {
    reads <- sample[indices.gc.percentages[[i]]]
    #If total number of reads in GC interval is not 0, determine average number of reads per bin, ignoring bins with no reads
    if(sum(reads != 0))
    {
      avg.reads.gc.set[i] <- sum(reads) /length(reads[reads > 0])
    }
  }
  return(avg.reads.gc.set)
}
#Corrects the number of reads in a given GC interval using the weights
CorrectGcWithWeights <- function(sample, indices.gc.percentages, weights)
{
  for (i in 1:length(indices.gc.percentages))
  {
    #For a given GC interval, subset matching bins and correct with weight
    sample[indices.gc.percentages[[i]]] <- sample[indices.gc.percentages[[i]]] * weights[i]
  }
  return(sample)
}

#
##Script
#

#
##Load input files
#

dir.create(path = args$o)

names.forward <- GetFiles(controlDir = args$c, strand = "forward")
names.reverse <- GetFiles(controlDir = args$c, strand = "reverse")

files.forward <- GetControlFiles(control.file.base.name = names.forward)
files.reverse <- GetControlFiles(control.file.base.name = names.reverse)
#Load matrix with GC percentages per bin
gc.percentages <- as.matrix(read.table(args$gc, sep ="\t"))
gc.percentages <- gc.percentages[,1:4985]

#.RDS file containing a list with indices of GC percentage bins in steps of 0.1%, pre-calculated
indices.gc.percentages <- readRDS(args$s)

corrected.samples.forward <- list()
corrected.samples.reverse <- list()


for (sample.current in 1:length(files.forward))
{
  #Load forward strand sample
  sample.forward <- files.forward[[sample.current]]
  #Use only autosomes
  sample.forward <- sample.forward[1:22,]

  #Same for reverse strand sample 
  sample.reverse <- files.reverse[[sample.current]]
  sample.reverse <- sample.reverse[1:22,]

  #remove reads which have no GC % count
  no.gc.count <- which(gc.percentages < 0)
  sample.forward[no.gc.count] <- 0
  sample.reverse[no.gc.count] <- 0

  #Sum forward and reverse strands
  sample <- sample.forward + sample.reverse

  #Gets mean number of reads per bin for a GC interval
  avg.reads.gc.interval <- GetMeanNumberOfReadsGcStep(indices.gc.percentages, sample)

  #Calculates mean global average of n. of reads per bin
  n.of.reads <- sample[(unlist(indices.gc.percentages))]
  n.of.bins <- length(n.of.reads[n.of.reads > 0]) 
  global.mean.n.of.bins <- sum(sample) / n.of.bins

  #Calculate weights
  weights <- global.mean.n.of.bins / avg.reads.gc.interval

  #Correct read counts using weights, forward and reverse apart
  weight.corrected.sample.forward <- CorrectGcWithWeights(sample.forward, indices.gc.percentages, weights)
  weight.corrected.sample.reverse <- CorrectGcWithWeights(sample.reverse, indices.gc.percentages, weights)
  
  output.name.reverse <- gsub(pattern = ".reverse.", replacement = ".reverse.binGC.corrected.", x = names.reverse[[sample.current]])
  output.name.forward <- gsub(pattern = ".forward.", replacement = ".forward.binGC.corrected.", x = names.forward[[sample.current]])
  
  
  corrected.samples.forward[[sample.current]] <- sample.forward
  corrected.samples.reverse[[sample.current]] <- sample.reverse
  
  #Write output files
  write.table(weight.corrected.sample.forward, paste(args$o, "/", output.name.forward, sep=""), quote = FALSE, sep ="\t", row.names = TRUE)
  write.table(weight.corrected.sample.reverse, paste(args$o, "/", output.name.reverse, sep=""), quote = FALSE, sep ="\t", row.names = TRUE)
}

output.names.reverse <- gsub(pattern = ".reverse.", replacement = ".reverse.binGC.corrected.", x = names.reverse)
output.names.forward <- gsub(pattern = ".forward.", replacement = ".forward.binGC.corrected.", x = names.forward)

rds.files.forward <- list(corrected.samples.forward, output.names.forward)
rds.files.reverse <- list(corrected.samples.reverse, output.names.reverse)

saveRDS(object = rds.files.forward, file = paste("/Users/dirkdeweerd/Graduation/Positive_Controls//rds_files/" , args$rds, "_Forward.rds", sep =""))
saveRDS(object = rds.files.reverse, file = paste("/Users/dirkdeweerd/Graduation/Positive_Controls//rds_files/" , args$rds, "_Reverse.rds", sep =""))
