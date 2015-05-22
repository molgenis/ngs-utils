#!/usr/bin/env Rscript
DOC = "This script does a GC correction as described in Fan, H. C., & Quake, S. R. (2010). Sensitivity of noninvasive prenatal detection of fetal aneuploidy from maternal plasma using shotgun sequencing is limited only by counting statistics"

# Retrieve command line parameters
suppressPackageStartupMessages(library("argparser"))

# Create parser object
parser <- arg.parser(DOC, name="Fan & Quake GC correct")

parser <- add.argument(parser, "-gc",  help = ".tsv file with GC percentage per bin.")
parser <- add.argument(parser, "-s",  help = ".RDS file containing a list with indices of GC percentage bins in steps of 0.1%, pre-calculated")
parser <- add.argument(parser, "-f",  help = ".tsv file with reads per bin, forward strand.")
parser <- add.argument(parser, "-r",  help = ".tsv file with reads per bin, reverse strand.")
parser <- add.argument(parser, "-of",  help = ".tsv output file with GC corrected reads per bin, forward strand.")
parser <- add.argument(parser, "-or",  help = ".tsv output file with GC corrected reads per bin, forward strand.")

#Parse command line arguments
args <- parse.args(parser, argv = commandArgs(trailingOnly = TRUE))


#
## Functions
#

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

#Load matrix with GC percentages per bin
gc.percentages <- as.matrix(read.table(args$gc, sep ="\t"))
gc.percentages <- gc.percentages[,1:4985]

#Load forward strand sample
sample.forward <- as.matrix(read.table(args$f))
#Use only autosomes
sample.forward <- sample.forward[1:22,]

#Same for reverse strand sample 
sample.reverse <- as.matrix(read.table(args$r))
sample.reverse <- sample.reverse[1:22,]

#remove reads which have no GC % count
no.gc.count <- which(gc.percentages < 0)
sample.forward[no.gc.count] <- 0
sample.reverse[no.gc.count] <- 0

#Sum forward and reverse strands
sample <- sample.forward + sample.reverse

#.RDS file containing a list with indices of GC percentage bins in steps of 0.1%, pre-calculated
indices.gc.percentages <- readRDS(args$s)

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

#Write output files
write.table(weight.corrected.sample.forward, args$of, quote = FALSE, sep ="\t", row.names = TRUE)
write.table(weight.corrected.sample.reverse, args$or, quote = FALSE, sep ="\t", row.names = TRUE)

