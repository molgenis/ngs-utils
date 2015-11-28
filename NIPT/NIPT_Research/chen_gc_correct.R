#!/usr/bin/env Rscript
DOC = "This script does a GC correction as described in Chen et al, (2011) Noninvasive prenatal diagnosis of fetal trisomy 18 and trisomy 13 by maternal plasma DNA sequencing."

span <- 0.75
# Retrieve command line parameters
suppressPackageStartupMessages(library("argparser"))

# Create parser object
parser <- arg.parser(DOC, name="Chen et al GC correct")

parser <- add.argument(parser, "-gc",  help = ".tsv file with GC percentage per bin.")
parser <- add.argument(parser, "-f",  help = ".tsv file with reads per bin, forward strand.")
parser <- add.argument(parser, "-r",  help = ".tsv file with reads per bin, reverse strand.")
parser <- add.argument(parser, "-of",  help = ".tsv output file with GC corrected reads per bin, forward strand.")
parser <- add.argument(parser, "-or",  help = ".tsv output file with GC corrected reads per bin, forward strand.")

#Parse command line arguments
args <- parse.args(parser, argv = commandArgs(trailingOnly = TRUE))

#
## Functions
#

GetCorrectionFactor <- function(sample, indices, gc.percentages)
{
  #vectorize matrices
  sample.vector <- as.vector(sample[indices])
  gc.vector <- as.vector(gc.percentages[indices])
  #fit a loess fit to bin counts ~ GC percentages
  fit.loess <-loess(sample.vector ~ gc.vector, span = span)
  #median bincounts
  m <- median(sample.vector)
  #Calculate correction factor
  f <- m / fit.loess$fitted
  
  return(f)
}

#
## Script
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

#Sum forward and reverse strands
sample <- sample.forward + sample.reverse
#Retrieve indices from bins which have a GC percentage and have reads
indices <- which(gc.percentages > 0 & sample > 0)
#Get correction factor
correction.factor <- GetCorrectionFactor(sample, indices, gc.percentages)

#Apply correction factor
sample.forward[indices] <- sample.forward[indices] * correction.factor
sample.reverse[indices] <- sample.reverse[indices] * correction.factor

#Write output files
write.table(sample.forward, args$of, quote = FALSE, sep ="\t", row.names = TRUE)
write.table(sample.reverse, args$or, quote = FALSE, sep ="\t", row.names = TRUE)

