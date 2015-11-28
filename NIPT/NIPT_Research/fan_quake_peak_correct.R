#!/usr/bin/env Rscript
DOC = "This script does a Peak Correction correction as described in Fan, H. C., & Quake, S. R. (2010). Sensitivity of noninvasive prenatal detection of fetal aneuploidy from maternal plasma using shotgun sequencing is limited only by counting statistics"

# Retrieve command line parameters
suppressPackageStartupMessages(library("argparser"))

# Create parser object
parser <- arg.parser(DOC, name="Fan & Quake Peak Correct correct")

parser <- add.argument(parser, "-f",  help = ".tsv file with reads per bin, forward strand.")
parser <- add.argument(parser, "-r",  help = ".tsv file with reads per bin, reverse strand.")
parser <- add.argument(parser, "-of",  help = ".tsv output file with GC corrected reads per bin, forward strand.")
parser <- add.argument(parser, "-or",  help = ".tsv output file with GC corrected reads per bin, forward strand.")

#Parse command line arguments
args <- parse.args(parser, argv = commandArgs(trailingOnly = TRUE))

#
## Functions
#

#Corrects bins to0 by row (chromosome)
GetIndices <-function(sample.row)
{
  # Defines mean of bins with reads
  mean.sample.row <- mean(sample.row[sample.row > 0])
  range <- 1.96 * sd(sample.row[sample.row > 0])
  #Gets indices of bins which are not in range of mean +- 1.96 * sd
  indices <- which(!sample.row > mean.sample.row - range | sample.row > mean.sample.row + range)
  return(indices)
}

#
## Script
#

#Load forward strand sample
sample.forward <- as.matrix(read.table(args$f))
#Use only autosomes
sample.forward <- sample.forward[1:22,]

#Same for reverse strand sample 
sample.reverse <- as.matrix(read.table(args$r))
sample.reverse <- sample.reverse[1:22,]

#Sum forward and reverse strands
sample <- sample.forward + sample.reverse

indices.row <- NULL
indices.row <- apply(sample,1,GetIndices)

# Sets bins to 0, by row (chromosome)
for (i in 1:nrow(sample))
{
  sample.forward[i,indices.row[[i]]] <- 0
  sample.reverse[i,indices.row[[i]]] <- 0
}

#Write output files
write.table(sample.forward, args$of, quote = FALSE, sep ="\t", row.names = TRUE)
write.table(sample.reverse, args$or, quote = FALSE, sep ="\t", row.names = TRUE)
