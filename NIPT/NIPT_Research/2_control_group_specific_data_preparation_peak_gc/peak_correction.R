#!/usr/bin/env Rscript
DOC = "This script does a Peak Correction correction as described in Fan, H. C., & Quake, S. R. (2010). Sensitivity of noninvasive prenatal detection of fetal aneuploidy from maternal plasma using shotgun sequencing is limited only by counting statistics"

# Retrieve command line parameters
suppressPackageStartupMessages(library("argparser"))

# Create parser object
parser <- arg.parser(DOC, name="Fan & Quake Peak Correct correct")

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
  return (files)
}
#Loads the .tsv binned files into a list and returns this list (controlFiles)
LoadFiles <- function(files)
{
  controlFiles <- list()
  for (l in 1:length(files) ){
    binfile <- read.delim(files[l], header = TRUE, sep = "\t", quote = "\"",
                          dec = ".", fill = TRUE)
    controlFiles[[l]]<- binfile[1:22,]
  }
  return (controlFiles)
}

#Corrects bins to0 by row (chromosome)
GetIndices <-function(sample.row)
{
  # Defines mean of bins with reads
  mean.sample.row <- mean(sample.row[sample.row > 0])
  range <- 1.96 * sd(sample.row[sample.row > 0])
  #Gets indices of bins which are not in range of mean +- 1.96 * sd
  indices <- which(sample.row > 0 & !sample.row > mean.sample.row - range | sample.row > mean.sample.row + range)
  return(indices)
}

#
## Script
#

names.forward <- GetFiles(controlDir = args$c, strand = "forward")
names.reverse <- GetFiles(controlDir = args$c, strand = "reverse")

files.forward <- LoadFiles(files = names.forward)
files.reverse <- LoadFiles(files = names.reverse)

bins.list <- list()

dir.create(path = args$o)
for (chromo in 1:length(files.forward))
{
  sample <- files.forward[[chromo]] + files.reverse[[chromo]]
  indices.row <- NULL
  indices.row <- apply(sample,1,GetIndices)
  bins.list[[chromo]] <- indices.row
}
chromosome.bins.list <- unlist(bins.list, recursive = FALSE)
bins.tabular <- list()

for (chromosome in 1:22)
{
  chromosome.current = seq(from=chromosome, by=22, to=(length(chromosome.bins.list) - (22-chromosome)))
  bins.tabular[[chromosome]] <- as.data.frame(table(unlist(chromosome.bins.list[chromosome.current])))
}

bins.correction.list <- list()
for (chromosome in 1:length(bins.tabular))
{
  peak.bins.chromosome <- bins.tabular[[chromosome]]
  bins.correction.list[[chromosome]] <- as.numeric(as.vector(peak.bins.chromosome[peak.bins.chromosome[,2] >= (length(names.forward) * 0.95), 1]))
}

corrected.samples.forward <- list()
corrected.samples.reverse <- list()

for (sample in 1:length(names.forward))
{
  sample.forward <- files.forward[[sample]]
  sample.reverse <- files.reverse[[sample]]
  for (chromosome in 1:length(bins.correction.list))
  {
    sample.forward[chromosome, bins.correction.list[[chromosome]]] <- 0
    sample.reverse[chromosome, bins.correction.list[[chromosome]]] <- 0
  }
  
  output.name.reverse <- gsub(pattern = ".reverse.", replacement = ".reverse.peak.", x = names.reverse[[sample]])
  output.name.forward <- gsub(pattern = ".forward.", replacement = ".forward.peak.", x = names.forward[[sample]])
  
  corrected.samples.forward[[sample]] <- sample.forward
  corrected.samples.reverse[[sample]] <- sample.reverse
  #Write output files
  write.table(sample.forward, paste(args$o, "/", output.name.forward, sep=""), quote = FALSE, sep ="\t", row.names = TRUE)
  write.table(sample.reverse, paste(args$o, "/", output.name.reverse, sep=""), quote = FALSE, sep ="\t", row.names = TRUE)
}

saveRDS(object = bins.correction.list, file = paste("/Users/dirkdeweerd/Graduation/Peak_Correction/", args$rds, ".rds", sep =""))


output.names.forward <- gsub(pattern = ".forward.", replacement = ".forward.peak.", x = names.forward)
output.names.reverse <- gsub(pattern = ".reverse.", replacement = ".reverse.peak.", x = names.reverse)

rds.files.forward <- list(corrected.samples.forward, output.names.forward)
rds.files.reverse <- list(corrected.samples.reverse, output.names.reverse)

saveRDS(object = rds.files.forward, file = paste("/Users/dirkdeweerd/Graduation/Control_Group_Final//rds_files/" , args$rds, "_Forward.rds", sep =""))
saveRDS(object = rds.files.reverse, file = paste("/Users/dirkdeweerd/Graduation/Control_Group_Final//rds_files/" , args$rds, "_Reverse.rds", sep =""))



