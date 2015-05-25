#!/usr/bin/env Rscript
DOC = "This script does a GC correction as described in Chen et al, (2011) Noninvasive prenatal diagnosis of fetal trisomy 18 and trisomy 13 by maternal plasma DNA sequencing."

span <- 0.75
# Retrieve command line parameters
suppressPackageStartupMessages(library("argparser"))

# Create parser object
parser <- arg.parser(DOC, name="Chen et al GC correct")

parser <- add.argument(parser, "-gc",  help = ".tsv file with GC percentage per bin.")
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
GetControlFiles = function(control.file.base.name, control.dir)
{
  control.file.bin = list()
  for (i in 1:length(control.file.base.name) ){
    binfile = as.matrix(read.delim(control.file.base.name[i], header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE))
    
    # Only select relevant chromosomes
    control.file.bin[[i]] = binfile
  }
  
  return (control.file.bin)
}
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

names.forward <- GetFiles(controlDir = args$c, strand = "forward")
names.reverse <- GetFiles(controlDir = args$c, strand = "reverse")

files.forward <- GetControlFiles(control.file.base.name = names.forward, control.dir = args$s)
files.reverse <- GetControlFiles(control.file.base.name = names.reverse, control.dir = args$s)

corrected.samples.forward <- list()
corrected.samples.reverse <- list()

dir.create(path = args$o)

for (sample.current in 1:length(files.forward))
  {
  #Load forward strand sample
  sample.forward <- as.matrix(files.forward[[sample.current]])
  #Use only autosomes
  sample.forward <- sample.forward[1:22,]

  #Same for reverse strand sample 
  sample.reverse <- as.matrix(files.reverse[[sample.current]])
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
  
  output.name.reverse <- gsub(pattern = ".reverse.", replacement = ".reverse.loessGC.corrected.", x = names.reverse[[sample.current]])
  output.name.forward <- gsub(pattern = ".forward.", replacement = ".forward.loessGC.corrected.", x = names.forward[[sample.current]])
  
  corrected.samples.forward[[sample.current]] <- sample.forward
  corrected.samples.reverse[[sample.current]] <- sample.reverse
  
  #Write output files
  write.table(sample.forward, paste(args$o, "/", output.name.forward, sep=""), quote = FALSE, sep ="\t", row.names = TRUE)
  write.table(sample.reverse, paste(args$o, "/", output.name.reverse, sep=""), quote = FALSE, sep ="\t", row.names = TRUE)
}

output.names.reverse <- gsub(pattern = ".reverse.", replacement = ".reverse.loessGC.corrected.", x = names.reverse)
output.names.forward <- gsub(pattern = ".forward.", replacement = ".forward.loessGC.corrected.", x = names.forward)

rds.files.forward <- list(corrected.samples.forward, output.names.forward)
rds.files.reverse <- list(corrected.samples.reverse, output.names.reverse)

saveRDS(object = rds.files.forward, file = paste("/Users/dirkdeweerd/Graduation/Positive_Controls/rds_files/" , args$rds, "_Forward.rds", sep =""))
saveRDS(object = rds.files.reverse, file = paste("/Users/dirkdeweerd/Graduation/Positive_Controls/rds_files/" , args$rds, "_Reverse.rds", sep =""))

