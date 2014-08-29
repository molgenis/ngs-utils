#################################################################################################################

#FUNCTIONS#
#Gets the filenames from the binned control files. The input is the directery where the .tsv files live and the strand
#(forward or reverse)
GetFiles <- function(controlDir, best50, strand)
{
  #Paste the strand + suffix bins.table.tsv to the control file name
  best50 <- paste(best50,".", strand, ".bins.table.tsv", sep = "")
  return (best50)
}
#Loads the .tsv binned files into a list and returns this list (controlFiles)
LoadFiles <- function(files, controlDir)
{
  setwd(controlDir)  
  controlFiles <- list()
  for (l in 1:length(files) ){
    binfile <- read.delim(files[l], header = TRUE, sep = "\t", quote = "\"",
                          dec = ".", fill = TRUE)
    controlFiles[[l]]<- binfile
  }
  return (controlFiles)
}
#Makes the control file by summing all the .tsv matrices into a single matrix and returns this matrix
MakeControlFile <- function(controlFiles, sample)
{
  #Makes an empty matrix. 
  expectFile <- matrix(0, nrow=24, ncol=4985)
  subtotal <- matrix(0, nrow=24, ncol=4985)
  autosomeReadsSampleReads <- sum(sample[1:22,]) / 57614
  expectFile <- sample / autosomeReadsSampleReads 
  #In this loop all the control files are summed into a single file
  for (l in 1:length(controlFiles) )
  {
    file <- controlFiles[[l]]  
    subtotal <- subtotal + file
    sampleReads <- sum(file[1:22]) / 57614
    expectFraction <- file / sampleReads
    expectFile <- expectFile + expectFraction
  }
  totals <<- subtotal
  expectFile <- expectFile / 50
  return (expectFile)
}
#This functions applies the filter for overdispersed bins. There is nothing returned, the results are written to disk
#in tsv format in an intermediate directory
FilterBins <- function(strand, sample, avgFile, chidir, total, sampleID)
{
  totalFile <- matrix(0, nrow=24, ncol=4985)
  filter <- matrix(FALSE, nrow=24, ncol=4985)
  sampleCorrected <- as.matrix(sample)
  chiFilesList <- list()
  chiFileNames <- vector(mode="character")
  #Iterates over the files, applies the filter and writes the output to disk
  for (l in 1:length(controlFiles)  )
  {
    binFile <- controlFiles[[l]]
    autosomeReads <- sum(binFile[1:22, ]) 
    autosomeReadsCorrected <- autosomeReads / 57614
    observed <- (binFile / autosomeReadsCorrected )  * total
    tussen <- (avgFile - observed) ^ 2 / avgFile
    totalFile <- totalFile + tussen
  }
  setwd(chidir)
  totalFile <- totalFile / 50
  filter[totalFile > 3.5] <- TRUE
  sampleCorrected[filter] <- as.matrix(sample[filter]) / as.matrix(totalFile[filter])
  write.table(sampleCorrected, paste(sampleID, ".", strand,  ".corrected.bins.table.tsv", sep=""), , quote = FALSE, sep ="\t", row.names = TRUE)
  for (l in 1:length(controlFiles)  )
  {
    binFile <- as.matrix(controlFiles[[l]])
    binFile[filter] <- as.matrix(binFile[filter]) / as.matrix(totalFile[filter])
    chiFilesList[[l]] <- binFile
    chiFileNames[l] <- sub(".bins.", ".corrected.bins.", files[l])
  }
  allChiData <- list(chiFilesList, chiFileNames)
  saveRDS(allChiData, paste(sampleID,".", strand, ".controlfiles.corrected.bins.rds", sep=""))
}
#################################################################################################################

#Script
#Stores the command line arguments in a vector 
args<-commandArgs(TRUE)
#Empty matrix designated for total reads per bin for all samples
totals <- matrix(0, nrow=24, ncol=4985)
#Loads the sample 
sample <- read.delim(args[2], header = TRUE, sep = "\t", quote = "\"",
                     dec = ".", fill = TRUE)
setwd(args[3])
#Reads the list of control files
best50list <- read.table(args[4])
#Converts it to a vector
filenames <- as.vector(best50list[,1])
#Selects only the best 50
best50 <- filenames[1:50]
#Gets the filenames of the best 50 control files
files <- GetFiles(controlDir = args[1], best50 = best50, strand = args[5])
#Loads the files
controlFiles <- LoadFiles(files = files, controlDir = args[1])
#Makes the (merged) control file
goodBinFilter<- MakeControlFile(controlFiles, sample = sample)
#Sums the total of control files + the sample
totals <- totals + sample
#Sums the total of autosomal reads of all 50 control files + sample
total <- sum(totals[1:22,])
#Divides the total by number of samples * number of bins
total <- total / (length(controlFiles) * 57614)
#print(goodBinFilter[3,3])
#Filters the samples with the goodBinFilter. Bins that have not passed the overdispersion test will be set to "NaN"
FilterBins(strand = args[5], sample = sample, avgFile = goodBinFilter, chidir = args[3], total = total, sampleID = args[6])
#Quits the R environment
quit(save = "no", status = 0)
