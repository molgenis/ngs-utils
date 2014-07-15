#################################################################################################################

                                          #FUNCTIONS#
#Gets the filenames from the binned control files. The input is the directery where the .tsv files live and the strand
#(forward or reverse)
GetFiles <- function(controlDir, strand, sampleID)
{
  setwd(controlDir)
  files <- list.files( pattern = paste("*", strand,".bins.table.tsv", sep = ""))
  sampleID <- paste(sampleID,".", strand, ".bins.table.tsv", sep ="")
  files <- files[!files %in% sampleID]
  return (files)
}
#Loads the .tsv binned files into a list and returns this list (controlFiles)
LoadFiles <- function(files)
{
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
  controlFile <- matrix(0, nrow=24, ncol=4985) 
  preCorrectedFiles <- list()
  correctedFiles <- list()
  autosomeReadsSampleReads <- sum(sample[1:22,])
  autosomeReadsTotalReads <- as.numeric(autosomeReadsSampleReads)
  preCorrectedFiles[[1]] <- sample / autosomeReadsSampleReads 
  #In this loop all the control files are summed into a single file
  for (l in 1:length(controlFiles) )
    {
       file <- controlFiles[[l]]  
       autosomeReadsTotalReads <- as.numeric(autosomeReadsTotalReads + sum(file[1:22,]))
       controlFraction <- file / sum(file[1:22,])
       preCorrectedFiles[[l+1]] <- controlFraction
       
     }
  autosomeReadsTotalReads <- autosomeReadsTotalReads / (length(controlFiles) + 1)
  for (i in 1:length(preCorrectedFiles))
  {
    test <- preCorrectedFiles[[i]] * autosomeReadsTotalReads
    correctedFiles[[i]] <- test
  }
  avgFile <- matrix(0, nrow=24, ncol=4985)
  for (i in 1:length(correctedFiles))
  {
   avgFile <- avgFile + correctedFiles[[i]]
  }
  avgFile <- avgFile / length(correctedFiles)
  endFile <- avgFile / (autosomeReadsTotalReads / 62217)
  test <- which.max(endFile[1,])
  goodBinFilter <- matrix(FALSE, nrow=24, ncol=4985) 
  goodBinFilter[endFile < 1.5] <- TRUE
  return(goodBinFilter)
}
#This functions applies the filter for overdispersed bins. There is nothing returned, the results are written to disk
#in tsv format in an intermediate directory
FilterBins <- function(strand, sample, filter, workdir, filename)
  {
  intermediateFilesList <- list()
  intermediateFileNames <- vector(mode="character")
  sample[!filter] <- 0
  setwd(workdir)
  write.table(sample, paste(args[4], ".", strand,  ".intermediate.corrected.bins.table.tsv", sep =""), quote = FALSE, sep ="\t", row.names = TRUE)
  #Iterates over the files, applies the filter and writes the output to disk
  for (l in 1:length(controlFiles)  )
    {
        binFile <- controlFiles[[l]]
        binFile[!filter] <- 0
        intermediateFilesList[[l]] <- binFile
        intermediateFileNames[l] <- sub(".bins.", ".intermediate.corrected.bins.", files[l])
  }
  allIntermediateData <- list(intermediateFilesList, intermediateFileNames)
  saveRDS(allIntermediateData, paste(filename,".", strand, ".controlfiles.intermediate.corrected.bins.rds", sep=""))
   }
#################################################################################################################

                                            #Script
#Stores the command line arguments in a vector 
args<-commandArgs(TRUE)
#Loads the sample 
sample <- read.delim(args[2], header = TRUE, sep = "\t", quote = "\"",
                     dec = ".", fill = TRUE)
#Gets the filenames from the control files
files <- GetFiles(controlDir = args[1], strand = args[5], sampleID = args[4])
#Loads the files
controlFiles <- LoadFiles(files = files)
#Makes the (merged) control file
goodBinFilter<- MakeControlFile(controlFiles, sample = sample)
#Filters the samples with the goodBinFilter. Bins that have not passed the overdispersion test will be set to "NaN"
FilterBins(strand = args[5], sample = sample, filter = goodBinFilter, workdir = args[3], filename = args[4])
#Quits the R environment
quit(save = "no", status = 0)
