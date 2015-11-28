#################################################################################################################
# This script identifies overdispersed bins. Overdispersed bins are set to 0 (reads) 
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
#Makes a template logical matrix. If a cell = true, bin did not pass overdispersion test.
MakeControlFile <- function(controlFiles, sample)
{
  #Makes an empty matrix. 
  controlFile <- matrix(0, nrow=numberOfRows, ncol=numberOfCols) 
  preCorrectedFiles <- list()
  correctedFiles <- list()
  #sums autosomal reads of sample 
  autosomeReadsSampleReads <- sum(sample[1:22,])
  #autosomal reads of all samples combined
  autosomeReadsTotalReads <- as.numeric(autosomeReadsSampleReads)
  preCorrectedFiles[[1]] <- sample / autosomeReadsSampleReads 
  #for every control file
  for (l in 1:length(controlFiles) )
    {
       file <- controlFiles[[l]]  
       #adds autosomal reads of current control file to total
       autosomeReadsTotalReads <- as.numeric(autosomeReadsTotalReads + sum(file[1:22,]))
       #Determines fractions of every bin of current control sample 
       controlFraction <- file / sum(file[1:22,])
       #stores fractionized current control file in list 
       preCorrectedFiles[[l+1]] <- controlFraction
       
     }
  #Determines average number of reads per sample 
  autosomeReadsAverageReads <- autosomeReadsTotalReads / (length(controlFiles) + 1)
  # For every control file multiply fraction with average number of autosomal reads per sample
  for (i in 1:length(preCorrectedFiles))
  {
    correctedFiles[[i]] <- preCorrectedFiles[[i]] * autosomeReadsAverageReads
  }
  #Empty matrix that will store the average normalized reads per sample 
  avgFile <- matrix(0, nrow=numberOfRows, ncol=numberOfCols)
  #Sum all corrected control files
  for (i in 1:length(correctedFiles))
  {
   avgFile <- avgFile + correctedFiles[[i]]
  }
  #Divide by number of control files 
  avgFile <- avgFile / length(correctedFiles)
  #Normalize by dividing corrected number of reads per bin by average number of reads per bin
  endFile <- avgFile / (autosomeReadsAverageReads / numberOfBins)
  #Logical matrix, all values are FALSE
  goodBinFilter <- matrix(FALSE, nrow=numberOfRows, ncol=numberOfCols) 
  #If endFile scores are above 1.5 set corresponding logical value to TRUE
  goodBinFilter[endFile < 1.5] <- TRUE
  return(goodBinFilter)
}
#This functions applies the filter for overdispersed bins. There is nothing returned, the results are written to disk
#in tsv format in an intermediate directory
FilterBins <- function(strand, sample, filter, workdir, filename)
  {
  #List for storing corrected samples
  intermediateFilesList <- list()
  #Vector where samplenames of corrected files are being stored 
  intermediateFileNames <- vector(mode="character")
  #Applying the logical vector, if value = TRUE corresponding bins are set to 0
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
  #Wrapper list, contains the list of corrected samples and vector containing the filenames of sample 
  allIntermediateData <- list(intermediateFilesList, intermediateFileNames)
  saveRDS(allIntermediateData, paste(filename,".", strand, ".controlfiles.intermediate.corrected.bins.rds", sep=""))
   }
#################################################################################################################

                                            #Script
#Stores the command line arguments in a vector 
#args[1] = Directory of control set
#args[2] = sample filename
#args[3] = temp directory where files produced during run of pipeline are stored
#args[4] = sample ID 
#args[5] = strand, either "forward" or "reverse" 
numberOfCols <- 4985
numberOfRows <- 24
numberOfBins <- 57614
args<-commandArgs(TRUE)
#Loads the sample 
sample <- read.delim(args[2], header = TRUE, sep = "\t", quote = "\"",
                     dec = ".", fill = TRUE)
#Gets the control filenames 
files <- GetFiles(controlDir = args[1], strand = args[5], sampleID = args[4])
#Loads the files control files
controlFiles <- LoadFiles(files = files)
#Determines which bins have chi score above the threshold
goodBinFilter<- MakeControlFile(controlFiles, sample = sample)
#Filters the samples with the goodBinFilter. Bins that have not passed the overdispersion test will be set to 0
FilterBins(strand = args[5], sample = sample, filter = goodBinFilter, workdir = args[3], filename = args[4])
#Quits the R environment
quit(save = "no", status = 0)
