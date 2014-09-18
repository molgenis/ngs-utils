#################################################################################################################
#This script calculates normalized chi square score per bin. If chi square score exceeds 3.5
#bins are corrected by dividing read counts by chi square score / degrees of freedom
#################################################################################################################

#FUNCTIONS#
#Gets the filenames from the binned control files. The input is the directery where the .tsv files live and the strand
#(forward or reverse)
GetFiles <- function(controlDir, best50, strand)
{
  setwd(controlDir)
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
    binfile <- as.matrix(read.delim(files[l], header = TRUE, sep = "\t", quote = "\"",
                          dec = ".", fill = TRUE))
    controlFiles[[l]]<- binfile[1:22,]
    #While file is loaded sum the autosomal reads to total
    totalReadsAllSamples <<- totalReadsAllSamples + sum(binfile)
  }
  return (controlFiles)
}
ChiCorrection  <- function(sample, controlFiles, averageReadsAllSamples)
{
  #List for normalized bin counts. Normalized bin counts are also the observed values for the Chi^2 formula
  normalizedBins <- list()
  #Sum autosomal sample reads
  autosomalReadsSample <- sum(as.matrix(sample))
  #Correction factor to normalize sample bin counts 
  correctionFactor <- averageReadsAllSamples / autosomalReadsSample
  #Normalize sample bin counts
  normalizedSample <- as.matrix(sample) * correctionFactor
  #Copy normalized sample reads, so all controls and samples can be added up to determine the mean 
  normalizedTotal <- normalizedSample
  #For every control sample, same as sample
  for (l in 1:length(controlFiles) )
  {
    #Load control file
    controlSample <- controlFiles[[l]]  
    #Sum autosomal reads of sample
    autosomalReadsControlSample <- sum(controlSample)
    #Get correction factor for normalization by dividing average reads of all samples
    #by reads of this sample
    correctionFactor <- averageReadsAllSamples / autosomalReadsControlSample
    #Multiply bin counts by correction factor
    normalizedControl <- controlSample * correctionFactor
    #Add to total to determine mean later
    normalizedTotal <- normalizedTotal + normalizedControl
    #Store normalized bin counts in list
    normalizedBins[[l]] <- normalizedControl
  }
  #Expected for chi square, expected = mean of observed
  expected <- normalizedTotal / (length(controlFiles) +1 )
  #Chi Square correction for sample
  chiScoreSample <- (expected - normalizedSample) ^ 2 / expected
  #Copy sample Chi square, so total chi square per bin can be summed
  chiScoreTotal <- chiScoreSample
  #For every control sample, same as sample
  for (i in 1:length(normalizedBins))
  {
    correctedObserved <- normalizedBins[[i]]
    chiScoreControlSample <- (expected - correctedObserved) ^ 2 / expected
    chiScoreTotal <- chiScoreTotal + chiScoreControlSample 
  }
  #return total summed chi square scores per bin 
   return (chiScoreTotal)
}

#This functions applies the filter for overdispersed bins. There is nothing returned, the results are written to disk
#in tsv format in an intermediate directory
FilterBins <- function(strand, sample, chiScoreTotal, degreesOfFreedom, chidir, sampleID)
{
  setwd(chidir)
  chiFilesList <- list()
  chiFileNames <- vector(mode = "character")
  #Convert chi squares to a normal distribution score 
  normalizedChiScoreTotal <- (chiScoreTotal - degreesOfFreedom) / (sqrt( 2 * degreesOfFreedom))
  #Total chi score divided by degrees of freedom is the correction factor for the readcounts
  chiCorrectionFactor <- chiScoreTotal / degreesOfFreedom
  #Intermediate matrix that holds the corrected readcounts
  correctedValues <- matrix(0, nrow=22, ncol=4985)
  #Logical matrix, value = true means cell is corrected
  filter <- matrix(FALSE, nrow=22, ncol=4985)
  #if normalized chi square is above 3.5 cells are changed to true
  filter[normalizedChiScoreTotal > 3.5] <- TRUE
  #Corrects cells using the logical filter matrix
  correctedValues[filter] <- as.matrix(sample[filter]) / as.matrix(chiCorrectionFactor[filter])
  #Replaces raw values with corrected values in sample
  sample[filter] <- correctedValues[filter]
  #Writes corrected sample to disk
  write.table(sample, paste(sampleID, ".", strand,  ".corrected.bins.table.tsv", sep=""), , quote = FALSE, sep ="\t", row.names = TRUE)
  #Same procedure as sample, now for every control file
  for (l in 1:length(controlFiles)  )
  {
    binFile <- as.matrix(controlFiles[[l]])
    correctedValues[filter] <- as.matrix(binFile[filter]) / as.matrix(chiCorrectionFactor[filter])
    binFile[filter] <- correctedValues[filter]
  
    chiFilesList[[l]] <- binFile
    #Vector for filenames, later to be used as a index
    chiFileNames[l] <- sub(".bins.", ".corrected.bins.", files[l])
  }
  #Wrapper list for filenames and corrected bin count matrices
  allChiData <- list(chiFilesList, chiFileNames)
  saveRDS(allChiData, paste(sampleID,".", strand, ".controlfiles.corrected.bins.rds", sep=""))
}
#################################################################################################################

#Script
#Stores the command line arguments in a vector 
#args[1] = Directory of control set
#args[2] = binned file of sample
#args[3] = temp directory where files produced during run of pipeline are stored
#args[4] = list of 50 best control samples
#args[5] = strand, either "forward" or "reverse" 
#args[6] = sample ID 
args<-commandArgs(TRUE)
#Empty matrix designated for total reads per bin for all samples

#Loads the sample 
sample <- read.delim(args[2], header = TRUE, sep = "\t", quote = "\"",
                     dec = ".", fill = TRUE)
#Only autosomal reads
sample <- as.matrix(sample[1:22, ])
#Copies autosomal reads sample. Every control sample's reads will later also be added to this
totalReadsAllSamples <- sum(sample)
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
averageReadsAllSamples <- totalReadsAllSamples / (length(controlFiles) + 1)
#Determines the chi square score per bin
chiScoreTotal <- ChiCorrection(sample = sample, controlFiles = controlFiles, averageReadsAllSamples = averageReadsAllSamples)
#Applies correction to bins
FilterBins(strand = args[5], sample = sample, chiScoreTotal = chiScoreTotal , degreesOfFreedom = length(controlFiles), 
           chidir = args[3], sampleID = args[6])
#Quits the R environment
quit(save = "no", status = 0)
