library(gridExtra)


##########FUNCTIONS#########

#This function makes a list of the intermediate corrected bin files produced in the step 'InterChiCorrection' 
#controlDir = The directory where the intermediate files are stored
#strand = The strand ( Forward or Reverse)
#sampleID = the sample ID, the sample is removed from the vector of possible control files
#
#returns a vector containing the files to be read
GetFiles <- function(controlDir, strand, sampleID)
{
  setwd(controlDir)
  files <- list.files( pattern = paste(strand, ".intermediate.corrected.bins.table.tsv", sep = ""))
  sampleID <- paste(sampleID,".", strand, ".intermediate.corrected.bins.table.tsv", sep ="")
  #Removes the sample from the possible control files
  files <- files[!files %in% sampleID]
  return (files)
}


#Loads the .tsv binned files into a list and returns this list (controlFiles)
#files = vector containing the filenames of all possible control files
#
#returns a list of bin files
LoadFiles <- function(files)
{
  #Empty list is created that will store the intermediate corrected bin files
  controlFiles <- list()
  #Eevery file is read and stored in the controlFiles list
  for (l in 1:length(files) ){
    binfile <- read.delim(files[l], header = TRUE, sep = "\t", quote = "\"",
                          dec = ".", fill = TRUE)
    #na's (not availables) are converted to 0 to make calulations on the bin files possibke
    binfile[is.na(binfile)] <- 0
    #stored the bin files in the list
    controlFiles[[l]]<- binfile
    
  }
  return (controlFiles)
}

#This function uses a regular expression to convert the file names into sample ID's
#file = a control filename
#
#returns a sample ID that will be used as a colname
GetCols <- function(file)
{
  colname <- sub("^([^.]*).*", "\\1", file)
  return(colname)
}
#This function determines the fraction of the chromosomes of the sample
#sample = filename of the sample
#
#retuns a matrix with the fraction of the chromosomes
GetSampleFracs <- function(sample)
{
  #Reads in the sample intermediate binfile
  sampleF <- read.delim(paste(sample, ".forward.intermediate.corrected.bins.table.tsv", sep = ""), header = TRUE, sep = "\t", quote = "\"",
                        dec = ".", fill = TRUE)
  
  sampleR <- read.delim(paste(sample, ".reverse.intermediate.corrected.bins.table.tsv", sep = ""), header = TRUE, sep = "\t", quote = "\"",
                        dec = ".", fill = TRUE)
  #NA's are converted to 0 for calculations
  sampleF[is.na(sampleF)] <- 0
  sampleR[is.na(sampleR)] <- 0
  #total reads 
  totalF <- sum(sampleF[1:22,])
  totalR <- sum(sampleR[1:22,])
  totalSample <- totalF + totalR
  sampleFractions <- vector(mode="numeric")
  for (i in 1:22)
  {
    sampleFractions[i] <- sum(sampleF[i,]) / totalSample
    sampleFractions[i+22] <- sum(sampleR[i,]) / totalSample
  }
  return(sampleFractions)
}

constructMatrix <- function(files, multiplier,totalReads, sampleFractions, chromosomes)
{
  for (j in 1:length(files))
  {
    file <- files[[j]]
    total <- totalReads[j]
    print(total)
    for (i in 1:22)
    {
      fraction <- sum(file[i,]) / total
      
      allmatrix[i +  multiplier,j] <<- (fraction - sampleFractions[i+multiplier]) ^ 2
      
    }
  }
  
}

GetTotalAutosomalReads <- function(file)
{
  return (sum(file[1:22,]))
}
GetTotals <- function(sample)
{
  return(sum(sample[chromosomes]))
}

############SCRIPT###############

args<-commandArgs(TRUE)

setwd(args[1])

forwardFileList <- readRDS(paste(args[2],".forward.controlfiles.intermediate.corrected.bins.rds", sep=""))
reverseFileList <- readRDS(paste(args[2],".reverse.controlfiles.intermediate.corrected.bins.rds", sep=""))

filesForward <- forwardFileList[[1]]
filesReverse <- reverseFileList[[1]]

totalForward <- sapply(filesForward, GetTotalAutosomalReads)
totalReverse <- sapply(filesReverse, GetTotalAutosomalReads)

colnames <- forwardFileList[[2]]

colnames <- sapply(colnames, GetCols)

totalReads <- totalForward + totalReverse

sampleFractions <- GetSampleFracs(args[2])

allmatrix <- matrix(0, nrow=45, ncol=length(filesForward))

colnames(allmatrix) <- colnames

chromosomesF <-c(1:12, 14:17, 19:20, 22) 
chromosomesR <-c(23:34,36:39, 41:42, 44)
chromosomes <- c(chromosomesF, chromosomesR)

constructMatrix(filesForward, multiplier = 0, totalReads, sampleFractions)
constructMatrix(filesReverse, multiplier = 22, totalReads, sampleFractions)

bestNames <- apply(allmatrix, 2, GetTotals)
bestNames <- round(bestNames, 8)
allmatrix[45,] <- bestNames 
bestNames <- sort(bestNames)
bestNames <- as.vector(names(bestNames[1:50]))

controlTable <- as.data.frame(allmatrix[45, bestNames])


colnames(controlTable) <-  "Sum of Squares"

pdf("best50Controls.pdf", height=15, width=7)
grid.table(controlTable)
dev.off()


write.table(allmatrix, args[3], quote = FALSE, sep ="\t", row.names = TRUE)
write(bestNames, file = args[4])
