################################################################################
# This script determines the 50 best matching control samples. Selection of these 
# samples is based on the lowest sum of squares.
# this is calculated by subtracting chromosome fraction of controls by cromosomal fraction of sample.
# The result is squared and summed. The 50 control samples with the lowest result are selected as best control

###############LIBRARIES########
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
  #total reads of the sample on autosomal chromosomes
  totalF <- sum(sampleF[1:22,])
  totalR <- sum(sampleR[1:22,])
  totalSample <- totalF + totalR
  #Vector for storing the chromosomal fractions of the autosomal reads 
  sampleFractions <- vector(mode="numeric")
  #Chromosomal fractions of forward and reverse strands are calculated and stored in the samplefraction Matrix
  for (i in 1:22)
  {
    sampleFractions[i] <- sum(sampleF[i,]) / totalSample
    sampleFractions[i+22] <- sum(sampleR[i,]) / totalSample
  }
  #The result is returned
  return(sampleFractions)
}
#This function calculates the chromosomal fractions of all the controls (forward and reverse apart) and return the 
#chromosomal fractions of all the control files
#files = the control files
#multiplier = This variable makes a distinction between forward and reverse. There are 22 autosomal chromosomes, the fractions 
#are stored in rows in a matrix. For example, chromosome 7F is the 7th row in the matrix, and 7R is the 7+22 row in the matrix
#Totalreads = Total number of reads in the control set
#sampleFractions = The chromosomal fractions of the sample 
#
#returns a matrix with the chromosomal fractions of all the control files
constructMatrix <- function(files, multiplier,totalReads, sampleFractions)
{
  #For each file
  for (j in 1:length(files))
  {
    #Takes the file out of the list as "file"
    file <- files[[j]]
    #Geths the corresponding total
    total <- totalReads[j]
    #For each chromosome the fractions is calculated
    for (i in 1:22)
    {
      fraction <- sum(file[i,]) / total
      
      allmatrix[i +  multiplier,j] <<- (fraction - sampleFractions[i+multiplier]) ^ 2  
    }
  }  
}
#Sums the total autosomal reads in a file and return this number
GetTotalAutosomalReads <- function(file)
{
  return (sum(file[1:22,]))
}

#Sums the total of the autosomal chromosomes minus chromosomes 13, 18 and 21
GetTotals <- function(sample)
{
  return(sum(sample[chromosomes]))
}

############SCRIPT###############
#Stores the command line arguments in a vector 
#args[1] = temp directory where files produced during run of pipeline are stored
#args[2] = sample ID
#args[3] = output, matrix with sum of squares of every chromosome of every control set, 
# and total per control set
# args[4] = output, list of 50 best control sampels
args<-commandArgs(TRUE)
setwd(args[1])
#Reads in the intermediate corrected control files data
forwardFileList <- readRDS(paste(args[2],".forward.controlfiles.intermediate.corrected.bins.rds", sep=""))
reverseFileList <- readRDS(paste(args[2],".reverse.controlfiles.intermediate.corrected.bins.rds", sep=""))
#Gets the files intermediate corrected bin files
filesForward <- forwardFileList[[1]]
filesReverse <- reverseFileList[[1]]
#Gets the total autosomal reads per sample
totalForward <- sapply(filesForward, GetTotalAutosomalReads)
totalReverse <- sapply(filesReverse, GetTotalAutosomalReads)
#Gets the filenames 
colnames <- forwardFileList[[2]]
#Converts the filenames to sample names
colnames <- sapply(colnames, GetCols)
#Gets total reads per control sample
totalReads <- totalForward + totalReverse
#Gets the fraction of the sample
sampleFractions <- GetSampleFracs(args[2])
#Makes an empty matrix for the chromosomal fractions of the control set
allmatrix <- matrix(0, nrow=45, ncol=length(filesForward))
#Sets the colnames 
colnames(allmatrix) <- colnames
#Forward chromosome row numbers except chromosomes 13, 18 and 21
chromosomesF <-c(1:12, 14:17, 19:20, 22) 
#Reverse chromosome row numbers except chromosomes 13, 18 and 21
chromosomesR <-c(23:34,36:39, 41:42, 44)
#Forward and reverse chromosomes combined
chromosomes <- c(chromosomesF, chromosomesR)
#Fills the matrix with chromosomal fractions
constructMatrix(filesForward, multiplier = 0, totalReads, sampleFractions)
constructMatrix(filesReverse, multiplier = 22, totalReads, sampleFractions)
#Gets the square of sums from the control samples
bestNames <- apply(allmatrix, 2, GetTotals)
#rounds of to a maximum of 8 digits
bestNames <- round(bestNames, 8)
#Inserts the sum of squares to the matrix
allmatrix[45,] <- bestNames 
#Sorts the sum of squares 
bestNames <- sort(bestNames)
bestNames <- as.vector(names(bestNames[1:177]))
#Adds the total sum of squares to the last row of the matrix
controlTable <- as.data.frame(allmatrix[45, bestNames])
#Colname for the sum of square tables 
colnames(controlTable) <-  "Sum of Squares"
#Writes the sum of squares table in PDF formar 
pdf("best50Controls.pdf", height=15, width=7)
grid.table(controlTable)
dev.off()
#writes the table with the chromosomal fractions  + sum of squares for the control set
write.table(allmatrix, args[3], quote = FALSE, sep ="\t", row.names = TRUE)
#writes list of best control files
write(bestNames, file = args[4])
