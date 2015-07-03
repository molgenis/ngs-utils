#####################################################################################################
# This script selects four sets of four best predictors each. Selection is based on adjusted R score
#



#################################################################################################################
                                #FUNCTIONS#

#This function calculates the chromosomal fractions of all the controls (forward and reverse apart) and return the 
#chromosomal fractions of all the control files
#files = the control files
#multiplier = This variable makes a distinction between forward and reverse. There are 22 autosomal chromosomes, the fractions 
#are stored in rows in a matrix. For example, chromosome 7F is the 7th row in the matrix, and 7R is the 7+22 row in the matrix
#
#returns a matrix with the chromosomal fractions of all the control files
constructMatrix <- function(files, multiplier)
{
  autoChromosomes <-c(1:12, 14:17, 19:20, 22) 
  for (j in 1:length(files))
  {
    file <- files[[j]]
    total <- sum(file[autoChromosomes,])
    
    for (i in 1:22)
    {
      fraction <- sum(file[i,]) / total
      allmatrix[i +  multiplier,j] <<- fraction
    
    }
  }
  
}
#This functions sets the row and column names of the dataframe. The input are the filenames.
#A regular expression is used to convert these to filenames to sampelnames
GetRowsAndCols <- function(filenamesForward)
{
  cols <- vector(length = length(filenamesForward))
  for (g in 1:length(filenamesForward))
  {
    cols[g] <- sub("^([^.]*).*", "\\1", filenamesForward[g])
    
  }
  rows<- c(paste("Chr", 1:22, "F", sep =""), paste("Chr", 1:22, "R", sep =""))
  rownames(allmatrix) <<- rows
  colnames(allmatrix) <<- cols
}
#This function removes a column from the dataframe. For instance, for removing a selected predictor
removeCols <- function(transposedMatrix, chromo)
{
  transposedMatrix[, which(colnames(transposedMatrix)==chromo)] <- NULL
  
  return(transposedMatrix)
}
#This function selects the next predictor. 
GetNextPredictor <- function(samples, chrFocusReads, predictors, step)
{
  #Vector to hold r square for predictions
  RSquared <- vector(mode = "numeric")
  #if step is 1, meaning first predictor is selected 
  if (step == 1)
  {
    #All possible candidate predictors are used and their adjusted r square value stored in RSquared vector
    for (i in 1:length(samples))
    {
    model <- lm(chrFocusReads ~ samples[,i])
    RSquared[i] <- summary(model)$adj.r.square
    }
  }
  #same as step1, but now for second predictor
  if (step == 2)
  {
    for (i in 1:length(samples))
    {
      model <- lm(chrFocusReads ~ transposedSamplesControl[, which(colnames(transposedSamplesControl)==predictors[1])] + samples[,i])
      RSquared[i] <- summary(model)$adj.r.square
    }
  }
  #same as step1, but now for third predictor
  if (step == 3)
  {
    for (i in 1:length(samples))
    {
      model <- lm(chrFocusReads ~ transposedSamplesControl[, which(colnames(transposedSamplesControl)==predictors[1])]
                  + transposedSamplesControl[, which(colnames(transposedSamplesControl)==predictors[2])]
                                             + samples[,i])
      RSquared[i] <- summary(model)$adj.r.square
    }
  }
  if (step == 4)
  {
    for (i in 1:length(samples))
    {
      model <- lm(chrFocusReads ~ transposedSamplesControl[, which(colnames(transposedSamplesControl)==predictors[1])]
                  + transposedSamplesControl[, which(colnames(transposedSamplesControl)==predictors[2])]
                  + transposedSamplesControl[, which(colnames(transposedSamplesControl)==predictors[3])]
                  + samples[,i])
      RSquared[i] <- summary(model)$adj.r.square
    }
  }
  
  if (step == 5)
  {
    for (i in 1:length(samples))
    {
      model <- lm(chrFocusReads ~ transposedSamplesControl[, which(colnames(transposedSamplesControl)==predictors[1])]
                  + transposedSamplesControl[, which(colnames(transposedSamplesControl)==predictors[2])]
                  + transposedSamplesControl[, which(colnames(transposedSamplesControl)==predictors[3])]
                  + transposedSamplesControl[, which(colnames(transposedSamplesControl)==predictors[4])]
                  + samples[,i])
      RSquared[i] <- summary(model)$adj.r.square
    }
  }
  #Orders the adjusted r squared values by index in decreasing order
  RSquaredOrdered <- order(RSquared, decreasing = TRUE)
  #Gets all chromosomes remaining in the candidate predictors (for instance, Chr1F, Chr2F etc)
  chromosomes <- names(samples)
  #Returns 
  return(chromosomes[RSquaredOrdered[1]])
}

#################################################################################################################

#Script
#Stores the command line arguments in a vector 
#args[1] = temp directory where files produced during run of pipeline are stored
#args[2] = output, table with predictor sets and relevant statistics
#args[3] = deprecated, cleaned up soon
#args[4] = output, table with chromosomal fraction of all 50 control samples
#args[5] = sample ID
#args[6] = chromosome of focus, either 13, 18 or 21
args<-commandArgs(TRUE)
#Gets the chromosome (13, 18 or 21)
chromo.focus = as.integer(args[6])
#Sets the workdir where the Chi2 corrected files live
setwd(args[1])
#Reads the Chi2 corrected files and filelist
forwardFileList <- readRDS(paste(args[5],".forward.controlfiles.corrected.bins.rds", sep=""))
reverseFileList <- readRDS(paste(args[5],".reverse.controlfiles.corrected.bins.rds", sep=""))
#Gets the Chi2 corrected files
filesForward <- forwardFileList[[1]]
filesReverse <- reverseFileList[[1]]
#Makes and empty matrix
allmatrix <- matrix(0, nrow=44, ncol=length(filesForward))
#Gets the filenames
GetRowsAndCols(forwardFileList[[2]])
#Fills the newly constructed matrix with chromosomal fractions
constructMatrix(filesForward, multiplier = 0)
constructMatrix(filesReverse, multiplier = 22)
#Writes the fraction table to the temporary directory
write.table(allmatrix,args[4], quote = FALSE, sep ="\t", row.names = TRUE,
            col.names = TRUE)
#Gets the reads for each of sample of the chromosome (13 18 or 21) and adds up forward and reverse
chrFocusReads <- allmatrix[chromo.focus,] + allmatrix[(chromo.focus + 22),]
#Makes an empty dataframe to store the predictors
output <- as.data.frame(matrix(0, nrow=4, ncol=15))
rownames(output)<- c("Set1", "Set2", "Set3", "Set4")
colnames(output) <- c("Pred1", "Pred2", "Pred3", "Pred4", "RSquared", "AdjRSquared", "FStatistic", 
                        "Sigma", "tStatIntercept", "tStatSlope", "Intercept", "Pred1Slope", "Pred2Slope",
                      "Pred3Slope","Pred4Slope")
#transposes allmatrix for easy column acces
transposedSamples <- as.data.frame(t(allmatrix))
#removes chromosomes 13, 18 and 21 from possible predictors
transposedSamples <- removeCols(transposedSamples, chromo = paste("Chr", 13, "F", sep=""))
transposedSamples <- removeCols(transposedSamples, chromo = paste("Chr", 13, "R", sep=""))
transposedSamples <- removeCols(transposedSamples, chromo = paste("Chr", 18, "F", sep=""))
transposedSamples <- removeCols(transposedSamples, chromo = paste("Chr", 18, "R", sep=""))
transposedSamples <- removeCols(transposedSamples, chromo = paste("Chr", 21, "F", sep=""))
transposedSamples <- removeCols(transposedSamples, chromo = paste("Chr", 21, "R", sep=""))
#copies the transposedsamples data frame
transposedSamplesControl <- transposedSamples
#this loops runs for 4 times, for 4 sets of predictors
for (i in 1:4)
{
  #Vector to hold predictors 
  predictors <- vector(mode = "numeric")
  #List to hold columns from the transposedsample dataframe
  predictorModel <- list()
    #This loop runs 4 times, for 4 predictors per set 
    for (j in 1:4)
    {  
      #Selects the next predictor
      predictors[j] <- GetNextPredictor(transposedSamples, chrFocusReads = chrFocusReads, predictors = predictors, step = j)
      #Stores the column in a list to build the model later
      predictorModel[[j]] <- transposedSamplesControl[, which(colnames(transposedSamplesControl)==predictors[j])]
      #Removes a selected predictor from the data frame 
      transposedSamples <- removeCols(transposedSamples, chromo = predictors[j])
    }
  #Builds the model
  model <- lm(chrFocusReads ~ predictorModel[[1]]
            + predictorModel[[2]]
            + predictorModel[[3]]
            + predictorModel[[4]])
  #coefficients are extracted
  cofs <- round(coef(model),3)
  #f statistics are extracted
  fstat <-round((summary(model)$fstatistic), 3)
  tvalues <- summary(model)$coefficients
  #Every iteration a row in the output data frame is added
  output[i,] <- c(predictors[1], predictors[2], predictors[3], predictors[4], round((summary(model)$r.square) , 3),
                  round((summary(model)$adj.r.square) ,3 ), fstat[1], (summary(model)$sigma), round(tvalues[1,3], 3),
                  round(tvalues[2,3], 3), cofs[1], cofs[2], cofs[3], cofs[4], cofs[5])          
             
}
#Writes the table with predictors and statistics to disk
write.table(output, paste( chromo.focus, args[2],  sep=""), quote = FALSE, sep ="\t", row.names = TRUE,
            col.names = TRUE)

