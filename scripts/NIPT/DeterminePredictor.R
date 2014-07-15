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
  for (j in 1:length(files))
  {
    file <- files[[j]]
    total <- sum(file[1:22,])
    
    for (i in 1:22)
    {
      fraction <- sum(file[i,]) / total
      allmatrix[i +  multiplier,j] <<- fraction
    
    }
  }
  
}

GetRowsAndCols <- function(filenamesForward)
{
  cols <- vector(length = 50)
  for (g in 1:length(filenamesForward))
  {
    cols[g] <- sub("^([^.]*).*", "\\1", filenamesForward[g])
    
  }
  rows<- c(paste("Chr", 1:22, "F", sep =""), paste("Chr", 1:22, "R", sep =""))
  rownames(allmatrix) <<- rows
  colnames(allmatrix) <<- cols
}

removeCols <- function(transposedMatrix, chromo)
{
  transposedMatrix[, which(colnames(transposedMatrix)==chromo)] <- NULL
  
  return(transposedMatrix)
}
GetNextPredictor <- function(samples, chrFocusReads, predictors, step)
{
  RSquared <- vector(mode = "numeric")
  if (step == 1)
  {
    for (i in 1:length(samples))
    {
    model <- lm(chrFocusReads ~ transposedSamplesControl[,i])
    RSquared[i] <- summary(model)$adj.r.square
    }
  }
  if (step == 2)
  {
    for (i in 1:length(samples))
    {
      model <- lm(chrFocusReads ~ transposedSamplesControl[, which(colnames(transposedSamplesControl)==predictors[1])] + samples[,i])
      RSquared[i] <- summary(model)$adj.r.square
    }
  }
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
  
  RSquaredSorted <- order(RSquared, decreasing = TRUE)
  chromosomes <- names(samples)
  return(chromosomes[RSquaredSorted[1]])
}

#################################################################################################################

#Script
#Stores the command line arguments in a vector 
args<-commandArgs(TRUE)
print(args)
chromo.focus = as.integer(args[6])

setwd(args[1])
forwardFileList <- readRDS(paste(args[5],".forward.controlfiles.corrected.bins.rds", sep=""))
reverseFileList <- readRDS(paste(args[5],".reverse.controlfiles.corrected.bins.rds", sep=""))

filesForward <- forwardFileList[[1]]
filesReverse <- reverseFileList[[1]]

allmatrix <- matrix(0, nrow=44, ncol=50)

constructMatrix(filesForward, multiplier = 0)
constructMatrix(filesReverse, multiplier = 22)

GetRowsAndCols(forwardFileList[[2]])

setwd(args[1])
write.table(allmatrix,args[4], quote = FALSE, sep ="\t", row.names = TRUE,
            col.names = TRUE)

chrFocusReads <- allmatrix[chromo.focus,] + allmatrix[(chromo.focus + 22),]

output <- as.data.frame(matrix(0, nrow=4, ncol=15))
rownames(output)<- c("Set1", "Set2", "Set3", "Set4")
colnames(output) <- c("Pred1", "Pred2", "Pred3", "Pred4", "RSquared", "AdjRSquared", "FStatistic", 
                        "Sigma", "tStatIntercept", "tStatSlope", "Intercept", "Pred1Slope", "Pred2Slope",
                      "Pred3Slope","Pred4Slope")

transposedSamples <- as.data.frame(t(allmatrix))

transposedSamples <- removeCols(transposedSamples, chromo = paste("Chr", 13, "F", sep=""))
transposedSamples <- removeCols(transposedSamples, chromo = paste("Chr", 13, "R", sep=""))
transposedSamples <- removeCols(transposedSamples, chromo = paste("Chr", 18, "F", sep=""))
transposedSamples <- removeCols(transposedSamples, chromo = paste("Chr", 18, "R", sep=""))
transposedSamples <- removeCols(transposedSamples, chromo = paste("Chr", 21, "F", sep=""))
transposedSamples <- removeCols(transposedSamples, chromo = paste("Chr", 21, "R", sep=""))
transposedSamplesControl <- transposedSamples

for (i in 1:4)
{
predictors <- vector(mode = "numeric")
predictorModel <- list()
  for (j in 1:4)
  {  
    predictors[j] <- GetNextPredictor(transposedSamples, chrFocusReads = chrFocusReads, predictors = predictors, step = j)
    predictorModel[[j]] <- transposedSamplesControl[, which(colnames(transposedSamplesControl)==predictors[j])]
    transposedSamples <- removeCols(transposedSamples, chromo = predictors[j])
  }
model <- lm(chrFocusReads ~ predictorModel[[1]]
            + predictorModel[[2]]
            + predictorModel[[3]]
            + predictorModel[[4]])

  cofs <- round(coef(model),3)
  fstat <-round((summary(model)$fstatistic), 3)
  tvalues <- summary(model)$coefficients
  
  output[i,] <- c(predictors[1], predictors[2], predictors[3], predictors[4], round((summary(model)$r.square) , 3),
                  round((summary(model)$adj.r.square) ,3 ), fstat[1], (summary(model)$sigma), round(tvalues[1,3], 3),
                  round(tvalues[2,3], 3), cofs[1], cofs[2], cofs[3], cofs[4], cofs[5])          
             
}
setwd(args[1])
write.table(output, paste( chromo.focus, args[2],  sep=""), quote = FALSE, sep ="\t", row.names = TRUE,
            col.names = TRUE)

transposedSamples <- as.data.frame(t(allmatrix))
transposedSamples <- removeCols(transposedSamples, chromo = paste("Chr", 13, "F", sep=""))
transposedSamples <- removeCols(transposedSamples, chromo = paste("Chr", 13, "R", sep=""))
transposedSamples <- removeCols(transposedSamples, chromo = paste("Chr", 18, "F", sep=""))
transposedSamples <- removeCols(transposedSamples, chromo = paste("Chr", 18, "R", sep=""))
transposedSamples <- removeCols(transposedSamples, chromo = paste("Chr", 21, "F", sep=""))
transposedSamples <- removeCols(transposedSamples, chromo = paste("Chr", 21, "R", sep=""))
transposedSamplesControl <- transposedSamples

output <- as.data.frame(matrix(0, nrow=4, ncol=17))
rownames(output)<- c("Set1", "Set2", "Set3", "Set4")
colnames(output) <- c("Pred1", "Pred2", "Pred3", "Pred4", "Pred5", "RSquared", "AdjRSquared", "FStatistic", 
                      "Sigma", "tStatIntercept", "tStatSlope", "Pred1Slope", "Pred2Slope",
                      "Pred3Slope","Pred4Slope, Pred5Slope")

for (i in 1:4)
{
  predictors <- vector(mode = "numeric")
  predictorModel <- list()
  for (j in 1:5)
  {  
    predictors[j] <- GetNextPredictor(transposedSamples, chrFocusReads = chrFocusReads, predictors = predictors, step = j)
    predictorModel[[j]] <- transposedSamplesControl[, which(colnames(transposedSamplesControl)==predictors[j])]
    transposedSamples <- removeCols(transposedSamples, chromo = predictors[j])
  }
   model <- lm(chrFocusReads ~ predictorModel[[1]]
              + predictorModel[[2]]
              + predictorModel[[3]]
              + predictorModel[[4]]
              + predictorModel[[5]])
  
  cofs <- round(coef(model),3)
  fstat <-round((summary(model)$fstatistic), 3)
  tvalues <- summary(model)$coefficients
  
  output[i,] <- c(predictors[1], predictors[2], predictors[3], predictors[4], predictors[5], round((summary(model)$r.square) , 3),
                  round((summary(model)$adj.r.square) ,3 ), fstat[1], (summary(model)$sigma), round(tvalues[1,3], 3),
                  round(tvalues[2,3], 3), cofs[1], cofs[2], cofs[3], cofs[4], cofs[5], cofs[6])          
}
setwd(args[1])
write.table(output, paste(chromo.focus, args[3], sep=""), quote = FALSE, sep ="\t", row.names = TRUE,
            col.names = TRUE)

