#Libs for nice pdf tables
library(gridExtra)
#Function that removes a column from a matrix
removeCols <- function(transposedMatrix, chromo)
{
  transposedMatrix[, which(colnames(transposedMatrix)==chromo)] <- NULL
  
  return(transposedMatrix)
}
#Loads the sample and determines the fractions of the chromosomes
GetSampleFracs <- function(sample)
{
  #Row numbers of the autosomal chromosomes without 13, 18 and 21
  autoChromosomes <-c(1:12, 14:17, 19:20, 22) 
  #Loads the sample bin counts
  sampleF <- read.delim(paste(sample, ".forward.corrected.bins.table.tsv", sep = ""), header = TRUE, sep = "\t", quote = "\"",
                        dec = ".", fill = TRUE)
  
  sampleR <- read.delim(paste(sample, ".reverse.corrected.bins.table.tsv", sep = ""), header = TRUE, sep = "\t", quote = "\"",
                        dec = ".", fill = TRUE)
  #Sums the total of autosomal reads without 13, 18 and 21, forward and reverse aparts
  totalF <- sum(sampleF[autoChromosomes,])
  totalR <- sum(sampleR[autoChromosomes,])
  #Sums forward and reverse reads
  totalReadsChrFocus <<- sum(sampleF[chromo.focus,]) + sum(sampleR[chromo.focus, ])
  #Empty matrix to hold the chromosomal fractions
  sampleFracs <- matrix(0, nrow=44, ncol=1)
  #For each row forward reads the fraction is calculated
  for (i in 1:22)
  {
    fraction <- sum(sampleF[i,]) / totalF
    sampleFracs[i,1] <- fraction
    
  }
  #same, but now for reverse reads. stored in the same matrix
  for (i in 1:22)
  {
    fraction <- sum(sampleR[i,]) / totalR
    sampleFracs[i + 22,1] <- fraction
    
  }
  #Proper row names
  rows<- c(paste("Chr", 1:22, "F", sep =""), paste("Chr", 1:22, "R", sep =""))
  rownames(sampleFracs) <- rows
  #Transposes the matrix and converts to dataframe
  sampleFracs <- as.data.frame(t(sampleFracs))
  return(sampleFracs)
}
#Functions that predicts chromosomal fractions using the model
predictSamples <- function(data)
{
 predict(model1, data)
}
#Builds the model
GetModel <- function(setNumber)
{
  pred1 <- transposedSamplesControl[, which(colnames(transposedSamplesControl)==predictors[setNumber,1])]
  pred2 <- transposedSamplesControl[, which(colnames(transposedSamplesControl)==predictors[setNumber,2])]
  pred3 <- transposedSamplesControl[, which(colnames(transposedSamplesControl)==predictors[setNumber,3])]
  pred4 <- transposedSamplesControl[, which(colnames(transposedSamplesControl)==predictors[setNumber,4])]
  return (model1 <- lm(chrFocusReads ~ pred1 + pred2 + pred3 + pred4))
}
#Sets up the data frame to be used for predicting 
GetPredictDataFrame <- function(setNumber, sourceFrame)
{
  return(data.frame(pred1 = sourceFrame[, which(colnames(sourceFrame)==predictors[setNumber,1])],
                    pred2 = sourceFrame[, which(colnames(sourceFrame)==predictors[setNumber,2])],
                    pred3 = sourceFrame[, which(colnames(sourceFrame)==predictors[setNumber,3])],
                    pred4 = sourceFrame[, which(colnames(sourceFrame)==predictors[setNumber,4])])) 
}
#Reads the .csv files contining the predictors
GetPredictorSets <- function(file)
{
  return(read.csv(file, header = TRUE, sep = "\t",
           dec = ".", fill = TRUE, comment.char = ""))
}
#Gets the mean of a row in a matrix
GetMean <- function(matrixRow, sampleRow)
{
  return(mean(as.matrix(matrixRow)))
}
#Gets the SD of a row in matrix
GetSD <- function(matrixRow)
{
  return(sd(as.matrix(matrixRow)))
}
#Make generic z score table 
MakeTable <- function(colls)
{
  table <- matrix(0, nrow=4, ncol=length(colls)) 
  rowws <- c("Set 1", "Set 2", "Set 3", "Set 4")
  colnames(table) <- colls
  rownames(table) <- rowws
  return (table)
}
#Makes table to store the regression factors
GetRegressionTable <- function(predictor)
{
  #Data frame as matrix
  extractMatrix <- as.matrix(predictor)
  #Empty matrix regression matrix
  regressionTable <- matrix("", nrow=4, ncol=2) 
  #sets col and row names
  cols <- c("Predictor Sets", "Regression Factor")
  rows <- c("Set 1", "Set 2", "Set 3", "Set 4")
  colnames(regressionTable) <- cols
  rownames(regressionTable) <- rows
  #for each predictor set
  for (i in 1:4)
    {
    #Gets predictors
    ps <- as.character(as.vector(extractMatrix[i,c(1:4)] ))
    #Gets regression factors
    regressionFactors <- as.vector(extractMatrix[i,12:15])
    #pastes regression factor to one string
    regressionFactors <- paste(regressionFactors, ps, sep ="")
    #puts a "+" sign between factors
    regressionFactors <- paste(regressionFactors, collapse=" + ")
    #Concatenates chromosomes to factors
    ps <- paste(ps, collapse=" ")
    #Gets the intercept
    intercept <- as.vector(extractMatrix[i,11])
    #Puts a "+" sign after the intercept
    intercept <- paste(intercept, " + ")
    #Concatenates the intercept and the regression tables
    regressionFactors <- paste(intercept, regressionFactors, sep="")
    #Stores result in regression table
    regressionTable[i,] <- c(ps, regressionFactors)
  }
  return(regressionTable)
}
#Converts the regression table to grob
PrepareRegressionTableForPrint <- function(regressionTable, chromo.focus)
{
  d <- head(regressionTable)
  table <- tableGrob(d)
  h <- grobHeight(table)
  w <- grobWidth(table)
  title <- textGrob(paste("Chromosome " , chromo.focus, sep = " "), y=unit(0.5,"npc") + 0.5*h, 
                    x=unit(0.35,"npc") + 0.5*h,
                    vjust=0, gp=gpar(fontsize=24))
  gt <- gTree(children=gList(table, title))
  return(gt)
}
#This function calculates the NCV (zscores)
DetermineNCV <- function(totalReadsChrFocus)
{
  #theoretical variation coefficent
  cv <<- (1.15 * (1 / sqrt(totalReadsChrFocus)))
  #if practical vc is lower than theoretical vc use theoretical
  if(sd(normControls) < (1.15 * (1 / sqrt(totalReadsChrFocus))))
  {
    Zvalues <- (normControls - mean(normControls)) / (1.15*(1 / sqrt(totalReadsChrFocus)))
    Zvalues[1:length(Zvalues)+1] <- (normSample - mean(normControls)) / (1.15*(1 / sqrt(totalReadsChrFocus)))
    corcof <<- "Yes"
  }
  # if not
  else
  {
   Zvalues <- (normControls - mean(normControls)) / sd(normControls)
   Zvalues[1:length(Zvalues)+1] <- (normSample - mean(normControls)) / sd(normControls)
   corcof <<- "No"
  }
  return (Zvalues)
}
#Writes pdf plots
WritePSSetsPDFs <- function(Zvalues, i)
{
  setwd(args[5])
  pdf(paste(args[3],"_Chromosome_", chromo.focus, "_Predictorset_4Predictors_Barplot_", i, ".pdf", sep = ""))
  barplot(Zvalues, col = ifelse(Zvalues < 6,'green','red'), axisnames = FALSE, 
          main = paste("Sample NCV = ", round(Zvalues[1:length(Zvalues)], digits = 2), sep=""), xlab = corcof)
  dev.off()
  pdf(paste(args[3],"_Chromosome_", chromo.focus, "_Predictorset_4Predictors_Boxplot_", i, ".pdf", sep = ""))
  boxplot(Zvalues[1:length(Zvalues) -1], col = ifelse(Zvalues < 6,'green','red'), main = paste("Mean control group = ", controlAv))
  dev.off()
  d <- density(Zvalues[1:length(Zvalues) -1])
  pdf(paste(args[3], "_Chromosome_", chromo.focus,"_Predictorset_4Predictors_Density_", i, ".pdf", sep = ""))
  plot(d)
  dev.off()
}
#################################################################################################################

#Script
#Stores the command line arguments in a vector 
args<-commandArgs(TRUE)

#Gets the chromosome (13, 18 or 21) as a numeric
chromo.focus = as.numeric(args[6])
colsNcvTable <- c("Sample PS", "Observed Sample", "Observed / Predicted", "Mean Control", "SD Control", "Practical VC",
                  "Reads", "Theoretical VC", "Theoretical Used", "P-Value Shapiro", "Zscore")
colsDiagnosticTable <- c("VC", "Normal Distributed", "Zscore")
#Reads the chromosomal fraction table
setwd(args[5])
allmatrix <- read.delim(args[1], header = TRUE, sep = "\t", quote = "\"",
                        dec = ".", fill = TRUE)

transposedSamples <- as.data.frame(t(allmatrix))
chrFocusF <- paste("Chr", chromo.focus, "F", sep="")
chrFocusR <- paste("Chr", chromo.focus, "R", sep="")
chrFocusReads <- transposedSamples[[chrFocusF]] + transposedSamples[[chrFocusR]]
#Removes the possible trisomy chromosomes from data
transposedSamples <- removeCols(transposedSamples, chromo = paste("Chr", 13, "F", sep=""))
transposedSamples <- removeCols(transposedSamples, chromo = paste("Chr", 13, "R", sep=""))
transposedSamples <- removeCols(transposedSamples, chromo = paste("Chr", 18, "F", sep=""))
transposedSamples <- removeCols(transposedSamples, chromo = paste("Chr", 18, "R", sep=""))
transposedSamples <- removeCols(transposedSamples, chromo = paste("Chr", 21, "F", sep=""))
transposedSamples <- removeCols(transposedSamples, chromo = paste("Chr", 21, "R", sep=""))
transposedSamplesControl <- transposedSamples
#Sets the filename for the file containing the predictor data 
filename <- paste(chromo.focus, args[2], sep="")
#Loads the predictor data
predictors <- GetPredictorSets(filename)
#Prepares the regression table, later printed as PDF
regressionTable <- GetRegressionTable(predictors)
#empty vector for focus reads
totalReadsChrFocus <- vector(mode = "numeric")
#Gets the chromosomal fractions of the sample
sampleFracs <- GetSampleFracs(args[3])
#Makes empty NCV table
NCVTable <- MakeTable(colsNcvTable)
# Makes empty diagnostics table
DiagnosticsTable <- MakeTable(colsDiagnosticTable)
for (i in 1:4)
{
  #Initializes cv value. a high value is used so if it goes wrong it will be detected
  cv <- 1  
  #Get model
  moden <- GetModel(i)
  
  corcof <- "Undetermined"
  #Makes data frame for prediction
  predictSample <- GetPredictDataFrame(i, sampleFracs)
  #Predicts the fraction of the sample
  predictedSample <-  predict(moden, predictSample)
  #Sums fractions of chromosome of focus (either 13, 18 or 21)
  observedSample <- sampleFracs[[chrFocusF]] + sampleFracs[[chrFocusR]]
  #Normalize chromosome of focus (observed fraction / predicted fraction)
  normSample <- (sampleFracs[[chrFocusF]] + sampleFracs[[chrFocusR]]) / predictedSample
  #Predict control group
  predictControl <- GetPredictDataFrame(i, transposedSamplesControl)
  #Sets row names
  rownames(predictControl) <- rownames(transposedSamplesControl)
  #predict control group fractions
  predictedControls <- predict(moden, predictControl)
  #Normalized control fractions (observed fraction / predicted fraction)
  normControls <- (chrFocusReads) / predictedControls
  #Determines Z scores
  Zvalues <- DetermineNCV(totalReadsChrFocus)
  #mean Z score control group
  controlAv <- round(mean(Zvalues[1:length(Zvalues) -1]), digits = 2)
  #Add Z score 
  predictors$ncv[i] <- round(Zvalues[length(Zvalues)], 3)
  #Add row to NCVtable
  NCVTable[i,] <- c(round(predictedSample,5), round(observedSample,5), round(observedSample / predictedSample ,5), round(mean(normControls),5), round(sd(normControls),6), 
                  round(sd(normControls) / mean(normControls),5), round(totalReadsChrFocus,0),  round(cv,4), corcof, round(shapiro.test(Zvalues[1:length(Zvalues) -1])$p.value,3), round(Zvalues[length(Zvalues)], 3))
  normalDistributed <- "Yes"
  #If normal group is normally distributed (tests with Shapiro - Wilk test)
  if(shapiro.test(Zvalues[1:length(Zvalues ) -1])$p.value < 0.05)
    {
      normalDistributed <- "No"
  }
  # practical vc
  cvUsed <- round(sd(normControls) / mean(normControls),5) * 100
  #If theoretical vc is used change vc value to theoretical
  if(corcof == "Yes")
  {
    cvUsed <- round(cv,4) * 100
  }
  #Adds vc to table
  DiagnosticsTable[i,] <- c(cvUsed, normalDistributed, round(Zvalues[length(Zvalues)], 3)) 
  #Writes PDFS
  WritePSSetsPDFs(Zvalues, i)
}
#Writes .csv tables
write.table(NCVTable,paste(args[3], "_Chromosome", chromo.focus, "_NCVTable.csv", sep =""), quote = FALSE, sep ="\t", row.names = TRUE,
            col.names = TRUE)
write.table(regressionTable,paste(args[3], "_Chromosome", chromo.focus, "_RegressionTable.csv", sep =""), quote = FALSE, sep ="\t", row.names = TRUE,
            col.names = TRUE)
write.table(t(sampleFracs),paste(args[3], "_samplefractions.csv", sep =""), quote = FALSE, sep ="\t", row.names = TRUE,
            col.names = TRUE)
write.table(DiagnosticsTable,paste(args[3], "_Chromosome", chromo.focus, "_DiagnostiekOutput.csv", sep =""), quote = FALSE, sep =",", 
            col.names = TRUE, row.names = FALSE)
write.table(predictors,paste(chromo.focus, args[4], sep =""), quote = FALSE, sep ="\t", row.names = TRUE,
            col.names = TRUE)
#Creates NCV table
nv <- PrepareRegressionTableForPrint(NCVTable, chromo.focus)
pdf(paste(args[3], "_Chromosome", chromo.focus, "_NCVTable.pdf", sep =""), height=5, width=14)
grid.draw(nv)
dev.off()
#Creates regression table pdf
gt <- PrepareRegressionTableForPrint(regressionTable, chromo.focus)
pdf(paste(args[3], "_Chromosome", chromo.focus, "_RegressionTable.pdf", sep=""), height=5, width=10)
grid.draw(gt)
dev.off()


