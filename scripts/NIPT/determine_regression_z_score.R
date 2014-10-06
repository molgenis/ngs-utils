####################################################################################################
# This script calculates the Z score between the predicted chromsomal fractions and 
# the observed chromosomal fraction

#Libs for nice pdf tables
library(gridExtra)
#Function that removes a column from a matrix
#################################################################################################################

#Script
#Stores the command line arguments in a vector
#args[1] = table with chromosomal fraction of all 50 control samples
#args[2] = table with predictor sets and relevant statistics
#args[3] = sample ID
#args[4] = output, table with predictor sets and relevant statistics with added Z score column
#args[5] = temp directory where files produced during run of pipeline are stored
#args[6] = chromosome of focus, either 13, 18 or 21
args<-commandArgs(TRUE)

#Gets the chromosome (13, 18 or 21) as a numeric
chromo.focus = as.numeric(args[6])
colsNcvTable <- c("Sample PS", "Observed Sample", "Observed / Predicted", "Mean Control", "SD Control", "Practical VC",
                  "Reads", "Theoretical VC", "Theoretical Used", "P-Value Shapiro", "Zscore")
colsDiagnosticTable <- c("VC", "Normal Distributed", "Zscore")
#Reads the chromosomal fraction table
setwd(args[5])
allmatrix <- read.delim(args[1], header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE)

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
  normControls <- chrFocusReads / predictedControls


  #Determines Z scores
  Zvalues <- DetermineNCV(totalReadsChrFocus)
  #mean Z score control group
  controlAv <- round(mean(Zvalues[1:length(Zvalues) -1]), digits = 2)
  #Add Z score 
  predictors$ncv[i] <- round(Zvalues[length(Zvalues)], 3)
  #Add row to NCVtable
  NCVTable[i,] <- c(round(predictedSample,5), round(observedSample,5), round(observedSample / predictedSample ,5), round(mean(normControls),5), round(sd(normControls),6), 
                  round(sd(normControls) / mean(normControls),5), round(totalReadsChrFocus,0),  round(cv,4), corcof, round(Shapiro.test(Zvalues[1:length(Zvalues) -1])$p.value,3), round(Zvalues[length(Zvalues)], 3))

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


