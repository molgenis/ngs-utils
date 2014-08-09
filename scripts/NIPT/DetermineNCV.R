library(gridExtra)

removeCols <- function(transposedMatrix, chromo)
{
  transposedMatrix[, which(colnames(transposedMatrix)==chromo)] <- NULL
  
  return(transposedMatrix)
}
GetSampleFracs <- function(sample)
{
  autoChromosomes <-c(1:12, 14:17, 19:20, 22) 
  sampleF <- read.delim(paste(sample, ".forward.corrected.bins.table.tsv", sep = ""), header = TRUE, sep = "\t", quote = "\"",
                        dec = ".", fill = TRUE)
  
  sampleR <- read.delim(paste(sample, ".reverse.corrected.bins.table.tsv", sep = ""), header = TRUE, sep = "\t", quote = "\"",
                        dec = ".", fill = TRUE)
  totalF <- sum(sampleF[autoChromosomes,])
  totalR <- sum(sampleR[autoChromosomes,])
  totalReadsChrFocus <<- sum(sampleF[chromo.focus,]) + sum(sampleR[chromo.focus, ])
  sampleFracs <- matrix(0, nrow=44, ncol=1)
  for (i in 1:22)
  {
    fraction <- sum(sampleF[i,]) / totalF
    sampleFracs[i,1] <- fraction
    
  }
  for (i in 1:22)
  {
    fraction <- sum(sampleR[i,]) / totalR
    sampleFracs[i + 22,1] <- fraction
    
  }
  rows<- c(paste("Chr", 1:22, "F", sep =""), paste("Chr", 1:22, "R", sep =""))
  rownames(sampleFracs) <- rows
  sampleFracs <- as.data.frame(t(sampleFracs))
  return(sampleFracs)
}

predictSamples <- function(data)
{
 predict(model1, data)
}

GetModel <- function(setNumber)
{
  pred1 <- transposedSamplesControl[, which(colnames(transposedSamplesControl)==predictors[setNumber,1])]
  pred2 <- transposedSamplesControl[, which(colnames(transposedSamplesControl)==predictors[setNumber,2])]
  pred3 <- transposedSamplesControl[, which(colnames(transposedSamplesControl)==predictors[setNumber,3])]
  pred4 <- transposedSamplesControl[, which(colnames(transposedSamplesControl)==predictors[setNumber,4])]
  return (model1 <- lm(chrFocusReads ~ pred1 + pred2 + pred3 + pred4))
}

GetPredictDataFrame <- function(setNumber, sourceFrame)
{
  return(data.frame(pred1 = sourceFrame[, which(colnames(sourceFrame)==predictors[setNumber,1])],
                    pred2 = sourceFrame[, which(colnames(sourceFrame)==predictors[setNumber,2])],
                    pred3 = sourceFrame[, which(colnames(sourceFrame)==predictors[setNumber,3])],
                    pred4 = sourceFrame[, which(colnames(sourceFrame)==predictors[setNumber,4])])) 
}

GetPredictorSets <- function(file)
{
  return(read.csv(file, header = TRUE, sep = "\t",
           dec = ".", fill = TRUE, comment.char = ""))
}

GetMean <- function(matrixRow, sampleRow)
{
  return(mean(as.matrix(matrixRow)))
}
GetSD <- function(matrixRow)
{
  return(sd(as.matrix(matrixRow)))
}
MakeTable <- function(colls)
{
  table <- matrix(0, nrow=4, ncol=length(colls)) 

  rowws <- c("Set 1", "Set 2", "Set 3", "Set 4")
  colnames(table) <- colls
  rownames(table) <- rowws
  return (table)
}

GetRegressionTable <- function(predictor)
{
  extractMatrix <- as.matrix(predictor)
  regressionTable <- matrix("", nrow=4, ncol=2) 

  cols <- c("Predictor Sets", "Regression Factor")
  rows <- c("Set 1", "Set 2", "Set 3", "Set 4")

  colnames(regressionTable) <- cols
  rownames(regressionTable) <- rows
  for (i in 1:4)
    {
    ps <- as.character(as.vector(extractMatrix[i,c(1:4)] ))
    regressionFactors <- as.vector(extractMatrix[i,12:15])
    regressionFactors <- paste(regressionFactors, ps, sep ="")
    regressionFactors <- paste(regressionFactors, collapse=" + ")
    ps <- paste(ps, collapse=" ")
    intercept <- as.vector(extractMatrix[i,11])
    intercept <- paste(intercept, " + ")
    regressionFactors <- paste(intercept, regressionFactors, sep="")
    regressionTable[i,] <- c(ps, regressionFactors)
  }
  return(regressionTable)
}

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

DetermineNCV <- function(totalReadsChrFocus)
{
  cv <<- (1.15 * (1 / sqrt(totalReadsChrFocus)))
  if(sd(normControls) < (1.15 * (1 / sqrt(totalReadsChrFocus))))
  {
    Zvalues <- (normControls - mean(normControls)) / (1.15*(1 / sqrt(totalReadsChrFocus)))
    Zvalues[51] <- (normSample - mean(normControls)) / (1.15*(1 / sqrt(totalReadsChrFocus)))
    corcof <<- "Yes"
  }
  else
  {
   Zvalues <- (normControls - mean(normControls)) / sd(normControls)
   Zvalues[51] <- (normSample - mean(normControls)) / sd(normControls)
   corcof <<- "No"
  }
  return (Zvalues)
}

WritePSSetsPDFs <- function(Zvalues, i)
{
  setwd(args[5])
  pdf(paste(args[3],"_Chromosome_", chromo.focus, "_Predictorset_4Predictors_Barplot_", i, ".pdf", sep = ""))
  barplot(Zvalues, col = ifelse(Zvalues < 6,'green','red'), axisnames = FALSE, 
          main = paste("Sample NCV = ", round(Zvalues[51], digits = 2), sep=""), xlab = corcof)
  dev.off()
  pdf(paste(args[3],"_Chromosome_", chromo.focus, "_Predictorset_4Predictors_Boxplot_", i, ".pdf", sep = ""))
  boxplot(Zvalues[1:50], col = ifelse(Zvalues < 6,'green','red'), main = paste("Mean control group = ", controlAv))
  dev.off()
  d <- density(Zvalues[1:50])
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
colsDiagnosticTable <- c("VC (%)", "Normal Distributed", "Zscore")
#Reads the chromosomal fraction table
setwd(args[5])
allmatrix <- read.delim(args[1], header = TRUE, sep = "\t", quote = "\"",
                        dec = ".", fill = TRUE)

transposedSamples <- as.data.frame(t(allmatrix))
chrFocusF <- paste("Chr", chromo.focus, "F", sep="")
chrFocusR <- paste("Chr", chromo.focus, "R", sep="")
chrFocusReads <- transposedSamples[[chrFocusF]] + transposedSamples[[chrFocusR]]
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

totalReadsChrFocus <- vector(mode = "numeric")

sampleFracs <- GetSampleFracs(args[3])

NCVTable <- MakeTable(colsNcvTable)
DiagnosticsTable <- MakeTable(colsDiagnosticTable)
for (i in 1:4)
{

cv <- 1  
  
moden <- GetModel(i)

corcof <- "Undetermined"

predictSample <- GetPredictDataFrame(i, sampleFracs)

predictedSample <-  predict(moden, predictSample)

observedSample <- sampleFracs[[chrFocusF]] + sampleFracs[[chrFocusR]]

normSample <- (sampleFracs[[chrFocusF]] + sampleFracs[[chrFocusR]]) / predictedSample

predictControl <- GetPredictDataFrame(i, transposedSamplesControl)
rownames(predictControl) <- rownames(transposedSamplesControl)

predictedControls <- predict(moden, predictControl)

normControls <- (chrFocusReads) / predictedControls

print(mean(normControls))

Zvalues <- DetermineNCV(totalReadsChrFocus)

controlAv <- round(mean(Zvalues[1:50]), digits = 2)
predictors$ncv[i] <- round(Zvalues[51], 3)

NCVTable[i,] <- c(round(predictedSample,5), round(observedSample,5), round(predictedSample / observedSample,5), round(mean(normControls),5), round(sd(normControls),6), 
                  round(sd(normControls) / mean(normControls),5), round(totalReadsChrFocus,0),  round(cv,4), corcof, round(shapiro.test(Zvalues[1:50])$p.value,3), round(Zvalues[51], 3))
normalDistributed <- "Yes"
if(shapiro.test(Zvalues[1:50])$p.value < 0.05)
{
  normalDistributed <- "No"
}

cvUsed <- round(sd(normControls) / mean(normControls),5) * 100

if(corcof == "Yes")
{
  cvUsed <- round(cv,4) * 100
}

DiagnosticsTable[i,] <- c(cvUsed, normalDistributed, round(Zvalues[51], 3)) 

WritePSSetsPDFs(Zvalues, i)
}

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

nv <- PrepareRegressionTableForPrint(NCVTable, chromo.focus)
pdf(paste(args[3], "_Chromosome", chromo.focus, "_NCVTable.pdf", sep =""), height=5, width=14)
grid.draw(nv)
dev.off()

gt <- PrepareRegressionTableForPrint(regressionTable, chromo.focus)
pdf(paste(args[3], "_Chromosome", chromo.focus, "_RegressionTable.pdf", sep=""), height=5, width=10)
grid.draw(gt)
dev.off()


