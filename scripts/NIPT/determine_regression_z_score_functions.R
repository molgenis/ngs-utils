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
  sampleFracs <- as.data.frame(t(sampleFracs)) # 1 row
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
# LAYOUT Prepare table to store the regression factors
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
  cv <<- (1.15*(1 / sqrt(totalReadsChrFocus)))
  #if practical vc is lower than theoretical vc use theoretical
  if(sd(normControls) < cv))
  {
    Zvalues <- (normControls - mean(normControls)) / cv
    Zvalues[length(Zvalues)+1] <- (normSample - mean(normControls)) / cv
    corcof <<- "Yes"
  }
  # if not
  else
  {
   Zvalues <- (normControls - mean(normControls)) / sd(normControls)
   Zvalues[length(Zvalues)+1] <- (normSample - mean(normControls)) / sd(normControls)
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
          main = paste("Sample NCV = ", round(Zvalues[length(Zvalues)], digits = 2), sep=""), xlab = corcof)
  dev.off()
  pdf(paste(args[3],"_Chromosome_", chromo.focus, "_Predictorset_4Predictors_Boxplot_", i, ".pdf", sep = ""))
  boxplot(Zvalues[1:length(Zvalues) -1], col = ifelse(Zvalues < 6,'green','red'), main = paste("Mean control group = ", controlAv))
  dev.off()
  d <- density(Zvalues[1:length(Zvalues) -1])
  pdf(paste(args[3], "_Chromosome_", chromo.focus,"_Predictorset_4Predictors_Density_", i, ".pdf", sep = ""))
  plot(d)
  dev.off()
}
