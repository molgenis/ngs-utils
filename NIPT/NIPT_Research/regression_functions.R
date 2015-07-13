#This function selects the next predictor. 
GetNextPredictor <- function(samples, frac.reads.chr.trisomy.observed, predictors, chromosomal.frac.control)
{
  #Vector to hold r square for predictions
  Fstatistics <- NULL
  chr.candidates <- colnames(samples)
  #if step is 1, meaning first predictor is selected 
  for (i in 1:length(colnames(samples)))
  {
    current.pred.set <- c(predictors, chr.candidates[i])
    chr.frac.df <- as.matrix(cbind(chromosomal.frac.control[ , current.pred.set]))
    model <- lm(frac.reads.chr.trisomy.observed ~ chr.frac.df)
    Fstatistics[i] <- summary(model)$fstatistic[1]
  }
  #Orders the adjusted r squared values by index in decreasing order
  FstatisticsOrdered <- order(Fstatistics, decreasing = TRUE)
  #Gets all chromosomes remaining in the candidate predictors (for instance, Chr1F, Chr2F etc)
  chromosomes <- colnames(samples)
  #Builds the 'best' model again
  best.pred.set <- c(predictors, chromosomes[FstatisticsOrdered[1]])
  best.chr.frac.df <- as.matrix(cbind(chromosomal.frac.control[ , best.pred.set]))
  best.model <- lm(frac.reads.chr.trisomy.observed ~ best.chr.frac.df)
  #Stores predictor and F statistic in a list
  predictor.fstat <- list(chromosomes[FstatisticsOrdered[1]], summary(best.model)$fstatistic)
  return(predictor.fstat)
}

GetPredictors <- function(chr.int, chromosomal.frac.control, frac.reads.chr.trisomy.observed)
{
  # Get #reads on chromosome with potential trisomy for each of the control samples
  frac.reads.chr.trisomy.observed = as.matrix(chromosomal.frac.control[,chr.int]) + as.matrix(chromosomal.frac.control[,(chr.int + 22)])
 #Define potentially trisomic chromosomes
  chr.potential.trisomic <- c(paste("Chr" ,chromosomes.trisomy, "F", sep=""),paste("Chr" ,chromosomes.trisomy, "R", sep=""))
  #Get predictor sets
  predictor.list.force.variable <- SelectModelsForceApproach(elements = max.predictors, min.predictors=1, 
                                                             frac.reads.chr.trisomy.observed=frac.reads.chr.trisomy.observed)
  predictor.list.force.fixed <- SelectModelsForceApproach(elements=max.predictors, min.predictors = max.predictors,
                                                          frac.reads.chr.trisomy.observed=frac.reads.chr.trisomy.observed)
  predictor.list.regression.variable <- SelectModelsForwardRegressionApproach(n.of.predictors = "variable", chr.potential.trisomic=chr.potential.trisomic, 
                                                                              frac.reads.chr.trisomy.observed=frac.reads.chr.trisomy.observed)
  predictor.list.regression.fixed <- SelectModelsForwardRegressionApproach(n.of.predictors = "fixed", chr.potential.trisomic=chr.potential.trisomic,
                                                                           frac.reads.chr.trisomy.observed=frac.reads.chr.trisomy.observed)
  #Converts colnumbers to chromosomes. 
  predictor.list.force.variable <- lapply(predictor.list.force.variable, ConvertColNumbersToChr, col.names= col.names)
  predictor.list.force.fixed <- lapply(predictor.list.force.fixed, ConvertColNumbersToChr, col.names= col.names)
  return (list(predictor.list.force.fixed, predictor.list.force.variable, predictor.list.regression.fixed, 
               predictor.list.regression.variable))
}

SelectModelsForceApproach <- function(elements, min.predictors, frac.reads.chr.trisomy.observed)
{
  predictor.list.force.try <- list()
  for (model in 1:n.models)
  {
    #Predictors that are already in use 
    used.predictors <- c(unlist(predictor.list.force.try), chromosomes.trisomy)
    #Get list with a possible subsets of candidate chromosomses
    f.stats <- NULL
    max.f.stat.model <- list()
    max.f.stat.value <- NULL
    #Build every model possible and retrieve F statistic
    for(elements in min.predictors:max.predictors)
    {
      possible.predictors <- as.list(set_combn(chromosomes.background[-used.predictors], m=elements))
      
      f.stats <- unlist(lapply(possible.predictors, FUN = BuildModel, frac.reads.chr.trisomy.observed=frac.reads.chr.trisomy.observed, 
                               chromosomal.frac.control=chromosomal.frac.control))
      
      max.f.stat.model[[elements]] <- as.numeric(possible.predictors[[which.max(unlist(f.stats))]])
      max.f.stat.value[elements]<-  max(unlist(f.stats))
    }
    #Select subset which has highest F statistic
    predictor.list.force.try[[model]] <- unlist(max.f.stat.model[[which.max(max.f.stat.value)]])
  }
  return(predictor.list.force.try)
}
#Selects models using forward regression approach
SelectModelsForwardRegressionApproach <- function(n.of.predictors, chr.potential.trisomic, frac.reads.chr.trisomy.observed)
{
  predictor.list <- list()
  for (model in 1:n.models)
  {
    #Define (or empty) variables 
    predictors <- NULL
    predictor.fstat <-NULL
    predictors.complementary <- list()
    fstat <- NULL
    #Select predictors in iterative way until max number of predictors is reached
    for (predictor in 1:max.predictors)
    {
      #Data frame with potential predictors, i.e. all chromosomes minus potential trisomic chromosomes, previously used chromosomes and their complementary strand
      potential.predictors <- chromosomal.frac.control[,colnames(chromosomal.frac.control)[!colnames(chromosomal.frac.control) %in% c(chr.potential.trisomic,  unlist(predictors.complementary), unlist(predictor.list))]]
      
      predictor.fstat <- GetNextPredictor(samples = potential.predictors, frac.reads.chr.trisomy.observed=frac.reads.chr.trisomy.observed, predictors = predictors, 
                                          chromosomal.frac.control=chromosomal.frac.control)
      predictors[predictor] <- predictor.fstat[[1]]
      predictors.complementary[[predictor]] <- unique(c(predictor.fstat[[1]], sub("F", "R", predictor.fstat[[1]]), sub("R", "F", predictor.fstat[[1]])))
      fstat[predictor] <- predictor.fstat[[2]][1]
    }
    #If the number of predictors is variable, select the set with the highest F statistic
    if (n.of.predictors == "variable")
    {
      predictor.list[[model]] <- predictors[1:which.max(fstat)]
    }
    #If the number of predictors use current ser
    if (n.of.predictors == "fixed")
    {
      predictor.list[[model]] <- predictors
    }
  }
  return(predictor.list)
}

ConvertColNumbersToChr <- function(col.numbers, col.names)
{
  return(col.names[col.numbers])
}

PredictTrisomy <- function(predictors, chr.int,  control.samples, sample.bins.f, 
                           sample.bins.r, sample.of.interest)
{
  #Determines chromosomal fraction of chromosome of interest controls
  frac.reads.chr.trisomy.observed = as.matrix(chromosomal.frac.control[,chr.int]) + as.matrix(chromosomal.frac.control[,(chr.int + 22)])
  #Practial VC
  vc.prac <- 1.15 * (1 / sqrt(sum(sample.bins.f[chr.int,]) + sum(sample.bins.r[chr.int,])))
  #Build model 
  mod <- BuildFullModel(control.samples[predictors], chr.of.interest.fractions =  frac.reads.chr.trisomy.observed)
  #Determine ratio observed/predicted fractions
  ratio <-  frac.reads.chr.trisomy.observed /  mod$fitted.values
  #Theoretical vc
  vc.theo <- sd(ratio) 
  #Select which VC to use
  vc <- max(c(vc.theo, vc.prac))
  #Predict sample
  sample.predicted <- predict(mod, sample.of.interest[predictors])
  #Ratio sample
  sample.ratio <- sum(sample.of.interest[c(chr.int,(chr.int+22))]) / sample.predicted  
  #Determine Z scores
  z.score.sample <- (sample.ratio - mean(ratio)) / vc
  z.score.controls <- (as.numeric(ratio) - mean(ratio)) / vc
  #Shapiro-Wilk test for normality
  shapiro.regression <- shapiro.test(z.score.controls)
  #Appens line to output file
  AppendLine(line=c(paste(predictors, collapse =" "), z.score.sample, vc,shapiro.regression$p.value), file = output.file)
  
}
BuildFullModel <- function(control.samples, chr.of.interest.fractions)
{
  return (lm(chr.of.interest.fractions ~ ., data=control.samples))
}

BuildModel <- function(predictors, frac.reads.chr.trisomy.observed,  chromosomal.frac.control)
{
  predictors = unlist(predictors)
  model<- lm(frac.reads.chr.trisomy.observed ~ chromosomal.frac.control[, predictors])
  return(summary(model)$fstatistic[1])
}