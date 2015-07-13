#Gets subset of size m from candidates
RetrieveSubsets <- function(candidates, number.of.elements)
{
  return(set_combn(candidates, m=number.of.elements))
}

CalculateVariation <- function(chromosomal.frac.control.reads, denominators, chr.of.interest)
{
  possible.denominators <- unlist(denominators)
  # If length of possible denominators= 1, so function rowSums does not work 
  if (length(possible.denominators) == 1)
  {
    mean.subset <- mean(chromosomal.frac.control.reads[,chr.of.interest] / chromosomal.frac.control.reads[,possible.denominators])
    sd.subset <- sd((chromosomal.frac.control.reads[,chr.of.interest] / (chromosomal.frac.control.reads[,possible.denominators])))
  }
  # Else, if function rowSums does work:
  else
  {
    mean.subset <- mean(chromosomal.frac.control.reads[,chr.of.interest] / rowSums(chromosomal.frac.control.reads[,possible.denominators]))
    sd.subset <- sd((chromosomal.frac.control.reads[,chr.of.interest] / (rowSums(chromosomal.frac.control.reads[,possible.denominators]))))
  }
  #Get vc
  coef.of.variation <- VC(sd.subset, mean.subset)
  return(coef.of.variation)
}

VC <- function(sd, mean)
{
  return((sd/mean)*100)
}

GetDenominators <- function(chr.of.interest, chromosomal.frac.control.reads)
{
  #Declare variables
  min.variation.subset <- list()
  min.variation.vc <- NULL
  
  for (n.of.elements in 1:max.elements)
  {
    #Possible denominators
    denominators <- as.list(RetrieveSubsets(candidates = control.chromosomes.single, number.of.elements = n.of.elements))
    #Determine variation for each subset
    variation <- unlist(lapply(FUN = CalculateVariation, X=denominators, 
                               chromosomal.frac.control.reads= chromosomal.frac.control.reads, chr.of.interest=chr.of.interest))
    #Select subset with least variation
    min.variation.subset[[n.of.elements]] <- as.numeric(denominators[[which.min(variation)]])
    #Also store variation of this subset
    min.variation.vc[n.of.elements] <- variation[which.min(variation)]
  }
  #Select subset with least variation
  denominators <- min.variation.subset[[which.min(min.variation.vc)]]
  return(denominators)
}
PredictTrisomyNCV <- function(chr.of.interest, denominators, chromosomal.frac.control.reads, sample.bins)
{
  # If length of possible denominators= 1, so function rowSums does not work 
  if (length(denominators) == 1)
  {
   tr.all <- chromosomal.frac.control.reads[,chr.of.interest] / chromosomal.frac.control.reads[,denominators]
  }
  # Else, if rowSums does work
  else
  {
    tr.all <-chromosomal.frac.control.reads[,chr.of.interest] / rowSums(chromosomal.frac.control.reads[,denominators])
  }
  tr.sample <- sum(sample.bins[chr.of.interest,]) / sum(sample.bins[denominators,])
  tr.mean <- mean(tr.all)
  tr.sd <- sd(tr.all)
  #Z score controls and sample
  control.Zscore <- (tr.all - tr.mean) / tr.sd
  sample.Zscore <- (tr.sample - tr.mean) / tr.sd
  #VC
  vc <- tr.sd / tr.mean
  #Shapiro-Wilk test for normality
  shapiro.NCV <- shapiro.test(control.Zscore)
  #Append line to output file
  AppendLine(line=c(paste(paste("Chr",denominators,sep=""), collapse =" "), sample.Zscore, vc,shapiro.NCV$p.value), file = output.file)
}