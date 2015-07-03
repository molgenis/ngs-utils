SumChiScores = function(bins.list) {
  # Scale each sample in bins.list so that all samples have the same number of reads
  bins.overall.mean = sum(as.numeric(unlist(bins.list))) / length(bins.list)
  mean.correction = function(m) m * bins.overall.mean / sum(m)
  bins.list.scaled = lapply(bins.list, mean.correction)
  
  # Calculate the expected value (= mean) per bin
  bins.sum.corrected = Reduce("+", bins.list.scaled)
  bins.scaled.expected = bins.sum.corrected / length(bins.list)
  
  # Calculate chi-score per bin per sample
  bins.chi.score = bins.list.scaled
  for (i in 1:length(bins.list)) {
    bins.chi.score[[i]] = (bins.scaled.expected - bins.chi.score[[i]])^2 / bins.scaled.expected
  }
  # Return sum of 
  return(Reduce("+", bins.chi.score))
}

# Corrects overdispersed bins in bins.list
CorrectBins = function(bins.list, chi.sum.bins, control.name, n.best.control.samples) {
  # About the input parameters:
  # bin.list = list(control sample 1, control sample 2, ..., sample of interest)
  # chi.sum.bins is determined only on control samples
  degrees.of.freedom = n.best.control.samples- 1 # number of control samples minus one
 
  # Convert chi squares to a normal distribution score 
  chi.sum.bins.normalized = (chi.sum.bins - degrees.of.freedom) / (sqrt( 2 * degrees.of.freedom))

  # Derive correction factor for the read counts
  chi.sum.bins.correction.factor = as.matrix(chi.sum.bins / degrees.of.freedom)
 
  # For each sample, correct extreme bins with high chi-square value
  index = which(chi.square.cut.off < chi.sum.bins.normalized) # Variation between these bins is considered too large
  for (i in 1:length(bins.list)) {
    m = as.matrix(bins.list[[i]])
    m[index] = m[index] / chi.sum.bins.correction.factor[index] # Therefore, we correct them here
    bins.list[[i]] = m
  }
  # Return corrected bin.list
  return(bins.list) 
}

ChiCorrect <- function(sample.forward, sample.reverse, control.bins.forward, control.bins.reverse)
{
    
  n.best.control.samples  = length(control.bins.forward)
  control.bins <- NULL
  for (p in 1:length(control.bins.forward))
  {
    control.bins[[p]] <- control.bins.forward[[p]] + control.bins.reverse[[p]]
  }
  # Calculate the chi square score per bin, based on only control samples
  chi.sum.bins = SumChiScores(bins.list = control.bins)
  
  control.bins.forward[[length(control.bins.forward) + 1]] <- sample.forward
  control.bins.reverse[[length(control.bins.reverse) + 1]] <- sample.reverse
  
  # Applies correction to bins
  correct.bins.list.forward = CorrectBins(bins.list = control.bins.forward, chi.sum.bins = chi.sum.bins, n.best.control.samples = n.best.control.samples)
  correct.bins.list.reverse = CorrectBins(bins.list = control.bins.reverse, chi.sum.bins = chi.sum.bins, n.best.control.samples = n.best.control.samples)
  
  sample.bins.corrected.forward <- correct.bins.list.forward[[length(correct.bins.list.forward)]]
  sample.bins.corrected.reverse <- correct.bins.list.reverse[[length(correct.bins.list.reverse)]]
  
  correct.bins.list.forward[[length(correct.bins.list.forward)]] <- NULL
  correct.bins.list.reverse[[length(correct.bins.list.reverse)]] <- NULL
  
  return (list(correct.bins.list.forward, sample.bins.corrected.forward, correct.bins.list.reverse, sample.bins.corrected.reverse))
}