ChromosomalFractionPerSample = function(control.file.names, bins.forward, bins.reverse)
{
  chr.frac = NULL
  for (i in 1:length(bins.forward))
  {
    # For sample 'i', get fwd and rev bins[chrom, bin.i] (1 ... bin.i ... number of bins)
    bins.fwd.i = bins.forward[[i]]
    bins.rev.i = bins.reverse[[i]]
    
    # Determine fraction of reads on each chromosome-direction
    chr.frac.fwd.i = rowSums(bins.fwd.i) / sum(bins.fwd.i)
    chr.frac.rev.i = rowSums(bins.rev.i) / sum(bins.rev.i)
    
    # Add these fractions as column
    chr.frac = rbind(chr.frac, c(chr.frac.fwd.i, chr.frac.rev.i))
  }
  
  # Set row and column names
  rownamen <- sub("^([^.]*).*", "\\1", control.file.names) # sample names
  rownames(chr.frac) = rownamen # sample names
  colnames(chr.frac) = c(paste("Chr", chromosomes.background.single, "F", sep =""), paste("Chr", chromosomes.background.single, "R", sep ="")) # 1F, 2F, ... nF, 1R, 2R, ... nR
  
  return(chr.frac)
}

ChromosomalReadsPerSample = function(control.file.names, bins.forward, bins.reverse)
{
  chr.frac = NULL
  for (i in 1:length(bins.forward))
  {
    # For sample 'i', get fwd and rev bins[chrom, bin.i] (1 ... bin.i ... number of bins)
    bins = bins.forward[[i]] + bins.reverse[[i]]
    
    
    # Determine fraction of reads on each chromosome-direction
    chr.frac.row = rowSums(bins) 
    
    
    # Add these fractions as column
    chr.frac = rbind(chr.frac, chr.frac.row)
  }
  
  # Set row and column names
  rownamen <- sub("^([^.]*).*", "\\1", control.file.names) # sample names
  rownames(chr.frac) = rownamen # sample names
  colnames(chr.frac) = paste("Chr" , chromosomes.background.single, sep="")
  return(chr.frac)
}

AppendLine <- function(line, file)
{
  FF <- as.matrix(t(line))
  write.table(FF, file = file, sep = ",", 
              col.names = FALSE, append=TRUE, row.names=FALSE)
}

AppendLineBon <- function(line)
{
  FF <- as.matrix(t(line))
  write.table(FF, file = "/Users/dirkdeweerd/Graduation/Data/voorbeeldbonferonni.csv", sep = ",", 
              col.names = FALSE, append=TRUE, row.names=FALSE)
}