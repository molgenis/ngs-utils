GetFiles <- function(controlDir, strand)
{
  setwd(controlDir)
  files <- list.files(  pattern = paste("*", strand,"*", sep = ""))
  return (files)
}

# Return list with bins of best control files
GetControlFiles = function(control.file.base.name)
{
  control.file.bin = list()
  for (i in 1:length(control.file.base.name) ){
    binfile = as.matrix(read.delim(control.file.base.name[i], header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE))
    
    # Only select relevant chromosomes
    control.file.bin[[i]] = binfile[chromosomes.background.single, ]
  }
  
  return (control.file.bin)
}
ChromosomalFractionPerSample = function(control.file.names, bins.forward, bins.reverse)
{
  chr.frac = NULL
  for (i in 1:length(bins.forward))
  {
    # For sample 'i', get fwd and rev bins[chrom, bin.i] (1 ... bin.i ... number of bins)
    bins.fwd.i = bins.forward[[i]]
    bins.fwd.i <- bins.fwd.i[chromosomes.background.single,]
    bins.fwd.i[which(is.na(bins.fwd.i))] <- 0
    bins.rev.i = bins.reverse[[i]]
    bins.rev.i <- bins.rev.i[chromosomes.background.single,]
    bins.rev.i[which(is.na(bins.rev.i))] <- 0
    
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
    bins <- bins[1:22,]
    bins[which(is.na(bins))] <- 0
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

AppendLineBon <- function(line, file)
{
  FF <- as.matrix(t(line))
  write.table(FF, file = "/Users/dirkdeweerd/Graduation/Data/voorbeeldbonferonni.csv", sep = ",", 
              col.names = FALSE, append=TRUE, row.names=FALSE)
}