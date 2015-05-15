#Determines the mean reads per bin for a given chromosome, only bins with reads are considered
MeanReadsChrom <- function(chromosome)
{
  return (sum(chromosome) / length(which(chromosome > 0)))
}
#Determines the square of the SD for the bincounts for a given chromosome, divided by bins that have reads
SdReadsChrom <- function(chromosome)
{
  return (sd(chromosome[which(chromosome > 0)])^2 / length(which(chromosome > 0)))
}

BonferonniCorrected <- function (sample.bins)
{
  sample.bins[which(is.na(sample.bins))] = 0
  #Gets means bin counts for every chromosome
  mean.reads.chromosome.bin <- apply(sample.bins, 1 , MeanReadsChrom)
  #Gets squared SD bincounts for a given chromosome divied by number of bins
  sd.reads.chromosome <- apply(sample.bins,1, SdReadsChrom)
  #comparison counter
  counter <- 0
  Z.scores <- NULL
  compared <- NULL
  #Outer loop to compare chromosomes
  for (first.chromosome in chromosomes.background.single)
  {
    #Sets numbers for inner loop
    comparison <- 1:22
    #Chromosome does not need to compared to itself
    comparison <- comparison[-first.chromosome]
    #Inner loop
    for (second.chromosome in comparison)
    {
      #Determines Z score comparison
      Z.scores[counter] <- (mean.reads.chromosome.bin[first.chromosome] - mean.reads.chromosome.bin[second.chromosome]) /
                                 (sqrt(sd.reads.chromosome[first.chromosome] + sd.reads.chromosome[second.chromosome]))
      #Stores which chromosomes are compared
      compared[counter] <- c(paste(first.chromosome, second.chromosome, collapse = " v " ))
      #Counts goes up 
      counter <- counter + 1
    }
  }
  #Converts Z-score to P value
  raw.p.values<- 2*pnorm(-abs(Z.scores))
  #Applies Boferronni correction
  corrected.p.values <- p.adjust(raw.p.values, method = "bonferroni")
  #Splits the corrected P values to chromosome
  p.values.chromosomes <- split(corrected.p.values, ceiling(seq_along(corrected.p.values)/(length(chromosomes.background.single)-1)))
  #For every P value
  for(chromosome in 1:length(p.values.chromosomes))
  {
    #Retrieve compared chromosomes
    current.chromosome <- p.values.chromosomes[[chromosome]]
    #IF length is more then one, i.e. if the chromosome has a significant call append result line to output file
    if(length(current.chromosome) != 0 )
    {
      AppendLine(line = c(chromosome, length(current.chromosome[current.chromosome < 0.001]), (length(chromosomes.background.single)-1),
                          mean(current.chromosome)), file = output.file)
    }
  }
  #Stores all significant comparisons in with p value to bonferronni corrected file
  #AppendLine(line = c("Compared chromosomes", "Corrected P Value"), file = output.file.bon)
   #for (line in which(corrected.p.values < 0.001))
  #{
   # AppendLine(line = c(compared[line], corrected.p.values[line]), file = output.file.bon)
   
 #}

}