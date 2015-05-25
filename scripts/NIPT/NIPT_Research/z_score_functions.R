#Gets Z scores control group
ZscoresControl <- function(chromosomal.frac.control, chr.int)
{
  return(scale(as.matrix(chromosomal.frac.control[,chr.int]) + as.matrix(chromosomal.frac.control[,(chr.int + 22)])))
}

ZscoreSample <- function(chr.int, chromosomal.frac.control,  sample.of.interest)
{
  #Determines chromosomal fraction of chromosome of interest for sample and controls
  sample.frac <- sample.of.interest[,chr.int] + sample.of.interest[,(chr.int+22)]
  control.fracs <- chromosomal.frac.control[,chr.int] + chromosomal.frac.control[,(chr.int+22)]
  #Z score sample
  z.score.sample <- (sample.frac - mean(control.fracs)) / sd(control.fracs)
  #Gets VC
  vc <-  sd(control.fracs) / mean(control.fracs)
  #Shapiro-Wilk test for normality
  shapiro.z <- shapiro.test(ZscoresControl(chromosomal.frac.control = chromosomal.frac.control, chr.int = chr.int))
  #Appens line to output file
  AppendLine(line=c(" ", z.score.sample, vc, shapiro.z$statistic), file = output.file)
}