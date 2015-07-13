DetermineMAD <- function(chr.int, chromosomal.frac.control, sample.of.interest)
{
  #Determines chromosomal fraction of chromosome of interest for sample and controls
  control.fracs <- chromosomal.frac.control[,chr.int] + chromosomal.frac.control[,(chr.int+22)]
  sample.frac <- sample.of.interest[,chr.int] + sample.of.interest[,(chr.int+22)]
  
  #Determines MAD score sample
  mad.score <- (sample.frac - median(control.fracs)) / mad(control.fracs)
  #Determines VC
  vc <- mad(control.fracs) / mean(control.fracs)
  #Determines MAD scores control groep
  mad.score.control <- (control.fracs - median(control.fracs)) / mad(control.fracs)
  #Shapiro-Wilk test for normality
  shapiro.mad <- shapiro.test(mad.score.control)
  #Appens line to output file
  AppendLine(line=c(" ", mad.score, vc, shapiro.mad$p.value), file = output.file)
}