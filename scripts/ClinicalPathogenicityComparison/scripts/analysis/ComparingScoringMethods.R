#############################
# Comparing Scoring Methods #
#############################
# This script can be used to calculate performance differences between three pathogenicity scoring methods (DANN, CADD and FitCon).
# With the use of Youden's index the performance of scoring pathogenic observed data with the scoring added by DANN, CADD or FitCon.
# 
# The results of the analysis shows is a CSV with statistical outcomes for each scoring method of the Youden's index, 
# the standard error of this index, fishers exact p-value and the prositive predicted values of sample taken from the pathogenicly expected regions having at least 20 variants.
# Also the amount of True positives, false positives, true negatives and false negatives are given for each scoring method.
# 
# 


#############
# libraries #
#############
#install.packages("ggplot2")
#install.packages("stringr")
#install.packages("getopt")
#install.packages("VariantAnnotation")

library("ggplot2")
library("stringr")
library("VariantAnnotation") 
library("getopt")


#################
# Youdens index #
#################
# a <- correctlyDiagnose (True Positives)
# b <- falseNegativeDiagnose (False Negatives)
# c <- falsePositiveDiagnose (False Positives)
# d <- correctlyDiagnose (True Negatives)
#
# Uses the proposed method by Youden et al (http://onlinelibrary.wiley.com/doi/10.1002/1097-0142%281950%293:1%3C32::AID-CNCR2820030106%3E3.0.CO;2-3/epdf)
# to calculate the index of performance. This index will be 0 when the performance is equal
# and will reach unitiy if there is less simmilarity between the methods tested.
#  
youdensIndex <- function(a, b, c, d){
  return (((a * d) - (b * c)) / (( a + b ) * ( c + d )))
}


##########################
# Youdens standard error #
##########################
# a <- correctlyDiagnose (True Positives)
# b <- falseNegativeDiagnose (False Negatives)
# c <- falsePositiveDiagnose (False Positives)
# d <- correctlyDiagnose (True Negatives)
#
# This function can be used to calculate the standard error of the youdens index.
# By doing so it is possible to calculate these values between different scoring methods.
#
# source: "(http://onlinelibrary.wiley.com/doi/10.1002/1097-0142%281950%293:1%3C32::AID-CNCR2820030106%3E3.0.CO;2-3/epdf)"
#
youdensErr <- function(a, b, c, d){
  return (sqrt(((( a * b ) / ( a + b )^ 3 ) + (( c * d ) / ( c + d )^3 ))))
}


################################
# youdens index Standard Error #
################################
# Calculates the standard error based of two youden indexes from different scoring methods
#
# Source: "(http://onlinelibrary.wiley.com/doi/10.1002/1097-0142%281950%293:1%3C32::AID-CNCR2820030106%3E3.0.CO;2-3/epdf)"
#
youdensStdErr <- function(youdensErr1, youdensErr2){
  return (sqrt( ((youdensErr1^2) + (youdensErr2^2)) ))
}


############################
# T.test for Youdens index #
############################
# This t test checks if index1 and index2 differ in their standard error. 
# If so it indicates a loss or gain in performance between these two different methods.
#
# Source: "(http://onlinelibrary.wiley.com/doi/10.1002/1097-0142%281950%293:1%3C32::AID-CNCR2820030106%3E3.0.CO;2-3/epdf)"
#
youdens.index.t.test = function(index1, index2,  standardError1, standardError2){
  diff <- abs( index1 - index2)

  return (diff / (sqrt(((standardError1^2) + (standardError2^2)))))
}


#############
# multiplot #
#############
# Allows for multiplotting based upon the number of columns specified
# Expects a number of plots to draw.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) 
  {
    print(plots[[1]])
    
  } 
  else
  {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


################
# analyse.Gene #
################
# performs an analysis by vcf object containing info fields derived from the MOLGENIS annotation framework
# input and a gene name and a boolean indicating if graphs to be plotted or not. 
#
#
analyse.Gene <- function(gene, vcf, plots, plotout){
  gene.datatable <- as.data.frame(info(vcf)[1:nrow(info(vcf)), 1:11][info(vcf)$GENE_SYMBOL == gene,])
  
  colours <- NULL
#get non pathogenic and pathogenic indici
  index.nonpathogenic = which( .02 < gene.datatable["GoNL_AF"] )
  index.nonpathogenic = c(index.nonpathogenic, which(.02 < gene.datatable["Thousand_Genomes_AF"]))
  index.nonpathogenic = c(index.nonpathogenic, which(.02 < gene.datatable["EXAC_AF"]))
  index.pathogenic = which( "5" == gene.datatable["CLINVAR_CLNSIG"] )
# set all colour labels to unknown 
  colours = rep("Unknown", nrow(gene.datatable))
  for (i in 1:nrow(gene.datatable["CLINVAR_CLNSIG"]))
  {
    elt = gene.datatable$CLINVAR_CLNSIG[[i]]
    if (0 != length(elt))
    {
      if (is.element("5", strsplit(elt, "\\|")[[1]])){
        index.pathogenic = c(index.pathogenic, i)  
      }
    }
  } 
# if there are any indexes found indicating pathogenic variants then set the pathogenic colour labels
  if(length(index.pathogenic) > 0)
  {
    colours[index.pathogenic] <- "Pathogenic"
  }
# if there are any indexes found indicating Benign variants then set the Benign colour labels
  if(length(index.nonpathogenic) > 0)
  {
    colours[index.nonpathogenic] <- "Benign"
  }
  
# set colour labels to the gene table
  gene.datatable$colour <- colours

# sort each dataset based on the scoring method
  gene.datatable.sorted.fitcon <- gene.datatable[order(unlist(gene.datatable["FITCON_SCORE"])),]
  gene.datatable.sorted.cadd <- gene.datatable[order(unlist(gene.datatable["CADD"])),]
  gene.datatable.sorted.cadd.scaled <- gene.datatable[order(unlist(gene.datatable["CADD"])),]
  gene.datatable.sorted.dann <- gene.datatable[order(unlist(gene.datatable["DANN_SCORE"])),]

# get the cutpoints for each scoring method
  cadd.scaled.cutpoint.index <- treshold.calculate(unlist(gene.datatable.sorted.cadd.scaled["colour"]))
  cadd.scaled.plot.unknown.turnover.cutoff <- unknown.treshold.calculate(unlist(gene.datatable.sorted.cadd.scaled["colour"]))
  
  cadd.cutpoint.index <- treshold.calculate(unlist(gene.datatable.sorted.cadd["colour"]))
  cadd.plot.unknown.turnover.cutoff <- unknown.treshold.calculate(unlist(gene.datatable.sorted.cadd["colour"]))
  
  dann.cutpoint.index <- treshold.calculate(unlist(gene.datatable.sorted.dann["colour"]))
  dann.plot.unknown.turnover.cutoff <- unknown.treshold.calculate(unlist(gene.datatable.sorted.dann["colour"]))
  
  fitcon.cutpoint.index <- treshold.calculate(unlist(gene.datatable.sorted.fitcon["colour"]))
  fitcon.plot.unknown.turnover.cutoff <- unknown.treshold.calculate(unlist(gene.datatable.sorted.fitcon["colour"]))
 
# draw graphs if plots = true
  if(plots == TRUE){
    fitcon.plot<-qplot(unlist(gene.datatable.sorted.fitcon["FITCON_SCORE"]), xlab="Fitcon Score", ylab="Counts", main=paste("FitCon score distribution for gene:", gene), fill=factor(unlist(gene.datatable.sorted.fitcon$colour)))
    cadd.plot<-qplot(unlist(gene.datatable.sorted.cadd["CADD"]), xlab="CADD Score", ylab="Counts", main=paste("CADD score distribution gene:",gene), fill=factor(unlist(gene.datatable.sorted.cadd$colour)))
    cadd.scaled.plot<-qplot(unlist(gene.datatable.sorted.cadd.scaled["CADD_SCALED"]), xlab="CADD scaled score", ylab="Counts", main=paste("Scaled CADD score distribution for gene:", gene), fill=factor(unlist(gene.datatable.sorted.cadd.scaled$colour)))
    dann.plot<-qplot(unlist(gene.datatable.sorted.dann["DANN_SCORE"]), xlab="DANN Score", ylab="Counts", main=paste("DANN score distribution for gene:", gene), fill=factor(unlist(gene.datatable.sorted.dann$colour)))
    
    cadd.scaled.plot <- cadd.scaled.plot + geom_vline(xintercept = unlist(gene.datatable.sorted.cadd.scaled$CADD_SCALED)[cadd.scaled.cutpoint.index], colour="purple", linetype = "longdash")
    cadd.scaled.plot <- cadd.scaled.plot + geom_vline(xintercept = unlist(gene.datatable.sorted.cadd.scaled$CADD_SCALED)[cadd.scaled.plot.unknown.turnover.cutoff], colour="orange", linetype = "longdash")
    cadd.scaled.plot <- cadd.scaled.plot + guides(fill=guide_legend(title="Classification"))
    
    dann.plot <- dann.plot + geom_vline(xintercept = unlist(gene.datatable.sorted.dann$DANN_SCORE)[dann.cutpoint.index], colour="purple", linetype = "longdash")
    dann.plot <- dann.plot + geom_vline(xintercept = unlist(gene.datatable.sorted.dann$DANN_SCORE)[dann.plot.unknown.turnover.cutoff], colour="orange", linetype = "longdash")
    dann.plot <- dann.plot + guides(fill=guide_legend(title="Classification"))
    
    cadd.plot <- cadd.plot + geom_vline(xintercept = unlist(gene.datatable.sorted.cadd$CADD)[cadd.cutpoint.index], colour="purple", linetype = "longdash")
    cadd.plot <- cadd.plot + geom_vline(xintercept = unlist(gene.datatable.sorted.cadd$CADD)[cadd.plot.unknown.turnover.cutoff], colour="orange", linetype = "longdash")
    cadd.plot <- cadd.plot + guides(fill=guide_legend(title="Classification"))
    
    pdf(paste(plotout, gene, ".pdf", sep=""),width=15, height=7)
    
    if(!all(is.nan(unlist(gene.datatable.sorted.fitcon$FITCON_SCORE)) )){
      fitcon.plot <- fitcon.plot + geom_vline(xintercept = unlist(gene.datatable.sorted.fitcon$FITCON_SCORE)[fitcon.cutpoint.index], colour="purple", linetype = "longdash")
      fitcon.plot <- fitcon.plot + geom_vline(xintercept = unlist(gene.datatable.sorted.fitcon$FITCON_SCORE)[fitcon.plot.unknown.turnover.cutoff], colour="orange", linetype = "longdash")
      fitcon.plot <- fitcon.plot + guides(fill=guide_legend(title="Classification"))
      
      multiplot(dann.plot, cadd.scaled.plot, cadd.plot, fitcon.plot, cols = 1)
    }
    else
    {
      multiplot(dann.plot, cadd.scaled.plot, cadd.plot, cols = 1)
    }
    
    dev.off()
  }

  # get contingency table per method according to colour classification
  cadd.scaled.contigency.table <- calculate.contigency.table(unlist(gene.datatable.sorted.cadd.scaled["colour"]))
  cadd.contigency.table <- calculate.contigency.table(unlist(gene.datatable.sorted.cadd["colour"]))
  fitcon.contigency.table <- calculate.contigency.table(unlist(gene.datatable.sorted.fitcon["colour"]))
  dann.contigency.table <- calculate.contigency.table(unlist(gene.datatable.sorted.dann["colour"]))
  
  # calculate Youden's index per scoring method according to the contigency table
  cadd.scaled.youden.index <- youdensIndex(cadd.scaled.contigency.table[1,1], cadd.scaled.contigency.table[2,1], cadd.scaled.contigency.table[1,2], cadd.scaled.contigency.table[2,2])
  cadd.youden.index <-  youdensIndex(cadd.contigency.table[1,1], cadd.contigency.table[2,1], cadd.contigency.table[1,2],cadd.contigency.table[2,2])
  fitcon.youden.index <-  youdensIndex(fitcon.contigency.table[1,1], fitcon.contigency.table[2,1],fitcon.contigency.table[1,2], fitcon.contigency.table[2,2])
  dann.youden.index <-  youdensIndex(dann.contigency.table[1,1], dann.contigency.table[2,1], dann.contigency.table[1,2], dann.contigency.table[2,2])
  
  # calculate Youden's index per scoring method according to the contigency table according to colour classification
  cadd.scaled.youden.error <- pathogenicity.method.youden.error(cadd.scaled.contigency.table)
  cadd.youden.error <- pathogenicity.method.youden.error(cadd.contigency.table)
  fitcon.youden.error <- pathogenicity.method.youden.error(fitcon.contigency.table)
  dann.youden.error <- pathogenicity.method.youden.error(dann.contigency.table)
  
  # get t-value of difference between two methods based on youden's statistics 
  #function(index1, index2,  standardError1, standardError2){
  cadd.scaled.vs.cadd.youden.tvalue <- youdens.index.t.test(cadd.scaled.youden.index, cadd.youden.index, cadd.scaled.youden.error, cadd.youden.error)
  cadd.scaled.vs.fitcon.youden.tvalue <- youdens.index.t.test(cadd.scaled.youden.index, fitcon.youden.index, cadd.scaled.youden.error, fitcon.youden.error)
  cadd.scaled.vs.dann.youden.tvalue <- youdens.index.t.test(cadd.scaled.youden.index, dann.youden.index, cadd.scaled.youden.error, dann.youden.error)
  
  cadd.vs.cadd.scaled.youden.tvalue <- youdens.index.t.test(cadd.youden.index, cadd.scaled.youden.index, cadd.youden.error, cadd.scaled.youden.error)
  cadd.vs.fitcon.youden.tvalue <- youdens.index.t.test(cadd.youden.index, fitcon.youden.index, cadd.youden.error, fitcon.youden.error)
  cadd.vs.dann.youden.tvalue <- youdens.index.t.test(cadd.youden.index, dann.youden.index, cadd.youden.error, dann.youden.error)
  
  fitcon.youden.tvalue.vs.cadd.scaled.youden.tvalue <- youdens.index.t.test(fitcon.youden.index, cadd.scaled.youden.index, fitcon.youden.error, cadd.scaled.youden.error)
  fitcon.youden.tvalue.vs.cadd.youden.tvalue <- youdens.index.t.test(fitcon.youden.index, cadd.youden.error, fitcon.youden.error, cadd.youden.error)
  fitcon.youden.tvalue.vs.dann.youden.tvalue <- youdens.index.t.test(fitcon.youden.index, dann.youden.index, fitcon.youden.error, dann.youden.error)
  
  dann.youden.vs.cadd.scaled.youden.tvalue <- youdens.index.t.test(dann.youden.index, cadd.scaled.youden.index, dann.youden.error, cadd.scaled.youden.error)
  dann.youden.vs.cadd.youden.tvalue <- youdens.index.t.test(dann.youden.index, cadd.youden.index, dann.youden.error, cadd.youden.error)
  dann.youden.vs.fitcon.youden.tvalue <- youdens.index.t.test(dann.youden.index, fitcon.youden.index, dann.youden.error, fitcon.youden.error)

# gets the positive predicted values for each scoring method according to the previously determined cutoff index
  cadd.ppv <- calculate.ppv(unlist(gene.datatable.sorted.cadd$colour)[cadd.plot.unknown.turnover.cutoff:nrow(gene.datatable.sorted.cadd)])
  fitcon.ppv <- calculate.ppv(unlist(gene.datatable.sorted.fitcon$colour)[fitcon.plot.unknown.turnover.cutoff:nrow(gene.datatable.sorted.fitcon)])
  dann.ppv <- calculate.ppv(unlist(gene.datatable.sorted.dann$colour)[dann.plot.unknown.turnover.cutoff:nrow(gene.datatable.sorted.dann)])
  cadd.scaled.ppv <- calculate.ppv(unlist(gene.datatable.sorted.cadd.scaled$colour)[cadd.scaled.plot.unknown.turnover.cutoff:nrow(gene.datatable.sorted.cadd.scaled)])

# calculate fishers exact values according to their colour label 
  cadd.scaled.fisher <- fisher.test(cadd.scaled.contigency.table)
  cadd.fisher <- fisher.test(cadd.contigency.table)
  fitcon.fisher <- fisher.test(fitcon.contigency.table)
  dann.fisher <- fisher.test(dann.contigency.table)

# calculate the pathogenic and non pathogenic pvalue according to a division made at the middle of the values expecting to have the pathogenic at the high end and Benign at the low end
  cadd.scaled.fisher.pathogenic <- calculate.fisher(unlist(gene.datatable.sorted.cadd.scaled$colour)[1:round(length(unlist(gene.datatable.sorted.cadd.scaled["colour"]))/2)])
  cadd.scaled.fisher.nonpathogenic <- calculate.fisher(unlist(gene.datatable.sorted.cadd.scaled$colour)[round(length(unlist(gene.datatable.sorted.cadd.scaled["colour"]))/2): length(unlist(gene.datatable.sorted.cadd.scaled["colour"]))])
  cadd.fisher.pathogenic <- calculate.fisher(unlist(gene.datatable.sorted.cadd$colour)[1:round(length(unlist(gene.datatable.sorted.cadd["colour"]))/2)])
  cadd.fisher.nonpathogenic <- calculate.fisher(unlist(gene.datatable.sorted.cadd$colour)[round(length(unlist(gene.datatable.sorted.cadd["colour"]))/2): length(unlist(gene.datatable.sorted.cadd["colour"]))])
  dann.fisher.pathogenic <- calculate.fisher(unlist(gene.datatable.sorted.dann$colour)[1:round(length(unlist(gene.datatable.sorted.dann["colour"]))/2)])
  dann.fisher.nonpathogenic <- calculate.fisher(unlist(gene.datatable.sorted.dann$colour)[round(length(unlist(gene.datatable.sorted.dann["colour"]))/2): length(unlist(gene.datatable.sorted.dann["colour"]))])
  fitcon.fisher.pathogenic <- calculate.fisher(unlist(gene.datatable.sorted.fitcon$colour)[1:round(length(unlist(gene.datatable.sorted.fitcon["colour"]))/2)])
  fitcon.fisher.nonpathogenic <- calculate.fisher(unlist(gene.datatable.sorted.fitcon$colour)[round(length(unlist(gene.datatable.sorted.fitcon["colour"]))/2): length(unlist(gene.datatable.sorted.fitcon["colour"]))])

  gene.output <- data.frame(
    genesymbol=gene,
    cadd.scaled.truepositives=cadd.scaled.contigency.table[1,1],
    cadd.scaled.falsepositives=cadd.scaled.contigency.table[1,2],
    cadd.scaled.falsenegatives=cadd.scaled.contigency.table[2,1],
    cadd.scaled.truenegatives=cadd.scaled.contigency.table[2,2],
    cadd.truepositives=cadd.contigency.table[1,1],
    cadd.falsepositives=cadd.contigency.table[1,2],
    cadd.falsenegatives=cadd.contigency.table[2,1],
    cadd.truenegatives=cadd.contigency.table[2,2],
    fitcon.truepositives=fitcon.contigency.table[1,1],
    fitcon.falsepositives=fitcon.contigency.table[1,2],
    fitcon.falsenegatives=fitcon.contigency.table[2,1],
    fitcon.truenegatives=fitcon.contigency.table[2,2],
    dann.truepositives=dann.contigency.table[1,1],
    dann.falsepositives=dann.contigency.table[1,2],
    dann.falsenegatives=dann.contigency.table[2,1],
    dann.truenegatives=dann.contigency.table[2,2],
    cadd.ppv=cadd.ppv,
    fitcon.ppv=fitcon.ppv,
    dann.ppv=dann.ppv, 
    cadd.scaled.ppv=cadd.scaled.ppv,
    cadd.scaled.youden.error=cadd.scaled.youden.error,
    dann.youden.error=dann.youden.error,
    fitcon.youden.error=fitcon.youden.error,
    cadd.youden.error=cadd.youden.error,
    cadd.scaled.fisher=cadd.scaled.fisher$p.value,
    cadd.fisher=cadd.fisher$p.value,
    fitcon.fisher=fitcon.fisher$p.value,
    dann.fisher=dann.fisher$p.value,
    cadd.scaled.youden.index=cadd.scaled.youden.index,
    cadd.youden.index=cadd.youden.index,
    fitcon.youden.index=fitcon.youden.index,
    dann.youden.index=dann.youden.index
  )

  return(gene.output)
}


#################
# calculate.ppv #
#################
# Given a array containing pathogenic, Benign and unknown labels.
# Calculates the positive predicted values in the array of labels.
#
# Returns number of pathogenic labels / total counts of occurences in the array
calculate.ppv <- function(population){
  ppv <- 0 

  # ppv = number of pathogenic variants in the population divided by the total population
  ppv = (length(population[population == "Pathogenic"]) / length(population[population != "Uknown"]))

  return (ppv)
}


#################
# calculate.contigency.table #
#################
# Given a array containing pathogenic, Benign and pathogenic labels.
# Calculates a contingecy table according to the order of values.
# pathogenic variants are counted as positive values and benign values are rated as negative values.
#
# Returns contingency table according to order of Benign and pathogenic label order in the input.
calculate.contigency.table<- function(population){
  falsenegative = 0
  truenegative = 0
  falsepositive = 0
  truepositive = 0
  
  #detect the true negatives and the false negatives by counting the inconsistency of order of pathogenic and Benign
  for(i in 1:length(population)){
    if(i < length(population)){
      if(!("Unknown" %in% population[i])){
        
        # distinguish between false and true negative values 
        if("Pathogenic" %in% population[1:i] && "Benign" %in% population[i] ){
          
          falsenegative = falsenegative + 1
        }
        
        else if(!("Pathogenic" %in% population[1:i]) && "Benign" %in% population[i] ){
          truenegativeIndex=i
          truenegative= truenegative + 1
        }
      }
    }
  }
  
  #detect the true positives and the false positives by counting the inconsistency of the order of pathogenic and Benign
  for(i in length(population):1){
    if(i > 1){
      if(!("Unknown" %in% population[i] )){
        # distinguish between false and true positive values 
        if( "Benign" %in% population[i:length(population)] && "Pathogenic" %in% population[i]  ){
          falsepositive = falsepositive + 1
        }
        
        else if(!("Benign" %in% population[i:length(population)]) && "Pathogenic" %in% population[i] ){
          truepositive = truepositive + 1
          
        }
      }
    }
  }
  contingency.table<-t(matrix( c(truepositive, falsenegative, falsepositive, truenegative), 2, 2))
  return(contingency.table)
}


#####################################
# pathogenicity.method.youden.error #
#####################################
# Given a array containing pathogenic, Benign and unknown labels 
# performs the fishers exact test on the false positives, true positives, 
# false negatives and true negatives for a population.
# The values for the fishers test are determined by the order of the array of labels.
# Benign are expected to be found at the lower indici and pathogenic are expected
# to be at the end of the array. 
# Given these assumptions the false negative and positives are determined. 
# Returns a fishers.test object
# 
calculate.fisher <- function(population){
  falsenegative = 0
  truenegative = 0
  falsepositive = 0
  truepositive = 0

  #detect the true negatives and the false negatives by counting the inconsistency of order of pathogenic and Benign
  for(i in 1:length(population)){
    if(i < length(population)){
      if(!("Unknown" %in% population[i])){
        
        # distinguish between false and true negative values 
        if("Pathogenic" %in% population[1:i] && "Benign" %in% population[i] ){
          
          falsenegative = falsenegative + 1
        }
        
        else if(!("Pathogenic" %in% population[1:i]) && "Benign" %in% population[i] ){
          truenegativeIndex=i
          truenegative= truenegative + 1
        }
      }
    }
  }

  #detect the true positives and the false positives by counting the inconsistency of the order of pathogenic and Benign
  for(i in length(population):1){
    if(i > 1){
      if(!("Unknown" %in% population[i] )){
        # distinguish between false and true positive values 
        if( "Benign" %in% population[i:length(population)] && "Pathogenic" %in% population[i]  ){
          falsepositive = falsepositive + 1
        }
        
        else if(!("Benign" %in% population[i:length(population)]) && "Pathogenic" %in% population[i] ){
          truepositive = truepositive + 1
          
        }
      }
    }
  }
  
  contingency.table<-t(matrix( c(truepositive, falsenegative, falsepositive, truenegative), 2, 2))

  return(fisher.test(contingency.table))
}


#####################################
# pathogenicity.method.youden.error #
#####################################
# Given a array containing pathogenic and Benign labels 
# calculates the Youdens standard error for a population.
# The Youdens standard error is calculated by assessing the false negative, false positive, true positive and true negative 
# 
pathogenicity.method.youden.error <- function(contingency.table){
  truepositive = contingency.table[1,1]
  falsepositive = contingency.table[1,2]
  falsenegative = contingency.table[2,1]
  truenegative = contingency.table[2,2]

  youdens.err<-youdensErr(truepositive, falsenegative, falsepositive, truenegative)
  return (youdens.err)
}


##############################
# unknown.treshold.calculate #
##############################
# Given a array containing pathogenic and unknown labels 
# calculates the index of the threshold for the ratio of pathogenic and unknown labels in a population.
# the threshold is set at 0.2 unknown labels
# The index is of the first point reaching this threshold in a population of atleast 20 is returned
# 
unknown.treshold.calculate <- function(population){
  index = 0
  pathogenic.count = 0
  unknown.count = 0
  counts = 0
  #check when the first uncertainty occurs
  for(i in length(population):1){
    counts = counts + 1
    if( "Unknown" %in% population[i]){
      unknown.count <- unknown.count + 1
    
    }else if( "Pathogenic"  %in% population[i]){
      pathogenic.count <- pathogenic.count + 1
    }
    if(pathogenic.count != 0 || unknown.count != 0 ){
   
      if((unknown.count / (pathogenic.count+unknown.count) ) > 0.20){
        index <- i
        if (counts > 20){
          index <- index - 1
          if(pathogenic.count == 0 ){
            index = length(population)
          }
          break
        }
      }
    }
  }

  #exclude faulty value so plus one to the index
  return(index)
}

######################
# treshold.calculate #
######################
# Given a array containing Benign and pathogenic labels
# calculates the index of the threshold for the last 
# occurrence of a Benign label and returns the index
# 
treshold.calculate <- function(population){
  index = 0
  #check when the first uncertainty occurs
  for(i in length(population):1){
    if( "Benign" %in% population[i]){
      index = i 
      break
    }
  }
  #exclude faulty value so plus one to the index
  index = index+1
  return(index)
}


##########################
# annotated VCF analysis #
##########################
spec <- matrix(c(
  'input'     , 'i', 1, "character", "input VCF file annotated using the molgenis annotation software, following fields are required in the info collumn        
   CADD                      Float 
   CADD_SCALED               Float 
   GoNL_GTC                  String
   GoNL_AF                   String
   Thousand_Genomes_AF       String
   EXAC_AF                   String
   FITCON_SCORE              Float 
   DANN_SCORE                Float 
   CLINVAR_CLNSIG            String
   CLINVAR_CLNALLE           String
   GENE_SYMBOL               String",
  'out'     , 'o', 1, "character", "ouput CSV path + filename of annotated vcf-file analysis",
  'plot'    , 'p', 1, "character", "T for plot output and F for no plot output, defaults to false",
  'pout'    , 'q', 1, "character", "Ouput directory to write png graphs of the analysis, every gene will generate 1 plot. gene name will be used as file name",
  'help'   , 'h', 0, "logical",   "this help"
),ncol=5,byrow=T)

# parse commandline option
opt = getopt(spec);

# check all arguments are present or help option is present, if needed show usage
if(!is.null(opt$help) || is.null(opt$input) || is.null(opt$out)){
  cat(paste(getopt(spec, usage=T),"\n"));
  q()
}

# set default plot options if needed
if (is.null(opt$plot)){
  opt$plot <- F
}

# check for plot output dir
if (opt$plot){
  if(is.null(opt$pout)){
    cat(paste(getopt(spec, usage=T),"\n"));
    q()
  }
}

# start analysis
vcf.file<-readVcf(opt$input, "hg19")
genes <- unique(info(vcf.file)$GENE_SYMBOL)
print(paste(c("## starting analysis of: ", opt$input)))
results.cgd.df<-lapply(genes, analyse.Gene, vcf=vcf.file, plots=opt$plot, plotout=opt$pout)
df.results<-do.call("rbind", results.cgd.df)
print(paste(c(opt$out,".csv"), sep=""))
write.table(df.results, file = opt$out, sep = "\t", col.names = NA, qmethod = "double")

print("#### done######")
