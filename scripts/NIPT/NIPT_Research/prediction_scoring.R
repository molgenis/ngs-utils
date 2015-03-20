#!/usr/bin/env Rscript
DOC = "Produces four different models. Each model predicts the fraction of reads mapping on the chromosome of interest. Each model contains a different number (based on best F-regression score) of predictors to do so. Predictors are unique within and between models."

# Constants
chromosomes.trisomy    = c(13, 18, 21, 35, 40, 43)      						# chromosomes with potential trisomy
chromosomes.background	= 1:44											# 'background' chromosomes used to predict trisomy
chromosomes.background.single <- 1:22
control.chromosomes 	= chromosomes.background[-chromosomes.trisomy]	# chroms we want to use to predict chromosomes with potential trisomy
n.models				= 4												# Number of models
max.predictors = 2
chr.of.interest <- list(13,18,21)
max.elements <- 4
prediction.methods <- c("Fixed number of predictors, brute force approach",
                        "Variable number of predictors, brute force approach",
                        "Fixed number of predictors, forward regression approach",
                        "Variable number of predictors, forward regression approach")
control.chromosomes.single = chromosomes.background.single[-unlist(chr.of.interest)]	



# Retrieve command line parameters
suppressPackageStartupMessages(library("argparser"))
suppressPackageStartupMessages(library("sets"))

# Create parser object
parser <- arg.parser(DOC, name="Create models absed on maximum F statistic score")
parser <- add.argument(parser, "-f",  help = ".RDS file with control file bin counts, forward strand.")
parser <- add.argument(parser, "-r",  help = ".RDS file with control file bin counts, reverse strand.")
parser <- add.argument(parser, "-sf",  help = ".RDS output object containing brute force approach predictor sets, variable numer of predictors")
parser <- add.argument(parser, "-sr",  help = ".RDS output object containing brute force approach predictor sets, fixed number of predictors")
parser <- add.argument(parser, "-rv",  help = ".RDS output object containing forward regression approach predictor sets, variable number of predictors")
parser <- add.argument(parser, "-rf",  help = ".RDS output object containing forward regression approach predictor sets, fixed number of predictors")
parser <- add.argument(parser, "-temp",       help = "Temp directory where files produced during pipeline run are stored")
parser <- add.argument(parser, "-o", help = "Output file, CSV format")
parser <- add.argument(parser, "-ob", help = "Output file for bonferronni corrected Z scores, CSV format")

args <- parse.args(parser, argv = commandArgs(trailingOnly = TRUE))

cat("LOADING GENERAL FUNCTIONS\n")
source("/Users/dirkdeweerd/Graduation/Scripts/prediction_scoring_general_functions.R")

output.file <- paste(args$temp,args$o, sep="")
output.file.bon <- paste(args$temp,args$ob, sep="")
getwd()
#Load data
forward.data <- readRDS(args$f)
reverse.data <- readRDS(args$r)

bins.forward = forward.data[[1]]
bins.reverse = reverse.data[[1]]
control.file.names <- forward.data[[2]]

sample.bins.f <- read.delim(args$sf)
sample.bins.r <- read.delim(args$sr)
sample.bins <- sample.bins.f + sample.bins.r

sample.of.interest <- data.frame(ChromosomalFractionPerSample("sample of interest", list(sample.bins.f), list(sample.bins.r)))
#Matrix with chromosomal fractions per control sample
chromosomal.frac.control <- ChromosomalFractionPerSample(control.file.names, bins.forward, bins.reverse)
#Colnames (chromosomes) from chromosomal matrix. Later used to convert col numbers to names
col.names <- colnames(chromosomal.frac.control)

cat("\nLOADING REGRESSION FUNCTIONS\n")
source("/Users/dirkdeweerd/Graduation/Scripts/regression_functions.R")
#Determine predictors
predictors.13 <- sapply(chr.of.interest[[1]], FUN=GetPredictors, chromosomal.frac.control=chromosomal.frac.control)
predictors.18 <- sapply(chr.of.interest[[2]], FUN=GetPredictors, chromosomal.frac.control=chromosomal.frac.control)
predictors.21 <- sapply(chr.of.interest[[3]], FUN=GetPredictors, chromosomal.frac.control=chromosomal.frac.control)
#Predict trisomy regression based
for( prediction in 1:n.models)
{
  cat(prediction.methods[prediction], "\n")
  AppendLine(line = c("Chromosome of interest: 13", prediction.methods[prediction]), file = output.file)
  AppendLine(line = c("Predictors", "Z score", "VC", "Shapiro"), file = output.file)
  lapply(predictors.13[[prediction]], FUN=PredictTrisomy, chr.int=13, 
         control.samples = data.frame(chromosomal.frac.control), sample.bins.f=sample.bins.f, sample.bins.r=sample.bins.r,
         sample.of.interest = sample.of.interest)
  AppendLine(line = c("Chromosome of interest: 18", prediction.methods[prediction]), file = output.file)
  AppendLine(line = c("Predictors", "Z score", "VC", "Shapiro"), file = output.file)
  lapply(predictors.18[[prediction]], FUN=PredictTrisomy, chr.int=18, 
         control.samples = data.frame(chromosomal.frac.control), sample.bins.f=sample.bins.f, sample.bins.r=sample.bins.r,
         sample.of.interest = sample.of.interest)
  AppendLine(line = c("Chromosome of interest: 21", prediction.methods[prediction]), file = output.file)
  AppendLine(line = c("Predictors", "Z score", "VC", "Shapiro"), file = output.file)
  lapply(predictors.21[[prediction]], FUN=PredictTrisomy, chr.int=21,
       control.samples = data.frame(chromosomal.frac.control), sample.bins.f=sample.bins.f, sample.bins.r=sample.bins.r,
       sample.of.interest = sample.of.interest)
}
AppendLine(line=" ", file = output.file)
cat("\n\nLOADING NCV FUNCTIONS\n")
source("/Users/dirkdeweerd/Graduation/Scripts/NCV_functions.R")
chromosomal.frac.control.reads <- ChromosomalReadsPerSample(control.file.names, bins.forward, bins.reverse)
#Determine Denominators 
denominator.list <- lapply(chr.of.interest, FUN=GetDenominators, chromosomal.frac.control.reads=chromosomal.frac.control.reads)

AppendLine(line = "NCV Scores:", file = output.file)
#Predict trisomy NCV
for(chr.trisomy in 1:length(chr.of.interest))
{
  AppendLine(line = paste("Chromosome of interest", chr.of.interest[[chr.trisomy]], collapse =" "), file = output.file)
  AppendLine(line = c("Denominators", "Z score", "VC", "Shapiro"), file = output.file)
  PredictTrisomyNCV(chr.of.interest=chr.of.interest[[chr.trisomy]], denominators = denominator.list[[chr.trisomy]],
                    chromosomal.frac.control.reads=chromosomal.frac.control.reads, sample.bins=sample.bins)
}
AppendLine(line=" ", file = output.file)
cat("\n\nLOADING NORMAL Z SCORE FUNCTIONS\n")
source("/Users/dirkdeweerd/Graduation/Scripts/z_score_functions.R")
AppendLine(line = "Normal Z Scores:", file = output.file)
#Determine 'normal' Z scores
for(z.score in 1:length(chr.of.interest))
{
  AppendLine(line = paste("Chromosome of interest", chr.of.interest[[z.score]], collapse =" "), file = output.file)
  AppendLine(line = c(" ", "Z score", "VC", "Shapiro"), file = output.file)
  ZscoreSample(chromosomal.frac.control = chromosomal.frac.control,chr.int = chr.of.interest[[z.score]], 
               sample.of.interest = sample.of.interest)
}
AppendLine(line=" ", file = output.file)
cat("\n\nLOADING MAD FUNCTIONS\n")
source("/Users/dirkdeweerd/Graduation/Scripts/mad_functions.R")
AppendLine("Mad Scores", file = output.file)
#Determine Median Absolute Deviation based Z -score
for(mad.score in 1:length(chr.of.interest))
{
 
  AppendLine(line = paste("Chromosome of interest", chr.of.interest[[mad.score]], collapse =" "), file = output.file)
  AppendLine(line = c(" ", "MAD score", "VC (MAD)", "Shapiro"), file = output.file)
  DetermineMAD(chromosomal.frac.control = chromosomal.frac.control,chr.int = chr.of.interest[[mad.score]], 
               sample.of.interest = sample.of.interest)
}
AppendLine(line=" ", file = output.file)
AppendLine("Bonferonni Corrected Scores", file = output.file)
cat("\n\nLOADING BONFERONNI CORRECTED FUNCTIONS\n")
source("/Users/dirkdeweerd/Graduation/Scripts/bonferonni_corrected_functions.R")
AppendLine(line = c("Chromosome", "Number of significant calls", "maximum number of calls", "Mean corrected P value"), file = output.file)
#Bonferronni corrected Z scores
BonferonniCorrected(sample.bins)