#!/usr/bin/env Rscript
DOC = "Produces four different models. Each model predicts the fraction of reads mapping on the chromosome of interest. Each model contains a different number (based on best F-regression score) of predictors to do so. Predictors are unique within and between models."

# Constants
chromosomes.trisomy    = c(13, 18, 21, 35, 40, 43)      						# chromosomes with potential trisomy
chromosomes.background	= 1:44											# 'background' chromosomes used to predict trisomy
chromosomes.background.single <- 1:22
control.chromosomes 	= chromosomes.background[-chromosomes.trisomy]	# chroms we want to use to predict chromosomes with potential trisomy
n.models				= 4												# Number of models
max.predictors = 2
fixed.predictors = 4
chr.of.interest <- list(13,18,21)
max.elements <- 19
prediction.methods <- c("Fixed number of predictors, forward regression approach",
                        "Variable number of predictors, forward regression approach")
control.chromosomes.single = chromosomes.background.single[-unlist(chr.of.interest)]
chi.square.cut.off    = 3.5	# explanation needed

# Retrieve command line parameters
suppressPackageStartupMessages(library("argparser"))
suppressPackageStartupMessages(library("sets"))

# Create parser object
parser <- arg.parser(DOC, name="Create models absed on maximum F statistic score")
parser <- add.argument(parser, "-temp",       help = "Temp directory where files produced during pipeline run are stored")
parser <- add.argument(parser, "-d",       help = "Ground directory where control group .rds files are stored")
parser <- add.argument(parser, "-ds",       help = "Ground directory where sample rds files are stored")
parser <- add.argument(parser, "-type", help = "Type of data, for instance Illumina_QNA_Peak_Bin, Solid_Loess etc")
parser <- add.argument(parser, "-bf", help="Best fit variable. selects how many samples to use in best fit. Set to 0 to not use best fit at all")
parser <- add.argument(parser, "-chi",       help = "Boolean. set to TRUE to use chi square correction. Set to FALSE to omit chi square correction")

args <- parse.args(parser, argv = commandArgs(trailingOnly = TRUE))

source("/gcc/groups/gcc/tmp01/ddeweerd/research_scripts/prediction_scoring_general_functions.R")
source("/gcc/groups/gcc/tmp01/ddeweerd/research_scripts/bonferonni_corrected_functions.R")
source("/gcc/groups/gcc/tmp01/ddeweerd/research_scripts/regression_functions.R")
source("/gcc/groups/gcc/tmp01/ddeweerd/research_scripts/NCV_functions.R")
source("/gcc/groups/gcc/tmp01/ddeweerd/research_scripts/mad_functions.R")
source("/gcc/groups/gcc/tmp01/ddeweerd/research_scripts/z_score_functions.R")
source("/gcc/groups/gcc/tmp01/ddeweerd/research_scripts/ChiCorrectionPipeline.R")
source("/gcc/groups/gcc/tmp01/ddeweerd/research_scripts/SumOfSquares.R")

cat("Currently analyzing", args$type, "\n")
#Load data
forward.data <- readRDS(paste(args$d, "/", args$type, "_Forward.rds", sep=""))
reverse.data <- readRDS(paste(args$d, "/", args$type, "_Reverse.rds", sep=""))

bins.forward.all = forward.data[[1]]
bins.reverse.all = reverse.data[[1]]
control.file.names.all <- forward.data[[2]]

forward.samples <- readRDS(paste(args$ds, "/", "Positive_Control_", args$type, "_Forward.rds", sep=""))
reverse.samples <- readRDS(paste(args$ds, "/", "Positive_Control_", args$type, "_Reverse.rds", sep=""))

names.forward <- forward.samples[[2]]
names.reverse <- reverse.samples[[2]]

files.forward <- forward.samples[[1]]
files.reverse <- reverse.samples[[1]]

best.fit <- as.numeric(args$bf)
best.fit.used <- "All"
if (args$bf > 0)
{
  best.fit.used <- paste("BestFit", best.fit, sep="")
}
chi <- "Uncorrected"
if (args$chi == TRUE)
{
  chi <- "Chi_Corrected"
}

output.dir <- paste(args$temp, "Output_", args$type, "_", best.fit.used,"_", chi, sep="")
dir.create(path = output.dir, recursive = TRUE)
setwd(output.dir)

for (current.sample in 1:length(names.forward))
{
  current.file.name <- sub("^([^.]*).*", "\\1",names.forward[current.sample])
  grep(pattern = current.file.name, x = control.file.names.all)
  sample.bins.f <- files.forward[[current.sample]]
  sample.bins.r <- files.reverse[[current.sample]]
  
  sample.bins.f[which(is.na(sample.bins.f))] = 0
  sample.bins.r[which(is.na(sample.bins.r))] = 0
  
  bins.forward <- bins.forward.all
  bins.reverse <- bins.reverse.all
  control.file.names <- control.file.names.all
  
  #Matrix with chromosomal fractions per control sample
  chromosomal.frac.control <- ChromosomalFractionPerSample(control.file.names, bins.forward, bins.reverse)
  
  sample.of.interest <- data.frame(ChromosomalFractionPerSample("sample of interest", list(sample.bins.f), list(sample.bins.r)))
  
  if (args$bf > 0)
  {
    best.fit.samples <- BestControlSet(chromosomal.fractions = chromosomal.frac.control, samplefractions = sample.of.interest, n.of.samples = args$bf)
    bins.forward <- bins.forward[best.fit.samples]
    bins.reverse <- bins.reverse[best.fit.samples]
    control.file.names<- control.file.names[best.fit.samples]
  }
  
  if (args$chi == TRUE)
  {
    chi.result <- ChiCorrect(sample.forward = sample.bins.f, sample.reverse = sample.bins.r, control.bins.forward = bins.forward, control.bins.reverse = bins.reverse)
    bins.forward <- chi.result[[1]]
    sample.bins.f <- chi.result[[2]]
    bins.reverse <- chi.result[[3]]
    sample.bins.r <- chi.result[[4]]
  }
  
  chromosomal.frac.control <- ChromosomalFractionPerSample(control.file.names, bins.forward, bins.reverse)
  sample.of.interest <- data.frame(ChromosomalFractionPerSample("sample of interest", list(sample.bins.f), list(sample.bins.r)))
  
  sample.bins <- sample.bins.f + sample.bins.r
  
  col.names <- colnames(chromosomal.frac.control)
  
  #Determine predictors
  predictors.13 <- sapply(chr.of.interest[[1]], FUN=GetPredictors, chromosomal.frac.control=chromosomal.frac.control)
  predictors.18 <- sapply(chr.of.interest[[2]], FUN=GetPredictors, chromosomal.frac.control=chromosomal.frac.control)
  predictors.21 <- sapply(chr.of.interest[[3]], FUN=GetPredictors, chromosomal.frac.control=chromosomal.frac.control)
  
  chromosomal.frac.control.reads <- ChromosomalReadsPerSample(control.file.names, bins.forward, bins.reverse)
  #Determine Denominators 
  denominator.list <- lapply(chr.of.interest, FUN=GetDenominators, chromosomal.frac.control.reads=chromosomal.frac.control.reads)
  
  output.name = sub(".bins.table.tsv", ".final.output.csv", names.forward[current.sample])
  output.name <- sub(pattern = ".forward", replacement = "", x = output.name)
  
  output.bon <- sub(pattern = ".final.", replacement = ".bonferroni.", x = output.name)
  output.file <- output.name
  output.file.bon <- output.bon
  
  #Predict trisomy regression based
  for( prediction in 1:2)
  {
    AppendLine(line = c("Chromosome of interest: 13", prediction.methods[prediction]), file = output.file)
    AppendLine(line = c("Predictors", "Z score", "VC", "Shapiro"), file = output.file)
    lapply(predictors.13[[prediction]], FUN=PredictTrisomy, chr.int=13, 
           control.samples = data.frame(chromosomal.frac.control), sample.bins.f=sample.bins.f, sample.bins.r=sample.bins.r,
           sample.of.interest = sample.of.interest)
    if (length(predictors.13[[prediction]]) < 4)
    {
      times <- 4 - length(predictors.13[[prediction]])
      for (blank.line in 1:times)
      {
        AppendLine(line = "", file = output.file)
      }
    }
    AppendLine(line = c("Chromosome of interest: 18", prediction.methods[prediction]), file = output.file)
    AppendLine(line = c("Predictors", "Z score", "VC", "Shapiro"), file = output.file)
    lapply(predictors.18[[prediction]], FUN=PredictTrisomy, chr.int=18, 
           control.samples = data.frame(chromosomal.frac.control), sample.bins.f=sample.bins.f, sample.bins.r=sample.bins.r,
           sample.of.interest = sample.of.interest)
    if (length(predictors.18[[prediction]]) < 4)
    {
      times <- 4 - length(predictors.18[[prediction]])
      for (blank.line in 1:times)
      {
        AppendLine(line = "", file = output.file)
      }
    }
    AppendLine(line = c("Chromosome of interest: 21", prediction.methods[prediction]), file = output.file)
    AppendLine(line = c("Predictors", "Z score", "VC", "Shapiro"), file = output.file)
    lapply(predictors.21[[prediction]], FUN=PredictTrisomy, chr.int=21,
           control.samples = data.frame(chromosomal.frac.control), sample.bins.f=sample.bins.f, sample.bins.r=sample.bins.r,
           sample.of.interest = sample.of.interest)
    if (length(predictors.21[[prediction]]) < 4)
    {
      times <- 4 - length(predictors.21[[prediction]])
      for (blank.line in 1:times)
      {
        AppendLine(line = "", file = output.file)
      }
    }
  }
  AppendLine(line=" ", file = output.file)
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
  AppendLine(line = c("Chromosome", "Number of significant calls", "maximum number of calls", "Mean corrected P value"), file = output.file)
  #Bonferronni corrected Z scores
  BonferonniCorrected(sample.bins)
}