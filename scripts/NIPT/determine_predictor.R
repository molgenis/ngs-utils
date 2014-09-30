#!/usr/bin/env Rscript
DOC = "Produces four different models. Each model predicts the fraction of reads mapping on the chromosome of interest. Each model contains four different predictors (e.g. 1F, 14R, 3F, 22F) to do so. Predictors are unique within and between models."

# Constants
chromosomes.trisomy		= c(13, 18, 21)									# chromosomes with potential trisomy
chromosomes.background	= 1:22											# 'background' chromosomes used to predict trisomy
control.chromosomes 	= chromosomes.background[-chromosomes.trisomy]	# chroms we want to use to predict chromosomes with potential trisomy
n.best.control.samples	= 50											# number of best control samples that is used for this analysis
n.models				= 4												# Number of models
n.predictors.per.model	= 4												# Number of predictors per model
round.n.decimals		= 3												# Number of digits on which results are rounded

# Retrieve command line parameters
suppressPackageStartupMessages(library("argparse"))

# Create parser object
parser = ArgumentParser(description = DOC)

# Define command line parameters
parser$add_argument("-n", "--chromosome",		type = "character",		metavar = "chromosome of interest",	required = T,	help = "The chromosome you want to predict.")
parser$add_argument("-n", "--name",				type = "character",		metavar = "sample name (id)", 		required = T,	help = "Sample's unique identifier.")
parser$add_argument("-w", "--workdir",			type = "character",		metavar = "work directory",			required = T,	help = "Intermediate results (needed by next steps in pipeline) are stored here.")
parser$add_argument("-i", "--models",			type = "character",		metavar = "prediction model", 		required = T,	help = "Output file with prediction models.")
parser$add_argument("-s", "--fractiontable",	type = "character",		metavar = "fraction table", 		required = T,	help = "Output file with per sample the fractions of each chromosome.")
parser$add_argument("-v", "--verbose",			action="store_true",	default=TRUE,									help = "Shows with which parameters this script is called [default].")
parser$add_argument("-q", "--quietly",			action="store_false",	dest="verbose",									help = "Print little or no output.")

#
## Parse command line parameters
#
args = parser$parse_args()

# Proceed only if cmnd-line parameters correct
if (args$verbose) {
	write("You have used the following arguments:", stdout())
	for (i in 1:length(args)) write(paste('--', names(args[i]), ": ", args[[i]], sep=''), stdout())
		
	write("\nStarting analysis... (may take minutes)\n", stdout())
}

#
## Load functions
#
# First determine the location of this R-script
LocationOfThisScript = function() # Function LocationOfThisScript returns the location of this .R script (may be needed to source other files in same dir)
{
	this.file = NULL
	# This file may be 'sourced'
	for (i in -(1:sys.nframe())) {
		if (identical(sys.function(i), base::source)) this.file = (normalizePath(sys.frame(i)$ofile))
	}

	if (!is.null(this.file)) return(dirname(this.file))

	# But it may also be called from the command line
	cmd.args = commandArgs(trailingOnly = FALSE)
	cmd.args.trailing = commandArgs(trailingOnly = TRUE)
	cmd.args = cmd.args[seq.int(from=1, length.out=length(cmd.args) - length(cmd.args.trailing))]
	res = gsub("^(?:--file=(.*)|.*)$", "\\1", cmd.args)

	# If multiple --file arguments are given, R uses the last one
	res = tail(res[res != ""], 1)
	if (0 < length(res)) return(dirname(res))

	# Both are not the case. Maybe we are in an R GUI?
	return(NULL)
}

# Load functions
source(paste(LocationOfThisScript(), "determine_predictor_functions.R", sep="/"))

# ---- File names in parameters!
# ---- You load list with 2 elements. Please save/load these individually.
# Load Chi^2 corrected files and file list
forward.data = readRDS(paste(args$workdir, "/", args$name,".forward.controlfiles.corrected.bins.rds", sep=""))
reverse.data = readRDS(paste(args$workdir, "/", args$name,".reverse.controlfiles.corrected.bins.rds", sep=""))

# Gets the Chi2 corrected files
bins.forward = forward.data[[1]]
bins.reverse = reverse.data[[1]]
control.file.names = forward.data[[2]]

# Get 'chromosomal fractions'
# chr.frac[sample, chrom] is the ratio '#reads on chrom-direction' / '#reads in sample-direction', where chrom is 1F, 2F, ... nF, 1R, 2R, ... or nR
# ---- is this correct? I can imagine that we want '#reads on chrom-direction' / '#reads in sample-TOTAL'?
chr.frac = ChromosomalFractionPerSample(control.file.names, bins.forward, bins.reverse)

# Save 'chromosomal fractions' (transposed)
# ---- Can't we just save chr.frac 'as is' (--> adapt later steps)
write.table(t(chr.frac), args$fractiontable, quote = FALSE, sep ="\t", row.names = TRUE, col.names = TRUE)

# Get #reads on chromosome with potential trisomy for each of the samples
n.reads.chr.trisomy = chr.frac[args$chromosome,] + chr.frac[args$chromosome + length(chromosomes.background),]

# Remove chromosomes that may have a trisomy so that they don't interfere with analysis
chr.frac = chr.frac[, -c(chromosomes.trisomy, chromosomes.trisomy + length(chromosomes.background))] # Remove {13, 18, 21}F and {13, 18, 21}R, respectively.

# Determine the n.models models with n.predictors.per.model predictors each (all predictors are 'unique' within and between the models)
model.list		= NULL # the models
predictor.list	= NULL # the predictors per model
for (i.model in 1:n.models)
{
	# Select the best n.predictors.per.model predictors
	chr.selected = NULL
	for (i.predictor in 1:n.predictors.per.model)
	{
		chr.selected = c(chr.selected, BestPredictor(n.reads.chr.trisomy, chr.frac, chr.selected))
	}
	
	model.list[i.model]		= lm(n.reads.chr.trisomy ~ chr.frac[, tail(chr.selected, n.predictors.per.model)])
	predictor.list[i.model]	= chr.selected
}

# Compose matrix with results we want to save
results = NULL
for (i.model in 1:n.models)
{
	m = model.list[i.model]

	# Predictors
	preds			= predictor.list[i.model]
	names(preds)	= paste("Pred", 1:n.predictors.per.model, sep="")
	
	# Coefficients
	coefs			= m$coef
	names(coefs)	= c("Intercept", paste("Coef", 2:n.predictors.per.model))
	
	# Custom round function
	Round = function(x) round(x, round.n.decimals)
		
	result = c(preds, Round(coefs), Round(unlist(summary(lm(as.data.frame(mat)))[c("r.squared", "adj.r.squared", "fstatistic", "sigma")])))
	
	# Add new row with 
	results = rbind(results, result)
}
rownames(results) = paste("Model", 1:n.models)

# Save results
# --- File names need to be composed outside of script (in parameters.csv or workflow.csv)
write.table(as.data.frame(results), paste( args$chromosome, args$models,  sep=""), quote = FALSE, sep ="\t", row.names = TRUE, col.names = TRUE)