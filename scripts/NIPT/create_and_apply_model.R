#!/usr/bin/env Rscript
DOC = "Produces four different models. Each model predicts the fraction of reads mapping on the chromosome of interest. Each model contains four different predictors (e.g. 1F, 14R, 3F, 22F) to do so. Predictors are unique within and between models."

# Constants
chromosomes.trisomy		= c(13, 18, 21)									# chromosomes with potential trisomy
chromosomes.background	= 1:22											# 'background' chromosomes used to predict trisomy
control.chromosomes 	= chromosomes.background[-chromosomes.trisomy]	# chroms we want to use to predict chromosomes with potential trisomy
n.models				= 4												# Number of models
n.predictors.per.model	= 4												# Number of predictors per model
round.n.decimals		= 3												# Number of digits on which results are rounded
cv.correction.factor	= 1.15											# (Arbitrary?) Factor inflates theoretical coefficient of variation
fetal.dna.perc.min		= 1												# Percentage of reads originating from fetus (under limit)
fetal.dna.perc.max		= 30											# Percentage of reads originating from fetus (upper limit)
p.shapiro.cutoff		= 0.05											# Below cutoff we assume data is not normally distributed; such data are ignored

# Retrieve command line parameters
suppressPackageStartupMessages(library("argparse"))

# Create parser object
parser = ArgumentParser(description = DOC)

# Define command line parameters
parser$add_argument("-n", "--chromosome",		type = "character",		metavar = "chromosome of interest",	required = T,	help = "The chromosome you want to predict.")
parser$add_argument("-s", "--strand",			type = "character",		metavar = "strand", 				required = T,	help = "Strand, either \"forward\" or \"reverse\" ")
parser$add_argument("-n", "--name",				type = "character",		metavar = "sample name (id)", 		required = T,	help = "Sample's unique identifier.")
parser$add_argument("-w", "--workdir",			type = "character",		metavar = "work directory",			required = T,	help = "Intermediate results (needed by next steps in pipeline) are stored here.")
parser$add_argument("-m", "--modeltable",		type = "character",		metavar = "info on model",	 		required = T,	help = "Output file with information on prediction models.")
parser$add_argument("-r", "--risktable",		type = "character",		metavar = "risk table",		 		required = T,	help = "Output file with results (incl. risk on trisomy) for this chromosome.")
parser$add_argument("-v", "--verbose",			action="store_true",	default=TRUE,										help = "Shows with which parameters this script is called [default].")
parser$add_argument("-q", "--quietly",			action="store_false",	dest="verbose",										help = "Print little or no output.")

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

# TODO ---- File names in parameters!
# TODO ---- You load list with 2 elements. Please save/load these individually.
# Load Chi^2 corrected files and file list
forward.data = readRDS(paste(args$workdir, "/", args$name,".forward.controlfiles.corrected.bins.rds", sep=""))
reverse.data = readRDS(paste(args$workdir, "/", args$name,".reverse.controlfiles.corrected.bins.rds", sep=""))

# Gets the Chi2 corrected files
bins.forward = forward.data[[1]]
bins.reverse = reverse.data[[1]]
control.file.names = forward.data[[2]]

# Get 'chromosomal fractions'
# chr.frac[sample, chrom] is the ratio '#reads on chrom-direction' / '#reads in sample-direction', where chrom is 1F, 2F, ... nF, 1R, 2R, ... or nR
# TODO ---- is this correct? I can imagine that we want '#reads on chrom-direction' / '#reads in sample-TOTAL'?
chr.frac = ChromosomalFractionPerSample(control.file.names, bins.forward, bins.reverse)

# Get #reads on chromosome with potential trisomy for each of the samples
frac.reads.chr.trisomy.observed = chr.frac[args$chromosome,] + chr.frac[args$chromosome + length(chromosomes.background),]

# Remove chromosomes that may have a trisomy so that they don't interfere with analysis
chr.frac = chr.frac[, -c(chromosomes.trisomy, chromosomes.trisomy + length(chromosomes.background))] # Remove {13, 18, 21}F and {13, 18, 21}R, respectively.

# Determine the n.models models with n.predictors.per.model predictors each (all predictors are 'unique' within and between the models)
model.list		= NULL # the models
chr.selected	= NULL
for (i.model in 1:n.models)
{
	# Select the best n.predictors.per.model predictors
	for (i.predictor in 1:n.predictors.per.model)
	{
		chr.selected = c(chr.selected, BestPredictor(frac.reads.chr.trisomy.observed, chr.frac, chr.selected, i.predictor))
	}
	
	# Add model to list
	model.list[i.model]	= lm(frac.reads.chr.trisomy.observed ~ chr.frac[, tail(chr.selected, n.predictors.per.model)])
}

#
## Load sample and determine the minimum value of the Coefficient of Variation, based on theory (multiplied with some factor)
#
# TODO ---- do not (re)construct file name here, but compose in parameters
sample.bins = read.table(paste(args$workdir, "/", args$name, ".", args$strand,  ".corrected.bins.table.tsv", sep=""), sep ="\t", row.names = TRUE)

# Get (chi-corrected) number of reads on chromosome of interest in sample of interest
n.reads.sample.chromosome = sample.bins[args$chromosome,] + sample.bins[args$chromosome + length(chromosomes.background),]

# TODO ---- Describe background! 
# I guess you mean sd.observed and sd.theoretical here (see below)
cv.theoretical = cv.correction.factor / sqrt(length(n.reads.sample.chromosome))

#
## Apply model and generate results
#
z.score.model							= NULL # z-score[[i.model]] is vector with z-scores for each of the control samples when using model i.model
p.shapiro.model							= NULL # p-value assuming data normally distributed, per model
trisomy.overrepresentation.mean.model	= NULL # mean of overrepresentation, per model # TODO ---- Should be 1? If so, remove
cv.model								= NULL # CV per model
for (i.model in 1:n.models)
{
	frac.predicted = predict(model.list[i.model], chr.frac)
	
	# number of times we observe what was expected
	trisomy.overrepresentation = frac.reads.chr.trisomy.observed / frac.predicted
	trisomy.overrepresentation.mean.model[i.model] = mean(trisomy.overrepresentation)

	# See http://en.wikipedia.org/wiki/Coefficient_of_variation
	# Z-score has sd in denominator, not CV, so I guess you mean sd.observed and sd.theoretical here
	# This can be correct, because mean(trisomy.overrepresentation) should be 1?
	cv.observed = (1 + 1 / (4 * length(trisomy.overrepresentation))) * sd(trisomy.overrepresentation) / mean(trisomy.overrepresentation)
	
	cv = max(cv.theoretical, cv.observed)
	cv.model[i.model] = cv
	
	z.score.model[[i.model]]	= (trisomy.overrepresentation - mean(trisomy.overrepresentation)) / cv
	p.shapiro.model[[i.model]]	= shapiro.test(z.score.model[[i.model]])
}

#
## Calculate risk on trisomy for chromosome of interest in our sample, given each of the models
#
risk = NULL # Predicted risk on trisomy, per model
for (i.model in 1:n.models)
{
	frac.predicted = predict(model.list[i.model], sample.bins)

	# Get observed fraction of reads on this chromosome compared to total number of reads
	frac.reads.chr.trisomy.observed = n.reads.sample.chromosome / sum(sample.bins)

	# number of times we observe what was expected in sample
	trisomy.overrepresentation.sample = frac.reads.chr.trisomy.observed / frac.predicted

	# z-score for sample
	z.score.sample = (trisomy.overrepresentation.sample - trisomy.overrepresentation.mean.model[i.model]) / cv.model[i.model]

	risk[i.model] = Risk(fetal.dna.perc.min, fetal.dna.perc.max, cv.model[i.model], z.score.sample, a.priori)
}

# Determine resulting risk, based on each of the models that have normally distributed z-scores
normally.distributed.z.scores = p.shapiro.cutoff < p.shapiro.model
risk.result = median( risk[which(normally.distributed.z.scores)] )

#
## Save results
#

# Save table with information about the models
	model.table = GetModelInfoTable(model.list)
	write.table(model.table, args$modeltable, quote = FALSE, sep ="\t", row.names = TRUE, col.names = TRUE)

# Save table with predicted risks
	risk.table = GetRiskInfoTable(cv.model, normally.distributed.z.scores, z.score.sample, risk, risk.result)
	write.table(risk.table, args$risktable, quote = FALSE, sep ="\t", row.names = TRUE, col.names = TRUE)