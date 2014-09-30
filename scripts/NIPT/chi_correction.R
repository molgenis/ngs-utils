#!/usr/bin/env Rscript
DOC = "Calculate normalized chi square score per bin."

# Constants
chi.square.cut.off		= 3.5	# explanation needed
chromosomes.focus		= 1:22	# only autosomal chromosomes #--Risk-> What if order changes?
n.best.control.samples	= 50	# number of best control samples that is used for this analysis

# Retrieve command line parameters
suppressPackageStartupMessages(library("argparse"))

# Create parser object
parser = ArgumentParser(description = DOC)

# Define command line parameters
parser$add_argument("-n", "--name",				type = "character",		metavar = "sample name (id)", 		required = T,	help = "Sample's unique identifier.")
parser$add_argument("-i", "--input",			type = "character",		metavar = "binned data", 			required = T,	help = ".tsv file with number of reads per bin.")
parser$add_argument("-s", "--strand",			type = "character",		metavar = "strand", 				required = T,	help = "Strand, either \"forward\" or \"reverse\" ")
parser$add_argument("-d", "--controldir",		type = "character",		metavar = "control directory ",		required = T,	help = "Directory with control set.")
parser$add_argument("-c", "--controlsamples",	type = "character",		metavar = "best control samples",	required = T,	help = "List with control samples, descending quality.")
parser$add_argument("-w", "--workdir",			type = "character",		metavar = "work directory",			required = T,	help = "Intermediate results (needed by next steps in pipeline) are stored here.")
parser$add_argument("-v", "--verbose",			action="store_true",	default=TRUE,										help = "Shows with which parameters this script is called [default].")
parser$add_argument("-q", "--quietly",			action="store_false",	dest="verbose",										help = "Print little or no output.")

#
## Parse command line parameters
#
args = parser$parse_args()

# Proceed only if cmnd-line parameters correct
if (args$verbose) {
	write("You have used the following arguments:", stdout())
	write(paste("\t--input:  ", args$input), stdout())
	write(paste("\t--output: ", args$output), stdout())
	write(paste("\t--pdf:    ", args$pdf), stdout())
	write("\nStart with binning... (may take minutes)\n", stdout())
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
source(paste(LocationOfThisScript(), "chi_correction_functions.R", sep="/"))

#
## Load data
#
# Load binned data of sample of interest and select chromosomes on which we focus
sample.bins = read.delim(args$input, header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE)
sample.bins = as.matrix(sample.bins[chromosomes.focus, ])

# Load best control samples
control.file.base.name = as.vector(read.table(args$controlsamples)[1:n.best.control.samples, 1])

# Loads the control bins
control.bins = GetControlFiles(control.file.base.name = control.file.base.name, control.dir = args$controldir, strand = args$strand, chromosomes.focus = chromosomes.focus))

#
## Start calculations
#
# Calculate the chi square score per bin, based on only control samples
chi.sum.bins = SumChiScores(bins.list = control.bins)

# Append sample of interest to list with control samples. Next correct all samples in list based on 'chi square score' in control samples.
bins.list = control.bins
bins.list[[ 1 + length(bin.list) + 1 ]] = sample.bins

# Applies correction to bins
correct.bins.list = CorrectBins(bins.list = bins.list, chi.sum.bins = chi.sum.bins, strand = args$strand, sample.name = args$name)

# Save sample as .tsv
sample.bins = correct.bins.list[[1 + length(correct.bins.list)]]
write.table(sample.bins, paste(args$workdir, "/", args$name, ".", args$strand,  ".corrected.bins.table.tsv", sep=""), , quote = FALSE, sep ="\t", row.names = TRUE)

# Remove sample of interest from list
correct.bins.list[[-length(correct.bins.list)]]

# Compose file names of control samples for storage
control.sample.files = sub(".bins.", ".corrected.bins.", args$controlsamples)

# Save control samples as RDS
saveRDS(list(correct.bins.list[[-length(correct.bins.list)]], control.sample.files), paste(args$workdir, "/", args$name, ".", args$strand, ".controlfiles.corrected.bins.rds", sep=""))

# Quit with normal return code
quit(save = "no", status = 0)