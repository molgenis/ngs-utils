#!/usr/bin/env Rscript
DOC = "Correct bins with high chi square score."

# Constants 
chi.square.cut.off  	= 3.5	# explanation needed
chromosomes.focus		= 1:22	# only autosomal chromosomes #--Risk-> What if order changes?


# Retrieve command line parameters
suppressPackageStartupMessages(library("argparser"))

# Create parser object
parser <- arg.parser(DOC, name="Chi Square correction")

# Define command line parameters
parser <- add.argument(parser, "-f",  help = ".tsv file with number of reads per bin forward strand.")
parser <- add.argument(parser, "-r",  help = ".tsv file with number of reads per bin reverse strand.")
parser <- add.argument(parser, "-temp",       help = "Temp directory where files produced during pipeline run are stored")
parser <- add.argument(parser, "-sample", 	help = "Sample ID ")
parser <- add.argument(parser, "-d", help = "Directory with control set.")
parser <- add.argument(parser, "-q", 		help = "Print little or no output.")

args <- parse.args(parser, argv = commandArgs(trailingOnly = TRUE))

GetFiles <- function(controlDir, strand, sampleID)
{
  setwd(controlDir)
  files <- list.files( pattern = paste("*", strand,".bins.table.tsv", sep = ""))
  sampleID <- paste(sampleID,".", strand, ".bins.table.tsv", sep ="")
  files <- files[!files %in% sampleID]
  return (sort(files))
}

# Return list with bins of best control files
GetControlFiles = function(control.file.base.name, control.dir, strand, chromosomes.focus)
{
  control.file.bin = list()
  for (i in 1:length(control.file.base.name) ){
    binfile = as.matrix(read.delim(control.file.base.name[i], header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE))
    
    # Only select relevant chromosomes
    control.file.bin[[i]] = binfile[chromosomes.focus, ]
  }
  
  return (control.file.bin)
}


SumChiScores = function(bins.list) {
  # Scale each sample in bins.list so that all samples have the same number of reads
  bins.overall.mean = sum(as.numeric(unlist(bins.list))) / length(bins.list)
  mean.correction = function(m) m * bins.overall.mean / sum(m)
  bins.list.scaled = lapply(bins.list, mean.correction)
  
  # Calculate the expected value (= mean) per bin
  bins.sum.corrected = Reduce("+", bins.list.scaled)
  bins.scaled.expected = bins.sum.corrected / length(bins.list)
  
  # Calculate chi-score per bin per sample
  bins.chi.score = bins.list.scaled
  for (i in 1:length(bins.list)) {
    bins.chi.score[[i]] = (bins.scaled.expected - bins.chi.score[[i]])^2 / bins.scaled.expected
  }
  
  # Return sum of 
  return(Reduce("+", bins.chi.score))
}

# Corrects overdispersed bins in bins.list
CorrectBins = function(bins.list, chi.sum.bins, strand, sample.name) {
  # About the input parameters:
  # bin.list = list(control sample 1, control sample 2, ..., sample of interest)
  # chi.sum.bins is determined only on control samples
  
  degrees.of.freedom = n.best.control.samples- 1 # number of control samples minus one
  
  # Convert chi squares to a normal distribution score 
  chi.sum.bins.normalized = (chi.sum.bins - degrees.of.freedom) / (sqrt( 2 * degrees.of.freedom))
  
  # Derive correction factor for the read counts
  chi.sum.bins.correction.factor = chi.sum.bins / degrees.of.freedom
  
  # For each sample, correct extreme bins with high chi-square value
  print(length(bins.list))
  for (i in 1:length(bins.list)) {
    m = bins.list[[i]]
    index = which(chi.square.cut.off < chi.sum.bins.normalized) # Variation between these bins is considered too large
    m[index] = m[index] / chi.sum.bins.correction.factor[index] # Therefore, we correct them here
    bins.list[[i]] = m
  }
  
  # Return corrected bin.list
  return(bins.list) 
}
#
## Parse command line parameters
#
#args = parser$parse_args()

# Proceed only if cmnd-line parameters correct
#if (args$verbose) {
 # write("You have used the following arguments:", stdout())
  #write(paste("\t--input:  ", args$input), stdout())
  #write(paste("\t--output: ", args$output), stdout())
  #write(paste("\t--pdf:    ", args$pdf), stdout())
  #write("\nStart with binning... (may take minutes)\n", stdout())
#}

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
print(LocationOfThisScript())
# Load functions
#source(paste(LocationOfThisScript(), "chi_correction_functions.R", sep="/"))

#
## Load data
#
# Load binned data of sample of interest and select chromosomes on which we focus
sample.bins.forward = read.delim(args$f, header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE)
sample.bins.forward = as.matrix(sample.bins.forward[chromosomes.focus, ])

sample.bins.reverse = read.delim(args$r, header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE)
sample.bins.reverse = as.matrix(sample.bins.reverse[chromosomes.focus, ])
# Load best control samples
control.file.base.name.forward = GetFiles(args$d, strand="forward", args$sample)
control.file.base.name.reverse = GetFiles(args$d, strand="reverse", args$sample)

n.best.control.samples  = length(control.file.base.name.forward)

# Loads the control bins
control.bins.forward = GetControlFiles(control.file.base.name.forward , control.dir = args$d, chromosomes.focus = chromosomes.focus)
control.bins.reverse = GetControlFiles(control.file.base.name.reverse , control.dir = args$d, chromosomes.focus = chromosomes.focus)

control.bins <- NULL
for (p in 1:length(control.bins.forward))
{
  control.bins[[p]] <- control.bins.forward[[p]] + control.bins.reverse[[p]]
}

print(length(control.bins))
#
## Start calculations
#
# Calculate the chi square score per bin, based on only control samples
chi.sum.bins = SumChiScores(bins.list = control.bins)

# Append sample of interest to list with control samples. Next correct all samples in list based on 'chi square score' in control samples.
bins.list.forward = control.bins.forward
bins.list.forward[[length(bins.list.forward) + 1 ]] = sample.bins.forward

bins.list.reverse = control.bins.reverse
bins.list.reverse[[length(bins.list.reverse) + 1 ]] = sample.bins.reverse
# Applies correction to bins
correct.bins.list.forward = CorrectBins(bins.list = bins.list.forward, chi.sum.bins = chi.sum.bins, strand = args$strand, sample.name = args$sample)
correct.bins.list.reverse = CorrectBins(bins.list = bins.list.reverse, chi.sum.bins = chi.sum.bins, strand = args$strand, sample.name = args$sample)
# Save sample as .tsv
# TODO ---- do not construct file name here, but compose in parameters
sample.bins.forward = correct.bins.list.forward[[length(correct.bins.list.forward)]]
write.table(sample.bins.forward, paste(args$temp, "/", args$sample, ".forward.corrected.bins.table.tsv", sep=""), quote = FALSE, sep ="\t", row.names = TRUE)

sample.bins.reverse = correct.bins.list.reverse[[length(correct.bins.list.reverse)]]
write.table(sample.bins.reverse, paste(args$temp, "/", args$sample, ".reverse.corrected.bins.table.tsv", sep=""), quote = FALSE, sep ="\t", row.names = TRUE)
# Remove sample of interest from list
correct.bins.list.forward[length(correct.bins.list.forward)] <- NULL
correct.bins.list.reverse[length(correct.bins.list.reverse)] <- NULL

# Compose file names of control samples for storage
control.sample.files.forward = sub(".bins.", ".corrected.bins.", control.file.base.name.forward)
control.sample.files.reverse = sub(".bins.", ".corrected.bins.", control.file.base.name.reverse)
# Save control samples as RDS
saveRDS(list(correct.bins.list.forward, control.sample.files.forward), paste(args$temp, "/", args$sample, ".forward.controlfiles.corrected.bins.rds", sep=""))
saveRDS(list(correct.bins.list.reverse, control.sample.files.reverse), paste(args$temp, "/", args$sample, ".reverse.controlfiles.corrected.bins.rds", sep=""))

# Quit with normal return code
quit(save = "no", status = 0)