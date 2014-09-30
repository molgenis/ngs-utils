# Return list with bins of best control files
getControlFiles = function(control.file.base.name, control.dir, strand, chromosomes.focus)
{
	files = paste(control.dir, "/", file.base.name, ".", strand, ".bins.table.tsv", sep = "")
	
	control.file.bin = list()
	for (i in 1:length(files) ){
		binfile = as.matrix(read.delim(files[i], header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE))

		# Only select relevant chromosomes
		control.file.bin[[i]] = binfile[chromosomes.focus, ]
	}
	
	return (control.file.bin)
}

sumChiScores = function(bins.list) {
	# Scale each sample in bins.list so that all samples have the same number of reads
	bins.overall.mean = sum(unlist(bins.list)) / length(bins.list)
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
correctBins = function(bins.list, chi.sum.bins, strand, sample.name) {
	# About the input parameters:
	# bin.list = list(control sample 1, control sample 2, ..., sample of interest)
	# chi.sum.bins is determined only on control samples
	
	degrees.of.freedom = n.best.control.samples - 1 # number of control samples minus one
	
	# Convert chi squares to a normal distribution score 
	chi.sum.bins.normalized = (chi.sum.bins - degrees.of.freedom) / (sqrt( 2 * degrees.of.freedom))

	# Derive correction factor for the read counts
	chi.sum.bins.correction.factor = chi.sum.bins / degrees.of.freedom

	# For each sample, correct extreme bins with high chi-square value
	for (i in 1:length(bins.list)) {
		m = bins.list[[i]]
		index = which(chi.square.cut.off < chi.sum.bins.normalized) # Variation between these bins is considered too large
		m[index] = m[index] / chi.sum.bins.correction.factor[index] # Therefore, we correct them here
		bins.list[[i]] = m
	}

	# Return corrected bin.list
	return(bins.list) 
}