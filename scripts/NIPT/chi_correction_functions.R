# Return list with bins of best control files
get.control.files = function(control.file.base.name, control.dir, strand, chromosomes.focus)
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

sum.chi.scores = function(bins.list) {
	n = length(bins.list)

	bins.list.mean = sum(unlist(bins.list)) / n
	
	mean.correction = function(m) m * bins.list.mean / sum(m)
	bins.list.corrected = lapply(bins.list, mean.correction)
	
	bins.sum.corrected = Reduce("+", bins.list.corrected)
	
	bins.expected = bins.sum.corrected / n
	
	bins.chi.score = bins.list.corrected
	for (i in 1:n) {
		bins.chi.score[[i]] = (bins.expected - bins.chi.score[[i]])^2 / bins.expected
	}
	
	return(Reduce("+", bins.chi.score))
}

#This functions applies the filter for overdispersed bins.
correct.bins.list = function(bins.list, chi.sum.bins, chi.dir, strand, sample.name) {
	degrees.of.freedom = n.best.control.samples # sample + controls - 1
	
	# Convert chi squares to a normal distribution score 
	chi.sum.bins.normalized = (chi.sum.bins - degrees.of.freedom) / (sqrt( 2 * degrees.of.freedom))

	# Derive correction factor for the read counts
	chi.sum.bins.correction.factor = chi.sum.bins / degrees.of.freedom
	


	n = length(bins.list)

	# For each sample, correct bins with high chi-square value
	for (i in 1:n) {
		m = bins.list[[i]]
		index = which(3.5 < chi.sum.bins.normalized) # variation between these bins is large
		m[index] = m[index] / chi.sum.bins.correction.factor[index]
		bins.list[[i]] = m
	}

	return(bins.list) # return corrected bin.list
}