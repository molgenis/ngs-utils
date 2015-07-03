# Create matrix chr.frac, where chr.frac[sample, chrom] is the ratio '#reads on chrom' / '#reads in sample', where chrom is 1F, 2F, ... nF, 1R, 2R, ... or nR
ChromosomalFractionPerSample = function(control.file.names, bins.forward, bins.reverse)
{
	chr.frac = NULL
	for (i in 1:length(bins.forward))
	{
		# For sample 'i', get fwd and rev bins[chrom, bin.i] (1 ... bin.i ... number of bins)
		bins.fwd.i = bins.forward[[i]]
		bins.rev.i = bins.reverse[[i]]
		
		# Determine fraction of reads on each chromosome-direction
		chr.frac.fwd.i = rowSums(bins.fwd.i) / sum(bins.fwd.i)
		chr.frac.rev.i = rowSums(bins.rev.i) / sum(bins.rev.i)
		
		# Add these fractions as column
		chr.frac = rbind(chr.frac, c(chr.frac.fwd.i, chr.frac.rev.i))
	}

	# Set row and column names
	rownames(chr.frac) = sub("^([^.]*).*", "\\1", control.file.names) # sample names
	colnames(chr.frac) = c(paste("Chr", chromosomes.background, "F", sep =""), paste("Chr", chromosomes.background, "R", sep ="")) # 1F, 2F, ... nF, 1R, 2R, ... nR

	chr.frac
}

BestPredictor = function(n.reads.chr.trisomy, chr.frac, chr.selected)
{
	# Define current sub set of predictors and 'candidate best predictors'
	chr.candidates = colnames(chr.frac)
	if (!is.null(chr.selected)) {
		chr.frac.currect.subset	= chr.frac[, chr.selected]
		chr.candidates			= chr.candidates[which(!(chr.candidates %in% chr.selected))]
	} else {
		chr.frac.currect.subset	= NULL
	}

	# Determine 'adjusted R-squared' for each of candidate predictors
	adj.r.squares = NULL
	for (i in 1:length(chr.candidates))
	{
		# Add candidate predictor to current subset
		chr.frac.test = cbind(chr.frac.currect.subset, chr.frac[, chr.candidates[i]])
		adj.r.squares[i] = summary(lm(n.reads.chr.trisomy ~ chr.frac.test))$adj.r.square
	}
	
	# Return name of 'best' predictor
	return(chr.candidates[which.max(adj.r.squares)])
}
