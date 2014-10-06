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

BestPredictor = function(frac.reads.chr.trisomy.observed, chr.frac, chr.selected, i.predictor)
{
	# Define current sub set of predictors and 'candidate best predictors'
	chr.candidates = colnames(chr.frac)
	if (!is.null(chr.selected))
	{
		if (1 == i.predictor)
		{
			chr.frac.currect.subset = NULL
		} else
		{
			chr.frac.currect.subset	= chr.frac[, tail(chr.selected, i.predictor - 1)]
		}

		chr.candidates = chr.candidates[which(!(chr.candidates %in% chr.selected))]
	} else {
		chr.frac.currect.subset	= NULL
	}

	# Determine 'adjusted R-squared' for each of candidate predictors
	adj.r.squares = NULL
	for (i in 1:length(chr.candidates))
	{
		# Add candidate predictor to current subset
		chr.frac.test = cbind(chr.frac.currect.subset, chr.frac[, chr.candidates[i]])
		adj.r.squares[i] = summary(lm(frac.reads.chr.trisomy.observed ~ chr.frac.test))$adj.r.square
	}
	
	# Return name of 'best' predictor
	return(chr.candidates[which.max(adj.r.squares)])
}

# Cumulative standard normal distribution
CumulativeNormalDistribution <- function(x)
{
        pnorm(x, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)
}

# Function risk Computes chance (value 0..1) on a trisomy, P(H | E)
# Where:
# H is is the hypothesis trisomy
# E is the observation of the mother's DNA
# P(H) is the a priori chance on a trisomy
# P(H | E) is the chance on a trisomy, given the extra information from observation E
#
# Bayes rule:  P(H | E) = P(E | H) * P(H) / P(E)
# Where:
# P(E | H) is the likelihood; the chance of E, given H.
# P(E) is the chance of E, which is equal to
# the chance of E, given H correct plus the chance of E, given H incorrect:
# P(E) = P( E | H ) * P (H) + P( E | ~H ) * (1 - P(H))
# See also: http://www.answers.com/topic/bayesian-inference
Risk <- function(lower, upper, cv, z.observed, a.priori)
{
        p = NULL # the risk
        
        # Validate input
        if (upper < lower | a.priori < 0 | 1 < a.priori | is.na(a.priori)) p = NA
        for (test in c(lower, upper, cv)) if (test < 0 | 100 < test) p = NA
        
        if (is.null(p))
        {
                lower = 0.5 * lower / cv
                upper = 0.5 * upper / cv

                if (upper - lower < 0.002)
                {       # This is done to handle the case that the upper and lower limit are (almost) identical
                        lower = lower - 0.001
                        upper = upper + 0.001
                }

                interval = upper - lower

                p = (CumulativeNormalDistribution(z.observed - upper) - CumulativeNormalDistribution(z.observed - lower)) / interval;
                p = p * a.priori / (p * a.priori + (1 - a.priori) * dnorm(z.observed))
        }
        
        return(p)
}

#
## Functions to construct output tables
#
# 
GetModelInfoTable = function(model.list)
{
	results = NULL
	for (i.model in 1:n.models)
	{
		m = model.list[i.model]

		# Custom round function
		Round = function(x) round(x, round.n.decimals)
		
		# Names of m$coefs are names of predictors (e.g. 1F or 3R)
		result = c(Round(m$coefs), Round(unlist(summary(m)[c("r.squared", "adj.r.squared", "fstatistic", "sigma")])))
	
		# Add new row with 
		results = rbind(results, result)
	}

	rownames(results) = paste("Model", 1:n.models)

	return(as.data.frame(results))
}

GetRiskInfoTable = function(cv.model, normally.distributed.z.scores, z.score.sample, risk, risk.result)
{
	results = cbind(CV = cv.model, nd = normally.distributed.z.scores, z = z.score.sample, PPV = risk, ppvEventual = c(risk.result, rep("", length(risk) - 1)))
	colnames(results)[2] = "Normally distributed"
	colnames(results)[3] = "Z-score"
	colnames(results)[5] = "Median PPV"
	
	return(as.data.frame(results))
}