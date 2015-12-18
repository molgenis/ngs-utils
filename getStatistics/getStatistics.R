#!/usr/bin/env Rscript
# author: mdijkstra
# modified: 20140130
library(stringr)

debug = F

# Get the location where this .R script was installed.
getMyLocation = function()
{
	thisfile = NULL
	# This file may be 'sourced'
	for (i in -(1:sys.nframe())) {
		if (identical(sys.function(i), base::source)) thisfile = (normalizePath(sys.frame(i)$ofile))
	}
	
	if (!is.null(thisfile)) return(dirname(thisfile))
	
	# But it may also be called from the command line
	cmdArgs <- commandArgs(trailingOnly = FALSE)
	cmdArgsTrailing <- commandArgs(trailingOnly = TRUE)
	cmdArgs <- cmdArgs[seq.int(from=1, length.out=length(cmdArgs) - length(cmdArgsTrailing))]
	res <- gsub("^(?:--file=(.*)|.*)$", "\\1", cmdArgs)
	
	# If multiple --file arguments are given, R uses the last one
	res <- tail(res[res != ""], 1)
	if (length(res) > 0) return(dirname(res))
	
	# Both are not the case. Maybe we are in an R GUI?
	return(NULL)
}

if (!is.null(getMyLocation())) {
	my.location=getMyLocation()
	source(str_c(my.location,.Platform$file.sep,'getStatisticsReadCommandLineArgs.R'))
	source(str_c(my.location,.Platform$file.sep,'getStatisticsFunctions.R'))
} else {
	source('getStatisticsReadCommandLineArgs.R')
	source('getStatisticsFunctions.R')
}

hsmetrics.col.selection = c('BAIT_SET', 'GENOME_SIZE', 'BAIT_TERRITORY', 'TARGET_TERRITORY', 'TOTAL_READS', 'PCT_PF_UQ_READS_ALIGNED', 'PF_UQ_BASES_ALIGNED', 'ON_BAIT_BASES', 'NEAR_BAIT_BASES', 'OFF_BAIT_BASES', 'MEAN_BAIT_COVERAGE', 'MEAN_TARGET_COVERAGE', 'ZERO_CVG_TARGETS_PCT', paste('PCT_TARGET_BASES_', c('2X', '10X', '20X', '30X'), sep=''), 'PCT_USABLE_BASES_ON_TARGET')
almetrics.col.selection = c('MEAN_READ_LENGTH', 'STRAND_BALANCE')
insert.col.selection = c('MEDIAN_INSERT_SIZE', 'MEAN_INSERT_SIZE', 'STANDARD_DEVIATION')
dedup.col.selection = c('READ_PAIR_DUPLICATES', 'PERCENT_DUPLICATION')
conc.col.selection = c('nSNPs', 'Overall_concordance')

hsmat = NULL
for (fn in hsmetrics.files)
{
	thismat = read.csv(fn, sep='\t', as.is=T, comment.char='#')[, hsmetrics.col.selection]
	hsmat = rbind(hsmat, thismat)
}

almat = NULL
for (fn in almetrics.files)
{
	thismat = read.csv(fn, sep='\t', as.is=T, comment.char='#')[3, almetrics.col.selection]
	almat = rbind(almat, thismat)
}

inmat = NULL
for (fn in insertmetrics.files)
{
	if (!file.exists(fn)) {
		inmat = rbind(inmat, matrix(rep(NA, length(insert.col.selection)), dimnames=list(NULL,insert.col.selection), nrow=1))
	} else {
		inmat = rbind(inmat, retrieve.table(fn, 6))
	}
}
inmat = inmat[, insert.col.selection, drop=F]

ddmat = NULL
for (fn in dedupmetrics.files)
{
	ddmat = rbind(ddmat, retrieve.table(fn, 0))
}
# According to Pieter vd Vlies, we should divide these number by two
ddmat[, 'PERCENT_DUPLICATION'] = as.numeric(ddmat[, 'PERCENT_DUPLICATION']) / 2
ddmat = ddmat[, dedup.col.selection]

concmat = NULL
for (fn in concordance.files)
{
	concmat = rbind(concmat, read.csv(fn))
}
concmat = concmat[, conc.col.selection]
concmat[, 'Overall_concordance'] = concmat[, 'Overall_concordance'] / 100

mat = cbind(samples,hsmat,almat,inmat)
if (0 < nrow(concmat)) mat = cbind(mat, concmat)

#mat = (cbind(samples, hsmat, almat, inmat, ddmat, concmat))
#mat = cbind(samples, hsmat, almat, inmat, concmat)

# Add 'derived columns'
MbMapped = mat[, 'MEAN_TARGET_COVERAGE'] * mat[, "TARGET_TERRITORY"] / 1e6
mat = add.col.after(mat, "PCT_PF_UQ_READS_ALIGNED", "On target bases (Mb)", MbMapped)

totalbases = mat[, "TOTAL_READS"] * as.numeric(mat[, "MEAN_READ_LENGTH"])
fractionBpOnBait = mat[, "ON_BAIT_BASES"] / totalbases

mat = add.col.after(mat, "MEAN_TARGET_COVERAGE", "Fraction bp on bait", fractionBpOnBait)
fractionBpNearBait = mat[, "NEAR_BAIT_BASES"] / totalbases

mat = add.col.after(mat, "Fraction bp on bait", "Fraction bp near bait", fractionBpNearBait)
fractionBpOffBait = mat[, "OFF_BAIT_BASES"] / totalbases

mat = add.col.after(mat, "Fraction bp near bait", "Fraction bp off bait", fractionBpOffBait)
fractionBpNotAligned = 1 - mat[, "PCT_PF_UQ_READS_ALIGNED"]

mat = add.col.after(mat, "Fraction bp off bait", "Fraction bp not aligned", fractionBpNotAligned)
captureSpecificity = 1 - mat[, "Fraction bp off bait"]

mat = add.col.after(mat, "Fraction bp not aligned", "Capture specificity", captureSpecificity)
fractionOnNearTarget = fractionBpOnBait + fractionBpNearBait

# Convert to Mb
mat[,  'PF_UQ_BASES_ALIGNED'] = mat[,  'PF_UQ_BASES_ALIGNED'] / 1e6


# remove columns
mat = mat[, -which(colnames(mat) %in% c("ON_BAIT_BASES", "NEAR_BAIT_BASES", "OFF_BAIT_BASES", "PCT_PF_UQ_READS_ALIGNED", "READ_PAIR_DUPLICATES", "ZERO_CVG_TARGETS_PCT"))]

# round on 2 digits
if (!precise) for (i in 3:ncol(mat)) {
	mat[,i] = round(as.numeric(as.vector(mat[,i])), 2)
}


# add nice column headers
mat = substitute.col.names(as.matrix(mat))

# pivote table
tmat = t(mat)
colnames(tmat) = NULL

# write table
write.table(tmat[, -ncol(tmat)], csvout, quote=F, col.names=F, sep=',')

# remove second line and add it to 'latex description' later
baitSet = tmat[ 2, ]
tmat 	= tmat[-2, ]

# convert table to latex
# as a _side effect_ 'lateXDescription' is also made global
# create multiple tables with 4 samples per table if there are > 4 samples

saveAsLatex = F
if (saveAsLatex)
{
	latexString = NULL
	ncols = 3
	for (i in 0 : ((ncol(tmat)-2) %/% ncols))
	{
		index = unique(pmin(i * ncols + 1 : ncols, ncol(tmat) - 1))
		index = c(index, ncol(tmat))
		latexString = str_c(latexString, mat2Latex(tmat[, index]))	
	}


	#latexDescription = str_c(latexDescription,  "\\\\ \n ",  "\\\\ \n ", baitSet[length(baitSet)], ": ", unique(baitSet[1:(length(baitSet)-1)]), "\\\\ \n ")

	# write statistics
	write(latexString, tableout)

	# write description of statistics
	write(latexDescription, descriptionout)
}

# write bait set
write(unique(baitSet[1:(length(baitSet)-1)]), baitsetout)

# dedup metrics are per flowcell lane and not per sample, so write them separately to a file for seperate display in QCReport
# write dedupmetrics
# add nice column headers
ddmat = substitute.col.names(as.matrix(ddmat))

# pivote table
tddmat = t(ddmat)
colnames(tddmat) = NULL
write(tddmat, qcdedupmetricsout)





