#!/usr/bin/env Rscript

library(stringr)
library(Rsamtools)

debug = F
if (!debug) {
	
	#Read script arguments
	cargs <- commandArgs(TRUE)
	args=NULL
	if(length(cargs)>0){
		flags = grep("^--.*",cargs)
		values = (1:length(cargs))[-flags]
		args[values-1] = cargs[values]
		if(length(args)<tail(flags,1)){
			args[tail(flags,1)] = NA
			}
			names(args)[flags]=cargs[flags]
			}
		
		## parse command line arguments:
		bamfile = args['--bam']
		chromosomes = if (is.element('--chromosome', names(args))) str_split(args['--chromosome'],',')[[1]] else 1:22
		interval_list = args['--interval_list']
		csvfile = args['--csv']
		pdffile = args['--pdf']
		Rcovlist = args['--Rcovlist']
} else {
	bamfile = 'old1.bam'
	chromosomes = "Y"
	csvfile = 'TESTold1.chrmY.csv'
	pdffile = 'old1.chrmY.pdf'
	interval_list = 'SureSelect_All_Exon_50MB_exons_hg19_human_g1k_v37.interval_list'
}

# read interval list
if (!exists("interval")) {
	interval = read.csv(interval_list, sep="\t", header=F, comment.char="@")
	interval = data.frame("chromosome"=interval$V1, "start"=interval$V2, "end"=interval$V3, row.names=NULL)
}

cov.lst.chr = NULL # this list holds coverage for all chromosomes
for (chromosome.i in 1:length(chromosomes)) {
	chromosome = chromosomes[chromosome.i]
	cat("Chromosome", chromosome)
	
	chr.rows = which(interval$chromosome == chromosome)
	
	cov.lst = list() # coverage of all bases in all intervals (appended after each other) in this chromosome

	## get reads with start position on this chromosome within interval
	which <- GRanges(chromosome, IRanges(interval$start[chr.rows], interval$end[chr.rows]))
	what <- c("pos", "qwidth")
	param <- ScanBamParam(which = which, what = what, flag = scanBamFlag(isDuplicate = FALSE))
	bam <- scanBam(bamfile, index=bamfile, param=param)

	for (int.i in 1:length(chr.rows)) {# for each interval
	int = chr.rows[int.i]
	# create coverage.interval vector (0 init)
	S = interval$start[int]
	E = interval$end[int]
	cov.interval = rep(0, E - S + 1)

#thisbam = bam[[int.i]]
#thispos = thisbam$pos
#thisqwidth = thisbam$qwidth
	attach(bam[[int.i]])
	n.reads = length(pos)
	if (0 < n.reads) for (read.i in 1:n.reads) {# for each read in this interval
		p = pos[read.i]
		if (!is.na(p)) { # update local coverage
			read.width = qwidth[read.i]
			index = which(S:E %in% p:(p + read.width - 1))
			cov.interval[index] = cov.interval[index] + 1
			}
			}
		detach()

		# update global coverage vector
		cov.lst[[ length(cov.lst) + 1 ]] = cov.interval
		}
	
	# update global coverage vector per chromosome
	cov.lst.chr[[chromosome]] = cov.lst
}

# write cov.lst.chr as R object to file
dput(cov.lst.chr, file=Rcovlist)

# calculate coverage statistics for whole chromosome

full.cov.vec = unlist(cov.lst.chr)
nbases       = length(full.cov.vec)

cov.matrix = matrix(NA, nc = 2)

current.cum.coverage = 0
for (this.coverage in 0:max(full.cov.vec)) {
	n = length(which(full.cov.vec == this.coverage))
	vec = c(this.coverage, n / nbases)
	cov.matrix = rbind(cov.matrix, vec)
}

# strip first line (with NA's)
cov.matrix = cov.matrix[-1,]

# add cummulative coverage column
cov.matrix = cbind(cov.matrix, cumsum = rev(cumsum(rev(cov.matrix[,2]))))

# add column names:
colnames(cov.matrix) = c("coverage", "count", "cumulative_count")

# write table
write.table(cov.matrix, file=csvfile, sep='\t', row.names=F)

# make plot
pdf(pdffile)
par(mai = c(1.36, 1.5, 1.093333, 0.56))
plot(cov.matrix[,1], cov.matrix[,3], xlim = c(0,80), t='l', lwd = 3, axes = F, xlab='',ylab='', ylim=c(0,1))
axis(1, lwd=0, cex.axis=3, line=1)
axis(1, lwd=2, labels=F)
axis(2, cex.axis=3,las=2)
axis(2, lwd=2, labels=F)
title(xlab='coverage >= x', ylab='fraction', main='', cex.lab=3, line=5)#bamfile)

plotdot = function(coverage) {
	x = coverage
	y = cov.matrix[coverage + 1, 3]
	points(x, y, pch=19, cex=1.3, col='gray20')
	text(x,   y, paste('', round(y, 2)), pos = 4, cex=1.3, col='gray20')
}

plotdot(1:3*10)

# plot lines to the point (20, .80)
y20 = cov.matrix[20 + 1, 3]
lines(c(20, 20), c(-.1, y20)) # vertical line
lines(c(-10, 20), c(y20, y20)) # horizontal line

#points(20, .80, pch=19, cex=1, col='darkgreen')
#legend('topright', col=c('darkgreen'), pch=19, legend=c('fraction 0.8 with coverage >= 20x'), box.lty=0, cex=1, bg = 'gray90')
dev.off()