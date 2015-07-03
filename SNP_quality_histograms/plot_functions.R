nice.label = function(x.seq) {

    step.size   = diff(range(x.seq)) / length(x.seq)
    ndigits     = - round( log(step.size, 10) ) + 1
    
    round(x.seq, ndigits)    
}

plot.histograms = function(bin.count, col.names.desc, col.type. = col.type) {
    par(mfrow = c(4,4))
    
    bin      = bin.count$bin
    count    = bin.count$count
	NA.count = bin.count$NA.count
    
    is.bool = which(col.type. == F)
    
    for (i in 1:nrow(count)) {
        x = bin[i,]
        y = count[i,]
        
        this.name = names(col.names.desc)[i]
        
        if (i %in% is.bool) {
            x = x[1:2]
            y = y[1:2]

            nfalse  = y[1]            
            ntrue   = y[2]
            
            ntrue.perc = round(100 * ntrue / (ntrue + nfalse))
    
            xargs   = c(paste('False (', 100 - ntrue.perc, '%)', sep=''), paste('True (', ntrue.perc, '%)', sep=''))
        } else {
            xargs = c(nice.label(x[-length(x)]), Inf)
            
            colnames.zero = c('Dels', 'HRun', 'MQ0')
            if (this.name %in% colnames.zero) {            
                nzeros  = y[1]
                y[1]    = 0
            }
            
        }

        bp.info = barplot(y, main = rownames(count)[i], col.lab='gray30', axes = F, names.arg = xargs)

        if (is.element(this.name, colnames.zero)) {
            text(tail(bp.info[,1],1), .9 * max(y), paste('zero: ', round(100 * nzeros / (nzeros + sum(y))), '%', sep=''), pos=2, col='darkgreen')
        }

		if (this.name == 'DP') {
            text(tail(bp.info[,1],1), .9 * max(y), paste('modus: ', round(x[which.max(y)]), sep=''), pos=2, col='darkblue')
		}
		
		if (0 < NA.count[i]) mtext(paste(NA.count[i], 'NA\'s found'), cex = .6, col='darkred')
        
        axis(2, las=2)
        title(xlab = col.names.desc[i], col='red')
    }
}

#plot.histograms(bin.count, col.names.desc)




















#
