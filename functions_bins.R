get.bins = function(info, col.names. = col.names, col.type. = col.type, nselect = 1000, nbins = 300) {
    stopifnot(3 < nbins) 
    
    nselect = min(nselect, length(info))
    index.random = sample(length(info))[1:nselect]

    # NB! This could have been implemented better. Practically it would not have been necessary to use the function 'info.as.matrix'.
    # Storing the smallest and biggest value, for each column, in the sample (of size nselect) would have been enough to generate the bins.
    mat = info.as.matrix(info[index.random], col.names., col.type.)
 
    is.bool = which(col.type. == F)

    bin = count = matrix(0, nr = length(col.names.), nc = nbins)
    rownames(bin) = rownames(count) = col.names. 
    for (i in 1:nrow(bin)) {
        if (is.element(i, is.bool)) {
            bin[i, ] = c(F, T, rep(Inf, nbins - 2))
        } else {
            # count[i, j] should contain the number of values v_i where j = max{j}{ v_i <= bin[i, j] }
            values      = as.numeric(mat[,i])
            if (all(is.na(values))) {
                bin[i,]    = rep(NA, nbins)
            } else {
                this.bin  = seq(min(values, na.rm=T), max(values, na.rm=T), length = nbins)
                this.bin  = c(this.bin[-1], Inf)
                bin[i,]   = this.bin            
            }
        }
    }

    NA.count = rep(0, nrow(bin))

    list(bin = bin, count = count, NA.count = NA.count)
}

get.counts = function(info, bin.count) {
    bin    	 = bin.count$bin
    count  	 = bin.count$count
	NA.count = bin.count$NA.count

    if (0 < length(info)) for (i in 1:length(info)) {
        if (i %% round(length(info) / 20) == 0) print(paste('Processed ', round(100 * i / length(info)), '%', sep = ''))
        this.line   = info[[i]] # data from file

        if (0 < length(this.line)) for (j in 1:length(this.line)) {
            # parse the info for this line
            item = strsplit(this.line[j], '=')[[1]]
    
            # determine and fill its corresponding column
            index = which(col.names == item[1])

           # only do st for known fields
            if (0 < length(index)) {
 
                bool = length(item) == 1
                
                if (bool) {
                    count[index, 2] = count[index, 2] + 1
                } else {
					# in case of multiple values: use median
                    #value = as.numeric(item[2])
					value = median(as.numeric(strsplit(item[2],',')[[1]]))

					# report values that are missing and also report values that are not numeric
					if (is.na(value)) {
						NA.count[index] = NA.count[index] + 1
						print(paste("NA found on chromosome ", col1[i], ", position ", col2[i], sep=''))
					} else {
                    	k = which.max(value <= bin[index,])
                    	count[index, k] = count[index, k] + 1
					}
                }
            }
        }
    }
    
    # Booleans are only in the info if they are TRUE. Therefore, we thus 'miss' the FALSE ones.
    # Total minus #TRUE -> #FALSE
    for (i in 1:nrow(bin)) {
        if (!any(is.na(bin[i,1:3]))) if (all(bin[i, 1:3] == c(F, T, Inf))) { # we deal with a boolean
            # Fill in the number of FALSE ones
            if (0 < count[i, 2]) { # only if we found at least one TRUE boolean
                count[i, 1] = length(info) - count[i, 2]
            }
        }
    }
    
    
    list(bin = bin, count = count, NA.count = NA.count)
}

info.as.matrix = function(info, col.names, col.type) {

    mat = matrix(NA, nr = 0, nc = length(col.names), dimnames=list(NULL, col.names))
    
    if (0 < length(info)) for (i in 1:length(info)) {
        input.vec   = info[[i]] # data from file
        output.vec  = col.type
        if (0 < length(input.vec)) for (j in 1:length(input.vec)) {
            # parse the info for this line
            item = strsplit(input.vec[j], '=')[[1]]
    
            # is the type of this item boolean?
            bool = length(item) == 1
           
            # only do st for known fields
            if (is.element(item[1], col.names)) {
                # determine and fill its corresponding column
                if (bool) {
                    index = which(col.names == item)
                    output.vec[index] = T
                } else {
                    # item has form: <key>=<data>[,data]
                    index = which(col.names == item[1])
					
					# in case of multiple values, use the median
                    output.vec[index] = median(as.numeric(strsplit(item[2],',')[[1]]))
                }
            }
        }
        mat = rbind(mat, output.vec)
    }

    mat
}























#
