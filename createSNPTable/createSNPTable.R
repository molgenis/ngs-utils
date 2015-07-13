#!/usr/bin/env Rscript
# author: mdijkstra
# modified: 20140130

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
	source(str_c(my.location,.Platform$file.sep,'readCommandLineArgs.R'))
	source(str_c(my.location,.Platform$file.sep,'createSNPTableFunctions.R'))
} else {
	source('readCommandLineArgs.R')
	source('createSNPTableFunctions.R')
}

# Is this lib stil needed?
library(stringr)

# type table
write(mat2Latex(addHeader(readTables(type), sample),'Functional type'), typetableout)

# class table
write(mat2Latex(addHeader(readTables(class), sample),'Functional class'), classtableout) 

# impact table
write(mat2Latex(addHeader(readTables(impact), sample),'Functional impact'), impacttableout)