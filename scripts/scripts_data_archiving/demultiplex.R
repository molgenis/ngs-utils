#!/usr/bin/env Rscript

#
##
### Setup environment
##
#
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(logging))
logging::basicConfig()

#
##
### Custom functions
##
#
usage <- function() {
    cat("Example usage:\n
    ./demultiplex.R --bcs 'AGAGAT,TAATTT,TCAGTT,TGACTT' \\
                     --mpr1 readOneFromPair.fq.gz \\
                     --mpr2 readTwoFromPair.fq.gz \\
                     --dmr1 'readOneFromPair_AGAGAT.fq.gz,readOneFromPair_TAATTT.fq.gz,readOneFromPair_TCAGTT.fq.gz,readOneFromPair_TGACTT.fq.gz' \\
                     --dmr2 'readTwoFromPair_AGAGAT.fq.gz,readTwoFromPair_TAATTT.fq.gz,readTwoFromPair_TCAGTT.fq.gz,readTwoFromPair_TGACTT.fq.gz' \\
                     --ukr1 readOneFromPair_UNKNOWN.fq.gz \\
                     --ukr2 readTwoFromPair_UNKNOWN.fq.gz \\
                     --ll   WARNING \n\n")
    
    cat("Explanation of flags:\n")
    cat("    --bcs:     List of nucleotide barcodes.\n")
    cat("    --mms:     Maximum number of allowed mismatches in the barcode (defaults to auto detection of max mismatches allowed for unambiguous detection of barcodes).\n")
    cat("    --force:   Force maximum number of allowed mismatches in the barcode.\n")
    cat("               Overrules check for max mismatches allowed to prevent unambiguous detection of barcodes.\n")
    cat("               Use only for reads with lot's of bad quality base calls in the barcodes and when you know what you are doing!\n")
    cat("    --mpr1:    Required input FastQ file containing MultiPlexed sequence Reads.\n")
    cat("               These are either the only reads in case of single read sequencing or all reads from one end of a paired end or mate pair library.\n")
    cat("               File may be compressed with gzip.\n")
    cat("    --mpr2:    Optional input FastQ file containing MultiPlexed sequence Reads.\n")
    cat("               These are all reads from the other end of a paired end or mate pair library.\n")
    cat("               Note: --mpr[1|2] FastQ files maybe in supplied in gzipped format if the filename extension is *.fq.\n")
    cat("               File may be compressed with gzip.\n")
    cat("    --nocheck: Disable read ID consistency check for paired end or mate pair libraries.\n")
    cat("               Do not throw an error when the sequence read IDs for both ends of a sequence pair are not identical.\n")
    cat("    --dmr1:    Required list of gzipped FastQ output file names for the DeMultiplexed sequence Reads.\n")
    cat("               These will contain either the only reads in case of single read sequencing or all reads from one end of a paired end or mate pair library.\n")
    cat("               Files will be compressed with gzip.\n")
    cat("    --dmr2:    Optional list of gzipped FastQ output file names for the DeMultiplexed sequence Reads.\n")
    cat("               These will contain all reads from the other end of a paired end or mate pair library.\n")
    cat("               Note: the order of the --dmr[1|2] output file names is important and shpould match the order of the barcodes supplied with --bcs.\n")
    cat("               Files will be compressed with gzip.\n")
    cat("    --ukr1:    Required gzipped FastQ output file name for the UnKnown sequence Reads in which none of the supplied barcodes was detected.\n")
    cat("               These are either the only UnKnown Reads in case of single read sequencing or all UnKnown Reads from one end of a paired end or mate pair library.\n")
    cat("               File will be compressed with gzip.\n")
    cat("    --ukr2:    Optional gzipped FastQ output file name for the UnKnown sequence Reads in which none of the supplied barcodes was detected.\n")
    cat("               These are all UnKnown Reads from the other end of a paired end or mate pair library.\n")
    cat("               File will be compressed with gzip.\n")
    cat("    --tm:      Threading mode: ST, MP (default).\n")
    cat("               ST: Single Thread. Decompression (in case the input was compressed), demultiplexing and compression of the output is all done in R on a single CPU/core.\n")
    cat("               MP: Multiple Pipes. In addition to a thread for demultiplexing in R, additional threads (one per file) are created.\n")
    cat("                   Decompressing input and compressing output is handled outside R by gzip and uncompressed data is streamed through a pipe resulting in additional threads.\n")
    cat("    --ll:      Log level (optional). One of FINEST, FINER, FINE, DEBUG, INFO (default), WARNING, ERROR or CRITICAL.\n\n")
    q()
}

#
# Create an object to handle IO.
#   
#   1. Either a classic "file handle" -> will work for both compressed and uncompressed input files.
#   2. Or a pipe to stream data from/to gzip for (de)compression in a separate thread -> works only with compressed data.
#
create.filehandle <- function(path, mode, threading_mode) {
    
    if (threading_mode == 'ST') {
        if (mode == 'r') {
            fh = file(path, mode)
        } else if (mode == 'w') {
            fh = gzfile(path, mode)
        } else {
            logging::levellog(loglevels[['FATAL']], 'Bad file access mode!')
            q()
        }
    } else if (threading_mode == 'MP') {
        if (mode == 'r') {
            fh = pipe(paste('gzip -cd', path, sep=' '), "r")
        } else if (mode == 'w') {
            fh = pipe(paste('gzip -c >', path, sep=' '), "w")
        } else {
            logging::levellog(loglevels[['FATAL']], 'Bad file access mode!')
            q()
        }
    } else {
        logging::levellog(loglevels[['FATAL']], 'Bad threading mode!')
        q()
    }
    
    return(fh)
    
}

compile.output_files <- function(mpr_path, read.n, barcodes, args, dmr_fhs, ukr_fhs, threading_mode) {
    
    logging::levellog(loglevels[['DEBUG']], 'Compiling paths and filehandles for output files...')
    #
    # Extract info from the file name/path with multiplexed reads (input) 
    # to create defaults for the output file names.
    #
    mpr_path.file = basename(mpr_path)[[1]]
    mpr_path.file.name = str_match(mpr_path.file, '^([^\\.]+)')[[1]][1]
    mpr_path.file.ext  = str_match(mpr_path.file, '(\\.[^\\.].*)$')[[1]][1]
    logging::levellog(loglevels[['FINE']], paste('MPR_', read.n, ' file = ', mpr_path.file, '; name = ', mpr_path.file.name, '; extension = ', mpr_path.file.ext, sep=''))
    
    #
    # Create output files for the demultiplexed reads in which we detected one of the barcodes.
    #  * One file per barcode.
    #
    dmr_list = character(0)
    dmr_flag = paste('--dmr', read.n, sep='')
    dmr_fhs[[read.n]] = list()
    barcodes.count = length(barcodes)
    dmr_and_ukr_extension = '.fq.gz'
    
    if (!is.na(args[dmr_flag])) {
        dmr_list = strsplit(args[dmr_flag], ',')[[1]]
    }
    
    for (i in 1:barcodes.count) {
        
        if (!is.null(dmr_list[i]) && !is.na(dmr_list[i])) { 
            #
            # Use the given output file names if supplied.
            #
            dmr_path = dmr_list[i]
            if (substr(dmr_path, nchar(dmr_path) - 5, nchar(dmr_path)) != dmr_and_ukr_extension) {
                logging::levellog(loglevels[['WARN']], paste('DMR_', read.n, ' file = ', dmr_path, ' does end with *', dmr_and_ukr_extension, ', which is the default for files created by this tool.', sep=''))
            }
        } else { 
            #
            # Generate output file names based on the names of the multiplexed FastQ input files supplemented with a barcode.
            #
            dmr_path  = paste(mpr_path.file.name, '_', barcodes[i], dmr_and_ukr_extension, sep='')
        }
        
        logging::levellog(loglevels[['FINE']], paste('DMR_', read.n, ' file for sample ', barcodes[i], ' = ', dmr_path, sep=''))
        dmr_fhs[[read.n]][[i]] = create.filehandle(dmr_path, 'w', threading_mode)
        
    }
    
    #
    # Create output files for reads in which we did not detect a barcode.
    #
    ukr_flag = paste('--ukr', read.n, sep='')
    if (!is.na(args[ukr_flag])) {
        # use the given file name
        ukr_path = args[ukr_flag]
        if (substr(ukr_path, nchar(ukr_path) - 5, nchar(ukr_path)) != dmr_and_ukr_extension) {
            logging::levellog(loglevels[['WARN']], paste('UKR_', read.n, ' file = ', ukr_path, ' does end with *', dmr_and_ukr_extension, ', which is the default for files created by this tool.', sep=''))
        }
    } else {
        # otherwise use standard file name
        ukr_path = paste(mpr_path.file.name, '_DISCARDED', dmr_and_ukr_extension, sep='')
    }
    logging::levellog(loglevels[['FINE']], paste('UKR_', read.n, ' file = ', ukr_path, sep=''))
    ukr_fhs[[read.n]] = create.filehandle(ukr_path, 'w', threading_mode)
    
    logging::levellog(loglevels[['DEBUG']], 'Created all filehandles for the output files.')
    
    fhs = list(dmr_fhs, ukr_fhs)
    return(fhs)
    
}

#
# Functions to process sequence barcodes.
#
verify.barcodes = function(barcodes.count, barcodes) {
    
    logging::levellog(loglevels[['DEBUG']], 'Checking if barcodes specified are sane...')
    #
    # Get the length of the first barcode.
    #
    first.barcode.length  = nchar(barcodes[1])

    #
    # Compare length of first barcode to that of the others and complain if they are not equal.
    #
    for (i in 2:barcodes.count) {
        other.barcode.length  = nchar(barcodes[i])
        if (first.barcode.length != other.barcode.length) {
            logging::levellog(loglevels[['FATAL']], 'Barcodes not of equal length!')
            q()
        }
    }
    
    #
    # Check if barcodes contain no illegal non-nucleotide-representing characters. 
    #
    for (i in 1:barcodes.count) {
        barcode = barcodes[i]
        illegal.characters = str_match(barcode, '([^ATCGU]+)')[[1]][1]
        if (!is.na(illegal.characters)) {
            logging::levellog(loglevels[['FATAL']], paste('Barcode # ', i, ' contains illegal characters! (Only nucleotides A, T, U, C and G are allowed.)', sep=''))
            logging::levellog(loglevels[['DEBUG']], paste('-----> Illegal character(s) found: ', illegal.characters, '.', sep=''))
            q()
        }
    }
    
    logging::levellog(loglevels[['DEBUG']], paste('All ', barcodes.count, ' are OK.', sep=''))
    
    return(first.barcode.length)
    
}

#
# Get sequence representing barcode from complete read.
#
get.barcode = function(read) {
    
    bc = strsplit(substr(read[2], 1, barcode.length), '')[[1]]
    
    mm = rep(0, barcodes.count)
    for (i in 1:barcodes.count) mm[i] = length(which(barcode.matrix[i,] != bc))
    
    index.match = which(mm == min(mm))
    if (length(index.match) == 1) if (mm[index.match] <= mismatches) return(index.match)
    
    return(NA)
    
}

get.max.mismatches = function(barcodes.count, barcode.length, barcode.matrix) {
    
    logging::levellog(loglevels[['DEBUG']], 'Comparing barcodes to each other to determine max similarity.')
    max.matches = 0
    
    for (i in 1:barcodes.count) {
        matches = rep(0, barcodes.count)
        barcode = barcode.matrix[i,]
        barcode.chars = strsplit(substr(barcode, 1, barcode.length), '')
        barcode.string = paste(barcode, sep='', collapse='')
        logging::levellog(loglevels[['FINE']], paste('Checking barcode # ', i, ': ', barcode.string, sep=''))
        for (j in 1:barcodes.count) {
            logging::levellog(loglevels[['FINE']], paste('--> Calculating similarty to barcode # ', j, sep=''))
            identical.positions = which(barcode.matrix[j,] == barcode.chars)
            logging::levellog(loglevels[['FINE']], paste('--> Found identical nucleotides on positions: ', paste(identical.positions, sep='', collapse=' '), sep=''))
            identical.positions.count = length(identical.positions)
            logging::levellog(loglevels[['FINE']], paste('--> Number of identical nucleotides: ', identical.positions.count, sep=''))
            matches[j] = identical.positions.count
        }
        #
        # Remove the barcode we are currently comparing to the others from the list
        # as that one will always be 100% identical to itself.
        #
        matches <- matches[-i]
        index.match = which(matches == max(matches))
        logging::levellog(loglevels[['FINE']], paste('--> Max matches: ', matches[index.match][[1]], sep=''))
        if (matches[index.match][[1]] > max.matches) {
            max.matches = matches[index.match][[1]]
            logging::levellog(loglevels[['FINER']], paste('--> Found new max matches: ', max.matches, sep=''))
        } else {
            logging::levellog(loglevels[['FINER']], paste('--> Max matches still the same: ', max.matches, sep=''))
        }
        
    }
    
    logging::levellog(loglevels[['DEBUG']], paste('Found max matches: ', max.matches, sep=''))
    mismatches = barcode.length - max.matches - 1
    logging::levellog(loglevels[['DEBUG']], paste('Setting max number of allowed mismatches in barcodes to: ', mismatches, sep=''))
    
    return(mismatches)
    
}

#
# Remove barcode.
#
remove.barcode = function(read) {
    #
    # Remove the barcode nucleotides from the sequence.
    #
    read[2] = substr(read[2], barcode.length + 1, nchar(read[2]))
    #
    # Remove quality scores for the nucleotides of the barcode too.
    #
    read[4] = substr(read[4], barcode.length + 1, nchar(read[4]))
    #
    # Return trimmed sequence read.
    #
    read
}

#
##
### Main.
##
#

#
# Read script arguments
#
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

arglist = c("--bcs", "--mms", "--force",  "--mpr1", "--mpr2", "--nocheck", "--dmr1", "--dmr2", "--ukr1", "--ukr2", "--tm", "--ll", NA)
#print(paste("!all(args %in% arglist)", !all(names(args) %in% arglist)))

#
# Handle arguments required to setup logging first.
#
if (is.element('--ll', names(args))) {
    log_level = args['--ll']
    log_level.position <- which(names(logging::loglevels) == log_level)
    if(length(log_level.position) < 1) {
        logging::levellog(loglevels[['WARN']], paste('Illegal log level ', log_level, ' specified.', sep=''))
        # Use default log level.
        log_level = 'INFO'
    }
    # Use the given log level.
} else {
    # Use default log level.
    log_level = 'INFO'
}

#
# Change the log lovel of both the root logger and it's default handler (STDOUT).
#
logging::setLevel(log_level)
logging::setLevel(log_level, logging::getHandler('basic.stdout'))
logging::levellog(loglevels[['INFO']], paste('Log level set to ', log_level, '.', sep=''))

#
# Check other arguments.
#
wrong.flags = length(args) == 0
if (!wrong.flags) wrong.flags = !all(names(args) %in% arglist)
if (!wrong.flags) wrong.flags = is.na(args['--bcs']) | is.na(args['--mpr1'])

if(wrong.flags) {
	if (!all(names(args) %in% arglist)) {
        logging::levellog(loglevels[['FATAL']], paste('Illegal parameter name or bad syntax for ', names(args)[!(names(args) %in% arglist)], '!', sep=''))
    }
    usage()
}

#
# Process other arguments.
#
mpr1_path	= args['--mpr1']
mpr2_path	= args['--mpr2']
nocheck		= is.element('--nocheck', names(args))
bcs			= args['--bcs']
barcodes	= strsplit(bcs, ',')[[1]]
barcodes.count = length(barcodes)
barcode.length = verify.barcodes(barcodes.count, barcodes)
barcode.matrix = matrix(unlist(strsplit(barcodes, '')), nc = barcode.length, byrow = T) # barcodes are below each other in a matrix (row number can be used to identify a barcode)
#
# Determine the amount of mismatches allowed when detecting barcodes in reads.
#
max.mismatches = get.max.mismatches(barcodes.count, barcode.length, barcode.matrix)
if (is.na(args['--mms'])) {
    #
    # Nothing specified by the user -> automatically determinine the max amount of mismatches allowed for detecting barcodes uniquely.
    #
    mismatches = max.mismatches
} else {
    if (args['--mms'] <= max.mismatches) {
        #
        # Set max mismatches to the value specified by the user, which is more stringent than strictly necessary for detecting barcodes uniquely.
        #
        mismatches = args['--mms']
    } else {
        logging::levellog(loglevels[['WARN']], 'Number of mismatches threshold specified is too lenient for the barcodes specified!')
        logging::levellog(loglevels[['WARN']], '--> This could result in ambiguous assignment of reads to barcodes/samples')
        if (is.element('--force', names(args))) {
            mismatches = args['--mms']
            logging::levellog(loglevels[['WARN']], paste('--> --force specified; Not overriding user input and allowing mismatches to be set to: ', mismatches, ', but you have been warned!', sep=''))
        } else {
            mismatches = max.mismatches
            logging::levellog(loglevels[['WARN']], paste('--> Overriding user input and setting mismatches to: ', mismatches, '.', sep=''))
        }
    }
}
logging::levellog(loglevels[['INFO']], paste('Maximum number of allowed mismatches in the barcodes: ', mismatches, '.', sep=''))

#
# Determine if we are dealing with Single Read, Paired End or Mate Pair sequencing.
#
seq_type_read_count = 1 # Single Read sequencing; default.
if (!is.na(args['--mpr2'])) {
    seq_type_read_count = 2 # Paired End or Mate Pair sequencing. 
}

#
# Determine Threading Mode
#
if (is.na(args['--tm'])) {
    #
    # Nothing specified by the user -> deafult to Single Thread (ST).
    #
    threading_mode = 'ST'
} else {
    threading_mode = args['--tm']
    if (threading_mode == 'ST' || threading_mode == 'MP') {
        logging::levellog(loglevels[['DEBUG']], paste('Threading mode set to: ', threading_mode, '.', sep=''))
    } else {
        logging::levellog(loglevels[['WARN']], paste('Illegal value specified for threading mode: ', threading_mode, '.', sep=''))
        threading_mode = 'ST'
        logging::levellog(loglevels[['WARN']], paste('Threading mode reset to default: ', threading_mode, '.', sep=''))
    }
}   
        
#
# Create filehandles.
#
# Option A: FileHandles (fh), commented out: easy to use, but single threaded:
#           decompressing input, demultiplexing and compressing output all on happen on the same single core.
# Option B: Using pipes, active: decompressing input and output in seperate threads for spead improvement.
#           Harder to use when submitting jobs to a PBS cluster as the optimal amount of cores to request 
#           depends on the amount of input & output files, which depends on sequencing technology and multiplexed samples.
#
dmr_fhs = list()
ukr_fhs = list()
mpr1_fh = create.filehandle(mpr1_path, 'r', threading_mode)
fhs = compile.output_files(mpr1_path, 1, barcodes, args, dmr_fhs, ukr_fhs, threading_mode)
if (seq_type_read_count == 2) {
    mpr2_fh = create.filehandle(mpr2_path, 'r', threading_mode)
    fhs = compile.output_files(mpr2_path, 2, barcodes, args, fhs[[1]], fhs[[2]], threading_mode)
}
dmr_fhs = fhs[[1]]
ukr_fhs = fhs[[2]]

#
# Variables for statistics.
#
n.total = n.bc.matched.zero = n.bc.matched.one = n.bc.matched.both.identical = n.bc.matched.both.different = 0
# how frequent is each of the barcodes found (index nr corresponds to rownumber in barcode.matrix):
bc.n.found = rep(0, barcodes.count)

#
# Process all reads.
#
ready = F
while (!ready) {
    
	read_1 = readLines(mpr1_fh, 4)
    
	if (length(read_1) == 0) {
		ready = T
        break
	}
    
	n.total = n.total + 1
    bc_1 = get.barcode(read_1) # Barcode number (corresponding to row number in barcode.matrix)

    if (seq_type_read_count == 1) {
       
        if (is.na(bc_1)) {
            
            # Write unidentified read.
            writeLines(read_1, ukr_fhs[[seq_type_read_count]])
            # Update statistics.
            n.bc.matched.zero = n.bc.matched.zero + 1
                    
        } else {
            
            # Remove the barcode.
            read_1 = remove.barcode(read_1)
            # Write identified read.
            writeLines(read_1, dmr_fhs[[seq_type_read_count]][[bc_1]])
            # Update statistics.
            n.bc.matched.one = n.bc.matched.one + 1
            bc.n.found[bc_1] = bc.n.found[bc_1] + 1
        }
        
    } else if (seq_type_read_count == 2) {
        
        read_2 = readLines(mpr2_fh, 4)
        bc_2   = get.barcode(read_2)
        
        if (!nocheck) {
            # Exit if the read IDs for a reads of a pair are not equal.
            # Works only for Illumina 1.7 read IDs:
            #stopifnot(substr(read_1[1], 1, nchar(read_1[1]) - 1) == substr(read_2[1], 1, nchar(read_2[1]) - 1))
            # Works with Illumina 1.7 and 1.8 read IDs:
            read_id_1 = strsplit(read_1[1], "[#[:space:]]")[[1]][1]
            read_id_2 = strsplit(read_2[1], "[#[:space:]]")[[1]][1]
            logging::levellog(loglevels[['DEBUG']], paste('Read1 ID: ', read_id_1 , sep=''))
            logging::levellog(loglevels[['DEBUG']], paste('Read2 ID: ', read_id_2 , sep=''))
            stopifnot(read_id_1 == read_id_2)
        }

        if (is.na(bc_1) && is.na(bc_2)) {
            
            # Write unidentified reads.
            writeLines(read_1, ukr_fhs[[1]])
            writeLines(read_2, ukr_fhs[[2]])
            # Update statistics.
            n.bc.matched.zero = n.bc.matched.zero + 1
            
        } else if (is.na(bc_1) || is.na(bc_2)) {
            
            if (is.na(bc_1)) bc_1 = bc_2
            # Remove the barcode.
            read_1 = remove.barcode(read_1)
            read_2 = remove.barcode(read_2)
            # Write identified read..
            writeLines(read_1, dmr_fhs[[1]][[bc_1]])
            writeLines(read_2, dmr_fhs[[2]][[bc_1]])
            # Update statistics
            n.bc.matched.one = n.bc.matched.one + 1
            bc.n.found[bc_1] = bc.n.found[bc_1] + 1
            
        } else {
            
            if (bc_1 != bc_2) {
                
                # Discard reads.
                writeLines(read_1, ukr_fhs[[1]])
                writeLines(read_2, ukr_fhs[[2]])
                # Update statistics
                n.bc.matched.both.different = n.bc.matched.both.different + 1

            } else {
                
                # Remove the barcode.
                read_1 = remove.barcode(read_1)
                read_2 = remove.barcode(read_2)
                # Write identified reads.
                writeLines(read_1, dmr_fhs[[1]][[bc_1]])    
                writeLines(read_2, dmr_fhs[[2]][[bc_1]])
                # Update statistics.
                n.bc.matched.both.identical = n.bc.matched.both.identical + 1
                bc.n.found[bc_1] = bc.n.found[bc_1] + 1
                
            }
        }
	}
    
    #logging::levellog(loglevels[['FINER']], paste('Processed:', n.total, 'reads/pairs.', sep=' '))
    
}

#
# Save statistics to log file.
#
nsmall = 2
numberwidth = floor(log(n.total,10))+1

if (seq_type_read_count == 1) {
    
    extraspacing = paste(rep(' ', numberwidth), collapse='')
    logging::levellog(loglevels[['INFO']], 'Identified reads:')
    for (i in 1:barcodes.count) {
        this.percentage = round(bc.n.found[i] / n.total * 100, 2)
        logging::levellog(loglevels[['INFO']], paste(':   * ', barcodes[i], ': ', format(bc.n.found[i], width = numberwidth), '  (', format(this.percentage, width = 5, nsmall = nsmall), '%)', sep=''))
    }
    logging::levellog(loglevels[['INFO']], paste(':   ======================================', paste(rep('=', numberwidth), collapse=''), sep=''))
    logging::levellog(loglevels[['INFO']], paste(':   * Total identified reads:               ', extraspacing, format(sum(bc.n.found), width = numberwidth), '  (', format(round(sum(bc.n.found) / n.total* 100, 2), width = 5, nsmall = nsmall), '%)', sep='')) 
    logging::levellog(loglevels[['INFO']], 'Unknown reads:')  
    logging::levellog(loglevels[['INFO']], paste(':   * Total unknown reads:                  ', extraspacing, format(n.bc.matched.zero, width = numberwidth), '  (', format(round(n.bc.matched.zero / n.total * 100, 2), width = 5, nsmall = nsmall), '%)', sep=''))
    logging::levellog(loglevels[['INFO']], paste(':                                          ', extraspacing, paste(rep('_', numberwidth + 2), collapse=''), '+', sep=''))    
    logging::levellog(loglevels[['INFO']], paste('Total number of reads:                      ', extraspacing, format(n.total, width = numberwidth), sep=''))
    
} else if (seq_type_read_count == 2) {
    
    extraspacing = paste(rep(' ', numberwidth), collapse='')
    logging::levellog(loglevels[['INFO']], 'Identified read pairs:')  
    for (i in 1:barcodes.count) {
        this.percentage = round(bc.n.found[i] / n.total * 100, 2)
        logging::levellog(loglevels[['INFO']], paste(':   * ', barcodes[i], ': ', format(bc.n.found[i], width = numberwidth), '  (', format(this.percentage, width = 5, nsmall = nsmall), '%)', sep=''))
    }
    logging::levellog(loglevels[['INFO']], paste(':   ===========================================', paste(rep('=', numberwidth), collapse=''), sep=''))
    logging::levellog(loglevels[['INFO']], paste(':   * Identical barcodes detected:   ', format(n.bc.matched.both.identical, width = numberwidth), '  (', format(round(n.bc.matched.both.identical / n.total * 100, 2), width = 5, nsmall = nsmall), '%)', sep=''))
    logging::levellog(loglevels[['INFO']], paste(':   * Only one barcode detected:     ', format(n.bc.matched.one, width = numberwidth), '  (', format(round(n.bc.matched.one / n.total * 100, 2), width = 5, nsmall = nsmall), '%)', sep=''))
    logging::levellog(loglevels[['INFO']], paste(':   * Total identified read pairs:               ', extraspacing, format(sum(bc.n.found), width = numberwidth), '  (', format(round(sum(bc.n.found) / n.total * 100, 2), width = 5, nsmall = nsmall), '%)', sep='')) 
    logging::levellog(loglevels[['INFO']], 'Unknown read pairs:')  
    logging::levellog(loglevels[['INFO']], paste(':   * Both barcodes not detected:    ', format(n.bc.matched.zero, width = numberwidth), '  (', format(round(n.bc.matched.zero / n.total * 100, 2), width = 5, nsmall = nsmall), '%)', sep=''))
    logging::levellog(loglevels[['INFO']], paste(':   * Conflicting barcodes detected: ', format(n.bc.matched.both.different, width = numberwidth), '  (', format(round(n.bc.matched.both.different / n.total * 100, 2), width = 5, nsmall = nsmall), '%)', sep=''))
    logging::levellog(loglevels[['INFO']], paste(':   * Total unknown read pairs:                  ', extraspacing, format(n.bc.matched.zero + n.bc.matched.both.different, width = numberwidth), '  (', format(round((n.bc.matched.zero + n.bc.matched.both.different) / n.total * 100, 2), width = 5, nsmall = nsmall), '%)', sep=''))
    logging::levellog(loglevels[['INFO']], paste(':                                               ', extraspacing, paste(rep('_', numberwidth + 2), collapse=''), '+', sep=''))    
    logging::levellog(loglevels[['INFO']], paste('Total number of read pairs:                      ', extraspacing, format(n.total, width = numberwidth), sep=''))
    
}

#
# Close file handles.
#
close(mpr1_fh)
if (seq_type_read_count == 2) {
    close(mpr2_fh)
}
for (i in 1:seq_type_read_count) {
    for (j in 1:barcodes.count) {
    	close(dmr_fhs[[i]][[j]])
    }
    close(ukr_fhs[[i]])
}

#
# We are done!
#
logging::levellog(loglevels[['INFO']], 'Finished!')

