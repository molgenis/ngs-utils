#!/usr/bin/env Rscript
DOC = "For each of the chromosomes 1-22 and X and Y, this script counts the reads from a .bed input file, per bin of size 50,000 bp. It produces a matrix chr vs. bins, which is saved as a .tsv-file. In addition it produces a .PDF that visualizes the matrix. Mitochondrial DNA and contigs are ignored."

# Define constants
bin.size  = 5e4				# Fix bin size at 50,000 base pairs
chr.focus	= c(1:22, "X", "Y")	# Apply analysis only to chromosomes 1-22, X, Y

# Retrieve command line parameters
suppressPackageStartupMessages(library("argparser"))

parser <- arg.parser(DOC, name="bins")

print(parser)

parser <- add.argument(parser, "--input",  help = "Bed file with reads.")
parser <- add.argument(parser, "-o",  help = ".tsv file with number of reads per bin.")
parser <- add.argument(parser, "-p", 			help = ".pdf file that visualizes distribution of reads accross bins (log-scale).")
parser <- add.argument(parser, "-v", 	help = "Shows with which parameters this script is called [default].")
parser <- add.argument(parser, "-q", 		help = "Print little or no output.")

args <- parse.args(parser, argv = commandArgs(trailingOnly = TRUE))

# Proceed only if cmnd-line parameters correct
#if (args$verbose) {
  #write("You have used the following arguments:", stdout())
  #write(paste("\t--input:  ", args$input), stdout())
 # write(paste("\t--output: ", args$output), stdout())
 # write(paste("\t--pdf:    ", args$pdf), stdout())
  #write("\nStart with binning... (may take minutes)\n", stdout())
#}

# Sets the column names for the sample for easy acces later
colnames = c("chrom", "startpos", "stoppos", "readname", "score", "strand")

# Reads the sample reads
sample = (read.delim(args$input, header = FALSE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, col.names = colnames))

# Gets the chromosomes from all the reads from the sample as a vector. This is used to assign the read to the correct chromosome
chr = as.vector(sample$chrom)

# Start position of read is used to determine the bin
pos = sample$startpos


##Added to rename chromosome names; add chr so output is recognized. Added 20150702 LJ##
# Naming of chromosomes may be 1..22,X,Y or chr1..chr22,chrX,chrY
# Make chr.focus 
if ("chr" == substr(chr[1],1,3))
{
	chr.focus = paste("chr", chr.focus, sep = "")
}
###


# Removes reads which are not in chr.focus (for chr.focus, see top of this script)
index.remove = which(! chr %in% chr.focus)
chr = chr[-index.remove]
pos = pos[-index.remove]

# Gets the maximium start position of the reads, and determine the number of bins
max.pos	= max(pos)
n.bins	= max.pos %/% bin.size + 1

# Makes an empty matrix to store the bin counts in 
bin = matrix(0, nrow = length(unique(chr)), ncol = n.bins, dimnames = list(unique(chr), NULL))

# Bins the reads per chromosome. Method used is integer division by binsize + 1, for instance a read with startposition 125,000 divided by 50,000 yields 2, + 1 is 3. The read belongs in bin 3.

# Fill the bins
for (i in 1:length(pos))
{
  i.bin				= pos[i] %/% bin.size + 1	# determine bin
  bin[chr[i], i.bin]	= bin[chr[i], i.bin] + 1	# increase bin with 1
}

# Save binned data as .tsv
write.table(bin, args$o, quote = FALSE, sep ="\t", row.names = TRUE)

# Plot matrix and save as PDF
pdf(args$p)
image(log(bin), axes = F)
axis(1,at=seq(0,1,l=length(chr.focus)),labels=chr.focus)
dev.off()

# Return with exit code 0
quit(save = "no", status = 0)
